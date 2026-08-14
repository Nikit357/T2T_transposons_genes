#!/usr/bin/env python
"""WP10 — split the over-wide Table 1 into Tables 1 and 2 (reviewer minor comment 4).

The reviewer's comment
----------------------
"Table 1 appears too wide... consider splitting into two tables, reducing column
widths, or moving less critical columns to supplementary materials."

The published Table 1 is 11 columns x 6 rows. This script reproduces it from the
stored data and emits the agreed split:

    Table 1. Enrichment of TE classes in gene TSS 10 kb neighbourhoods.  (5 cols)
        Class | TEs in TSS windows | TEs total | Observed OR | Observed/random OR

    Table 2. Statistical support for TE class enrichment.                (5 cols)
        Class | Adjusted Fisher p | Random OR (mean +- SD) | Empirical p | Adjusted empirical p

Plus `TableS_class_enrichment_full.csv`, which keeps every column including the
RAW Fisher p that the main tables drop — G3's statistics policy requires raw
p-values to be available in the supporting material so readers can apply their own
correction (WP9).

**Values do not change.** Under decision D1 the permutation background stays at
N = 500, so every OR, random OR, fold change and empirical p is exactly as
published. This is a reformat, not a regeneration — and the script asserts that by
checking each reproduced number against the values read out of the submitted
manuscript's Table 1 (`PUBLISHED_TABLE1`).

How the class-level numbers are derived
---------------------------------------
There is no stored class-level enrichment CSV — only the family- and subfamily-level
ones were saved — so the class rows are recomputed with the **exact contingency
table the frozen `TEs_mapped_on_TSS_analysis.ipynb` cell 55 used** (that cell is
commented "THIS IS TABLE 1"):

    N_Genome = 3,117,275,501       # T2T-CHM13v2.0 assembly length in bp
    N_TSS    =   272,233,268       # bp covered by the 10 kb TSS neighbourhoods
    K = TEs of this class genome-wide
    k = TEs of this class inside TSS neighbourhoods

    table = [[k,     N_TSS - k],
             [K - k, (N_Genome - N_TSS) - (K - k)]]
    fisher_exact(table, alternative="two-sided")

This mixes element counts with base-pair counts, which is unusual, but it is what
was published and D1 requires the values to stand, so it is reproduced verbatim
rather than corrected. The random OR uses the identical formula with `k` replaced by
each permutation's own class count.

Two things had to be got right, and both are verified rather than assumed:

  * `K` and `k` come from the **class column directly**, not from aggregating the
    44 curated families. The curated list omits some LTR and DNA families, so
    aggregating it under-counts those two classes (49,741 vs 51,103 for LTR).
  * `k` is read from `output/observed_TEs_in_windows_10kb.tsv.zst`, built by
    `05a_build_windows.sh`. Its six class counts sum to 582,540 and match the
    published Table 1 row for row, which is what licenses using the same pipeline
    for the 5 kb and 20 kb tables.

The random OR mean and SD come from the compacted permutation store
(`revision_lib.permutation_totals(by="class_name")`), the losslessly compacted form
of the same N = 500 background the published table used.

Usage
-----
    python revision_G3/10_tables.py
    python revision_G3/10_tables.py --no-verify   # skip the published-value check
"""

from __future__ import annotations

import argparse
import os
import shlex
import subprocess
import sys

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import revision_lib as rl  # noqa: E402

OUTPUT_DIR = rl.OUTPUT_DIR

# Read out of the submitted manuscript's Table 1 (word/document.xml, the single
# table in the document). Used only to assert that nothing moved.
#   class: (TE count in TSS, TE count total, OR, random OR mean, random OR SD,
#           observed/random fold change)
PUBLISHED_TABLE1 = {
    "LINE":     (169930, 1005214, 2.13, 2.43, 0.009, 0.877),
    "LTR":      (51103,   531410, 1.11, 1.67, 0.010, 0.667),
    "SINE":     (302480, 1706485, 2.25, 1.53, 0.005, 1.468),
    "DNA":      (57684,   458177, 1.51, 1.61, 0.010, 0.938),
    "SVA":      (1170,      6274, 2.40, 1.75, 0.094, 1.368),
    "Helitrons":(173,       1869, 1.07, 1.61, 0.163, 0.661),
}

# The manuscript labels the Retroposon class "SVA" and the RC class "Helitrons".
CLASS_LABEL = {
    "LINE": "LINE",
    "LTR": "LTR",
    "SINE": "SINE",
    "DNA": "DNA",
    "Retroposon": "SVA",
    "RC": "Helitrons",
}
CLASS_ORDER = ["LINE", "LTR", "SINE", "DNA", "Retroposon", "RC"]


def format_scientific(value: float) -> str:
    """G3-friendly scientific notation, fixing the malformed `9.3*10-133` style (G8).

    The submitted table writes `9.3*10-133` and `<10-200`, which are neither
    superscripted nor unambiguous. This emits `9.3 x 10^-133` and `< 10^-200`,
    which the journal can typeset, and keeps the `<10^-200` floor where the Fisher
    p underflows double precision.
    """
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return ""
    if value == 0:
        return "< 10^-200"
    if value >= 0.001:
        return f"{value:.3g}"
    exponent = int(np.floor(np.log10(value)))
    mantissa = value / 10**exponent
    return f"{mantissa:.1f} x 10^{exponent}"


# Verbatim from the frozen TEs_mapped_on_TSS_analysis.ipynb cell 55 ("THIS IS
# TABLE 1"). N_TSS is the base-pair span covered by the 10 kb TSS neighbourhoods.
N_GENOME = 3_117_275_501
N_TSS = 272_233_268


def enrichment_2x2(k: float, K: float):
    """The published contingency table: element counts against base-pair counts.

    Reproduced exactly as the frozen notebook built it, including the mixed units,
    because D1 requires the published values to stand unchanged.
    """
    return [[k, N_TSS - k], [K - k, (N_GENOME - N_TSS) - (K - k)]]


def genome_class_counts() -> dict[str, int]:
    """K per class: every annotated TE in the genome, by class."""
    out = subprocess.run(
        f"awk -F'\\t' '{{ n[$7]++ }} END {{ for (k in n) print k \"\\t\" n[k] }}' "
        f"{shlex.quote(rl.REPEATS_BED)}",
        shell=True, check=True, capture_output=True, text=True,
        env=dict(os.environ, LC_ALL="C"),
    ).stdout
    return {
        line.split("\t")[0]: int(line.split("\t")[1])
        for line in out.strip().split("\n") if line
    }


def observed_class_counts(window: str = "10kb") -> dict[str, int]:
    """k per class: TEs intersecting the TSS neighbourhoods, from 05a's output."""
    path = os.path.join(OUTPUT_DIR, f"observed_TEs_in_windows_{window}.tsv.zst")
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"{path} missing — run 05a_build_windows.sh first"
        )
    out = subprocess.run(
        f"zstd -dc {shlex.quote(path)} "
        f"| awk -F'\\t' '{{ n[$4]++ }} END {{ for (k in n) print k \"\\t\" n[k] }}'",
        shell=True, check=True, capture_output=True, text=True,
        env=dict(os.environ, LC_ALL="C"),
    ).stdout
    return {
        line.split("\t")[0]: int(line.split("\t")[1])
        for line in out.strip().split("\n") if line
    }


def class_level_enrichment(window: str = "10kb") -> pd.DataFrame:
    """Rebuild the class-level enrichment table with the published 2x2."""
    genome_counts = genome_class_counts()
    observed_counts = observed_class_counts(window)

    rows = []
    for class_name in CLASS_ORDER:
        K = genome_counts[class_name]
        k = observed_counts[class_name]
        odds_ratio, fisher_p = stats.fisher_exact(
            enrichment_2x2(k, K), alternative="two-sided"
        )
        rows.append(
            {
                "class_name": class_name,
                "class_label": CLASS_LABEL[class_name],
                "observed_tes": k,
                "total_tes": K,
                "observed_or": odds_ratio,
                "fisher_p_raw": fisher_p,
            }
        )

    table = pd.DataFrame(rows)
    # BH across the six classes — the stored family-level adjusted p is corrected
    # across 44 families, so the class-level value has to be recomputed.
    table["fisher_p_adjusted"] = multipletests(
        table["fisher_p_raw"], method="fdr_bh"
    )[1]
    return table


def random_or_distribution(table: pd.DataFrame, window: str = "10kb") -> dict:
    """The N=500 random OR per class, using the identical 2x2 with k -> k_random."""
    totals = rl.permutation_totals(window=window, by="class_name")
    wide = totals.pivot(index="seed", columns="class_name", values="n").fillna(0)
    genome_totals = table.set_index("class_name")["total_tes"].to_dict()

    random_ors = {}
    for class_name in CLASS_ORDER:
        K = genome_totals[class_name]
        k_random = wide[class_name].to_numpy(dtype=float)
        a = k_random
        b = N_TSS - k_random
        c = K - k_random
        d = (N_GENOME - N_TSS) - (K - k_random)
        with np.errstate(divide="ignore", invalid="ignore"):
            random_ors[class_name] = (a * d) / (b * c)
    return random_ors


def random_or_background(table: pd.DataFrame, window: str = "10kb") -> pd.DataFrame:
    random_ors = random_or_distribution(table, window)
    return pd.DataFrame(
        {
            "class_name": CLASS_ORDER,
            "random_or_mean": [np.nanmean(random_ors[c]) for c in CLASS_ORDER],
            "random_or_sd": [np.nanstd(random_ors[c], ddof=1) for c in CLASS_ORDER],
            "n_permutations": [int(np.isfinite(random_ors[c]).sum()) for c in CLASS_ORDER],
        }
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--no-verify", action="store_true")
    args = parser.parse_args()

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    table = class_level_enrichment()
    background = random_or_background(table)
    table = table.merge(background, on="class_name")
    table["or_observed_to_random"] = table["observed_or"] / table["random_or_mean"]

    # Empirical p, in the exact convention the published tables use. Decoded from
    # the stored family table, where `Random_Odds_Ratio_Quantile_Of_Observed` = q and
    # `p_raw_empirical` = 2 * min(q, 1-q) — e.g. Merlin q = 0.998 -> p = 0.004, and
    # CR1 q = 0.004 -> p = 0.008. When the observed OR falls outside the whole random
    # distribution q is 0 or 1 and the raw p is 0, which the published pipeline
    # clamps to 2/(N+1) = 2/501 = 0.003992, quoted throughout as 0.004. That clamp is
    # why every class in Table 1 shows the same empirical p, and it is unchanged by
    # this revision (D1).
    random_ors = random_or_distribution(table)
    quantiles, raw, clamped = [], [], []
    for _, row in table.iterrows():
        null_or = random_ors[row["class_name"]]
        null_or = null_or[np.isfinite(null_or)]
        q = float((null_or < row["observed_or"]).mean())
        p = 2 * min(q, 1 - q)
        quantiles.append(q)
        raw.append(p)
        clamped.append(max(p, 2 / (len(null_or) + 1)))
    table["random_or_quantile_of_observed"] = quantiles
    table["empirical_p_raw_unclamped"] = raw
    table["empirical_p_raw"] = clamped
    table["empirical_p_adjusted"] = multipletests(
        table["empirical_p_raw"], method="fdr_bh"
    )[1]

    # ---------------- verification against the submitted Table 1 --------------
    print("Reproduced class-level enrichment vs the submitted Table 1:")
    print(f"{'class':10s} {'TSS count':>10s} {'total':>10s} {'OR':>6s} "
          f"{'randOR':>7s} {'SD':>7s} {'fold':>6s}   match")
    mismatches = []
    for _, row in table.iterrows():
        label = row["class_label"]
        pub = PUBLISHED_TABLE1[label]
        got = (
            int(row["observed_tes"]), int(row["total_tes"]), row["observed_or"],
            row["random_or_mean"], row["random_or_sd"], row["or_observed_to_random"],
        )
        ok = (
            got[0] == pub[0]
            and got[1] == pub[1]
            and abs(got[2] - pub[2]) <= 0.005
            and abs(got[3] - pub[3]) <= 0.005
            and abs(got[4] - pub[4]) <= 0.0005
            and abs(got[5] - pub[5]) <= 0.005
        )
        if not ok:
            mismatches.append((label, pub, got))
        print(
            f"{label:10s} {got[0]:>10,} {got[1]:>10,} {got[2]:>6.2f} "
            f"{got[3]:>7.2f} {got[4]:>7.3f} {got[5]:>6.3f}   {'OK' if ok else 'MISMATCH'}"
        )

    if mismatches and not args.no_verify:
        print("\nMISMATCH against the published Table 1:")
        for label, pub, got in mismatches:
            print(f"  {label}: published {pub}")
            print(f"  {label}: computed  {tuple(round(g, 4) if isinstance(g, float) else g for g in got)}")
        print(
            "\nD1 requires these values to be unchanged. Resolve before emitting the "
            "split tables.",
            file=sys.stderr,
        )
        return 1

    # ---------------- emit the split tables ----------------------------------
    table1 = pd.DataFrame(
        {
            "Class": table["class_label"],
            "TEs in TSS windows": table["observed_tes"].map("{:,}".format),
            "TEs total": table["total_tes"].map("{:,}".format),
            "Observed OR": table["observed_or"].map("{:.2f}".format),
            "Observed/random OR": table["or_observed_to_random"].map("{:.3f}".format),
        }
    )
    table2 = pd.DataFrame(
        {
            "Class": table["class_label"],
            "Adjusted Fisher p": table["fisher_p_adjusted"].map(format_scientific),
            "Random OR (mean +- SD)": [
                f"{m:.2f} +- {s:.3f}"
                for m, s in zip(table["random_or_mean"], table["random_or_sd"])
            ],
            "Empirical p": table["empirical_p_raw"].map("{:.3f}".format),
            "Adjusted empirical p": table["empirical_p_adjusted"].map("{:.3f}".format),
        }
    )

    full = table.copy()
    full["fisher_p_raw_formatted"] = full["fisher_p_raw"].map(format_scientific)
    full["fisher_p_adjusted_formatted"] = full["fisher_p_adjusted"].map(format_scientific)

    table1.to_csv(os.path.join(OUTPUT_DIR, "Table1.csv"), index=False)
    table2.to_csv(os.path.join(OUTPUT_DIR, "Table2.csv"), index=False)
    full.to_csv(
        os.path.join(OUTPUT_DIR, "TableS_class_enrichment_full.csv"), index=False
    )

    print("\n--- Table 1. Enrichment of TE classes in gene TSS 10 kb neighbourhoods ---")
    print(table1.to_string(index=False))
    print("\n--- Table 2. Statistical support for TE class enrichment ---")
    print(table2.to_string(index=False))
    print(
        "\nWrote Table1.csv, Table2.csv and TableS_class_enrichment_full.csv "
        "(the latter keeps the raw Fisher p, per G3's statistics policy)."
    )
    print(
        "\nNote for the caption: the empirical p is at its floor of "
        f"2/(500+1) = {rl.EMPIRICAL_P_FLOOR:.6f}, quoted as 0.004, for every class — "
        "the observed OR lies outside the entire N = 500 random distribution."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
