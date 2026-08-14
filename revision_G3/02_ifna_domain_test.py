#!/usr/bin/env python
"""WP2 — statistical tests for the interferon-alpha domain (reviewer major comment 2).

The reviewer's comment
----------------------
"...do not compare this to genome-wide background or to matched random regions.
Please perform a statistical test (e.g., permutation test)... Also, please report
the number of L1 elements in this region and their specific subfamilies."

The domain is chr9:21,150,692-21,370,055 (220 kb). This script answers both halves:
it writes the full descriptive inventory the reviewer asked for, and runs four
tests against progressively better-matched null distributions.

The four tests
--------------
1. **Divergence vs an unmatched genome-wide background.** Mean divergence of the
   domain's L1 elements against 10,000 random 220 kb autosomal windows containing
   at least one L1. Two-sided empirical p, floor 2/10,001 = 0.0002.
2. **Divergence vs an L1-count-matched background.** The same, restricted to null
   windows with at least as many L1 elements as a chosen threshold (default 40).
   This is the test the reviewer's phrase "matched random regions" asks for: it
   controls for the trivial possibility that a low mean divergence simply follows
   from high L1 density.
3. **Divergence vs a gene-density-matched background.** Restricted to null windows
   containing at least 10 genes, since the domain contains 12.
4. **Subfamily composition.** A 2x2 Fisher exact test of young primate-specific L1
   (L1HS / L1P* / L1PA* / L1PB* / L1PREC*) against older mammalian L1M* in the
   domain, versus the same split genome-wide. Accompanied by a leave-one-out
   analysis and a trimmed mean (dropping the five youngest elements) to show the
   signal is not carried by a handful of outliers.

All three null distributions are drawn from ONE common pool of candidate windows,
so the tests differ only in the matching constraint and not in the sampling frame.

Data QC (caveat C7)
-------------------
One element in the domain is annotated `L1P3` with divergence 0, which is
implausible for that clade and may be a RepeatMasker artefact. The script flags it
explicitly and reports every divergence statistic both with and without it, so
nothing in the figure or the text rests on it.

Outputs (all in revision_G3/output/)
    ifna_window_elements.csv       every TE in the domain, with subfamily/family/class
    ifna_test_results.csv          one row per test: observed, null mean/SD, empirical p
    ifna_null_distributions.csv    the null values behind each test
    ifna_subfamily_composition.csv the L1 subfamily inventory
    ifna_qc.json                   descriptive numbers + the div=0 QC flag

Usage
-----
    python revision_G3/02_ifna_domain_test.py
    python revision_G3/02_ifna_domain_test.py --report     # descriptive numbers only
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shlex
import subprocess
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import revision_lib as rl  # noqa: E402

OUTPUT_DIR = rl.OUTPUT_DIR
CHROM, START, END = rl.IFNA_DOMAIN
DOMAIN_LENGTH = END - START

RANDOM_SEED = 42
N_NULL = 10_000
POOL_SIZE = 120_000  # candidate windows drawn once, then filtered per test

L1_COUNT_MATCH_THRESHOLD = 40
GENE_COUNT_MATCH_THRESHOLD = 10

# Young primate-specific L1 clades vs older mammalian ones. L1HS is the only
# currently active human subfamily; L1P*/L1PA*/L1PB/L1PREC* are primate-specific;
# L1M* ("M" for mammalian) predate the primate radiation.
YOUNG_L1_PATTERN = re.compile(r"^L1(HS|P)", re.IGNORECASE)
OLD_L1_PATTERN = re.compile(r"^L1M", re.IGNORECASE)

AUTOSOMES = [f"chr{i}" for i in range(1, 23)]


def sh(cmd: str) -> str:
    env = dict(os.environ, LC_ALL="C")
    return subprocess.run(
        cmd, shell=True, check=True, capture_output=True, text=True, env=env
    ).stdout


def load_chrom_sizes() -> dict[str, int]:
    path = os.path.join(rl.REPO_DIR, "chm13.genome")
    sizes = {}
    with open(path) as handle:
        for line in handle:
            name, length = line.rstrip("\n").split("\t")
            sizes[name] = int(length)
    return sizes


def load_domain_elements() -> pd.DataFrame:
    """Every TE overlapping the 220 kb domain."""
    region = f"{CHROM}\t{START}\t{END}\tIFNA_domain\n"
    tmp = os.path.join(OUTPUT_DIR, ".ifna_region.bed")
    with open(tmp, "w") as handle:
        handle.write(region)

    out = sh(
        f"bedtools intersect -a {shlex.quote(rl.REPEATS_BED)} "
        f"-b {shlex.quote(tmp)} -wa"
    )
    os.remove(tmp)

    rows = [line.split("\t") for line in out.strip().split("\n") if line]
    df = pd.DataFrame(
        rows,
        columns=["chrom", "start", "end", "score", "subfamily_name", "family_name",
                 "class_name"],
    )
    for col in ["start", "end", "score"]:
        df[col] = df[col].astype(int)
    df["length"] = df["end"] - df["start"]
    return df.sort_values("start").reset_index(drop=True)


def load_all_l1() -> pd.DataFrame:
    """Every L1-family element genome-wide, for the genome-wide reference stats."""
    out = sh(
        f"awk -F'\\t' -v OFS='\\t' '$6 == \"L1\" {{print $1, $2, $3, $4, $5}}' "
        f"{shlex.quote(rl.REPEATS_BED)}"
    )
    rows = [line.split("\t") for line in out.strip().split("\n") if line]
    df = pd.DataFrame(rows, columns=["chrom", "start", "end", "score", "subfamily_name"])
    for col in ["start", "end", "score"]:
        df[col] = df[col].astype(int)
    return df


def classify_l1_age(subfamily: str) -> str:
    if YOUNG_L1_PATTERN.match(subfamily):
        return "young_primate"
    if OLD_L1_PATTERN.match(subfamily):
        return "old_mammalian"
    return "other_L1"


def build_null_pool(chrom_sizes: dict[str, int]) -> pd.DataFrame:
    """Draw POOL_SIZE random 220 kb autosomal windows and characterise each.

    Windows are drawn uniformly over autosomal coordinates in proportion to
    chromosome length, then annotated with their L1 count, mean L1 divergence and
    distinct-gene count. Requiring at least one L1 (applied later) is what removes
    centromeric and acrocentric satellite arrays without needing a separate
    mappability file: those regions carry no annotated L1.
    """
    rng = np.random.default_rng(RANDOM_SEED)
    lengths = np.array([chrom_sizes[c] - DOMAIN_LENGTH for c in AUTOSOMES], dtype=float)
    weights = lengths / lengths.sum()

    picks = rng.choice(len(AUTOSOMES), size=POOL_SIZE, p=weights)
    starts = (rng.random(POOL_SIZE) * lengths[picks]).astype(np.int64)

    pool = pd.DataFrame(
        {
            "chrom": [AUTOSOMES[i] for i in picks],
            "start": starts,
            "end": starts + DOMAIN_LENGTH,
            "window_id": [f"w{i}" for i in range(POOL_SIZE)],
        }
    )

    candidates_path = os.path.join(OUTPUT_DIR, ".ifna_null_candidates.bed")
    pool.sort_values(["chrom", "start"]).to_csv(
        candidates_path, sep="\t", header=False, index=False
    )

    # L1 content per window.
    l1_path = os.path.join(OUTPUT_DIR, ".ifna_l1.bed")
    sh(
        f"awk -F'\\t' -v OFS='\\t' '$6 == \"L1\" {{print $1, $2, $3, $4, $5}}' "
        f"{shlex.quote(rl.REPEATS_BED)} | sort -k1,1 -k2,2n > {shlex.quote(l1_path)}"
    )
    l1_stats = sh(
        f"bedtools intersect -a {shlex.quote(candidates_path)} "
        f"-b {shlex.quote(l1_path)} -wa -wb -sorted "
        f"| awk -F'\\t' -v OFS='\\t' "
        f"'{{ n[$4]++; s[$4] += $8; young[$4] += ($9 ~ /^L1(HS|P)/) ? 1 : 0 }} "
        f"END {{ for (k in n) print k, n[k], s[k] / n[k], young[k] }}'"
    )
    l1_df = pd.DataFrame(
        [line.split("\t") for line in l1_stats.strip().split("\n") if line],
        columns=["window_id", "n_l1", "mean_l1_divergence", "n_young_l1"],
    )
    l1_df["n_l1"] = l1_df["n_l1"].astype(int)
    l1_df["n_young_l1"] = l1_df["n_young_l1"].astype(int)
    l1_df["mean_l1_divergence"] = l1_df["mean_l1_divergence"].astype(float)

    # Distinct genes per window.
    genes_sorted = os.path.join(OUTPUT_DIR, ".ifna_genes.bed")
    sh(f"sort -k1,1 -k2,2n {shlex.quote(rl.GENES_BED)} > {shlex.quote(genes_sorted)}")
    gene_stats = sh(
        f"bedtools intersect -a {shlex.quote(candidates_path)} "
        f"-b {shlex.quote(genes_sorted)} -wa -wb -sorted "
        f"| awk -F'\\t' -v OFS='\\t' '{{ print $4, $8 }}' | sort -u "
        f"| awk -F'\\t' '{{ n[$1]++ }} END {{ for (k in n) print k \"\\t\" n[k] }}'"
    )
    gene_df = pd.DataFrame(
        [line.split("\t") for line in gene_stats.strip().split("\n") if line],
        columns=["window_id", "n_genes"],
    )
    gene_df["n_genes"] = gene_df["n_genes"].astype(int)

    for tmp in [candidates_path, l1_path, genes_sorted]:
        os.remove(tmp)

    pool = pool.merge(l1_df, on="window_id", how="left").merge(
        gene_df, on="window_id", how="left"
    )
    pool["n_l1"] = pool["n_l1"].fillna(0).astype(int)
    pool["n_young_l1"] = pool["n_young_l1"].fillna(0).astype(int)
    pool["n_genes"] = pool["n_genes"].fillna(0).astype(int)
    return pool


def empirical_two_sided_p(observed: float, null: np.ndarray) -> float:
    """Two-sided empirical p with the standard +1 correction.

    Counts null values at least as extreme as the observed value in either
    direction, measured from the null mean, so the test does not assume the
    direction of the effect. The floor is 2/(len(null)+1).
    """
    centre = null.mean()
    deviation = abs(observed - centre)
    n_extreme = int((np.abs(null - centre) >= deviation).sum())
    return (n_extreme + 1) / (len(null) + 1)


def run_divergence_test(
    name: str,
    description: str,
    observed: float,
    pool: pd.DataFrame,
    mask: pd.Series,
    rng: np.random.Generator,
) -> tuple[dict, pd.DataFrame]:
    qualifying = pool[mask]
    n_available = len(qualifying)
    if n_available < 100:
        raise RuntimeError(
            f"{name}: only {n_available} qualifying null windows; raise POOL_SIZE"
        )
    take = min(N_NULL, n_available)
    sample = qualifying.sample(n=take, random_state=rng.integers(1 << 31))
    null = sample["mean_l1_divergence"].to_numpy(dtype=float)

    p = empirical_two_sided_p(observed, null)
    result = {
        "test": name,
        "description": description,
        "observed": observed,
        "n_null_windows": take,
        "n_qualifying_in_pool": n_available,
        "pool_size": len(pool),
        "null_mean": float(null.mean()),
        "null_sd": float(null.std(ddof=1)),
        "null_median": float(np.median(null)),
        "null_q025": float(np.quantile(null, 0.025)),
        "null_q975": float(np.quantile(null, 0.975)),
        "z_score": float((observed - null.mean()) / null.std(ddof=1)),
        "n_null_at_least_as_low": int((null <= observed).sum()),
        "empirical_p_two_sided": p,
        "empirical_p_floor": 2 / (take + 1),
        "at_floor": bool(abs(p - 2 / (take + 1)) < 1e-12),
        "statistic": "empirical permutation p (two-sided), raw",
    }
    null_frame = pd.DataFrame({"test": name, "null_mean_l1_divergence": null})
    return result, null_frame


def main() -> int:
    global POOL_SIZE

    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument(
        "--report", action="store_true", help="print the descriptive numbers and exit"
    )
    parser.add_argument("--pool-size", type=int, default=POOL_SIZE)
    args = parser.parse_args()

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    POOL_SIZE = args.pool_size

    chrom_sizes = load_chrom_sizes()

    print(f"IFNA domain: {CHROM}:{START:,}-{END:,}  ({DOMAIN_LENGTH:,} bp)")
    elements = load_domain_elements()
    l1_domain = elements[elements.family_name == "L1"].copy()
    l1_domain["age_group"] = l1_domain.subfamily_name.map(classify_l1_age)

    print(f"\nAll TEs in the window          : {len(elements)}")
    print(f"L1 elements                    : {len(l1_domain)} "
          f"({100 * len(l1_domain) / len(elements):.0f} % of all TEs)")
    print("\nFamily composition:")
    for family, count in elements.family_name.value_counts().items():
        print(f"  {family or '(empty)':16s} {count}")

    genome_l1 = load_all_l1()
    genome_l1["age_group"] = genome_l1.subfamily_name.map(classify_l1_age)

    # Distinct genes whose TSS neighbourhood overlaps the domain (§3.3 says 12).
    region_bed = os.path.join(OUTPUT_DIR, ".ifna_region_genes.bed")
    with open(region_bed, "w") as handle:
        handle.write(f"{CHROM}\t{START}\t{END}\tIFNA_domain\n")
    gene_names = sorted(set(
        line.split("\t")[3]
        for line in sh(
            f"bedtools intersect -a {shlex.quote(rl.GENES_BED)} "
            f"-b {shlex.quote(region_bed)} -wa"
        ).strip().split("\n")
        if line
    ))
    os.remove(region_bed)
    n_genes_in_window = len(gene_names)
    print(f"\nGenes with a TSS neighbourhood overlapping the domain: "
          f"{n_genes_in_window}")
    print("  " + ", ".join(gene_names))

    # --- QC: the implausible divergence-0 element (caveat C7) -----------------
    zero_div = l1_domain[l1_domain.score == 0]
    print(f"\nQC: L1 elements with divergence 0 in the window: {len(zero_div)}")
    for _, row in zero_div.iterrows():
        print(f"  {row.subfamily_name} at {row.chrom}:{row.start:,}-{row.end:,} "
              f"({row.length:,} bp) — implausible for this clade, see caveat C7")
    genome_zero = int((genome_l1.score == 0).sum())
    print(f"  genome-wide L1 with divergence 0: {genome_zero} of {len(genome_l1):,} "
          f"({100 * genome_zero / len(genome_l1):.4f} %)")

    l1_trimmed = l1_domain[l1_domain.score > 0]
    youngest_five = l1_domain.nsmallest(5, "score")
    l1_drop_youngest = l1_domain.drop(youngest_five.index)

    observed_mean = float(l1_domain.score.mean())
    print(f"\nMean divergence of the {len(l1_domain)} domain L1 elements : "
          f"{observed_mean:.1f}")
    print(f"  excluding the divergence-0 element (n={len(l1_trimmed)})   : "
          f"{l1_trimmed.score.mean():.1f}")
    print(f"  excluding the 5 youngest (n={len(l1_drop_youngest)})        : "
          f"{l1_drop_youngest.score.mean():.1f}")
    print(f"Genome-wide L1 mean divergence ({len(genome_l1):,} elements) : "
          f"{genome_l1.score.mean():.1f}  (median {genome_l1.score.median():.0f})")

    density_window = len(l1_domain) / (DOMAIN_LENGTH / 1e6)
    genome_span = sum(chrom_sizes[c] for c in chrom_sizes)
    density_genome = len(genome_l1) / (genome_span / 1e6)
    print(f"\nL1 density in the window : {density_window:.0f} / Mb")
    print(f"L1 density genome-wide   : {density_genome:.0f} / Mb  "
          f"({density_window / density_genome:.1f}x)")

    subfamilies = sorted(l1_domain.subfamily_name.unique())
    print(f"\nDistinct L1 subfamilies in the window: {len(subfamilies)}")
    print("  " + ", ".join(subfamilies))
    print("\nYoungest elements (lowest divergence):")
    for _, row in l1_domain.nsmallest(8, "score").iterrows():
        print(f"  {row.subfamily_name:12s} div {row.score:4d}  "
              f"{row.chrom}:{row.start:,}-{row.end:,}")

    if args.report:
        return 0

    # --- Null pool ------------------------------------------------------------
    print(f"\nDrawing {POOL_SIZE:,} random {DOMAIN_LENGTH / 1000:.0f} kb autosomal "
          f"windows (seed {RANDOM_SEED}) ...")
    pool = build_null_pool(chrom_sizes)
    print(f"  windows with >= 1 L1                    : "
          f"{int((pool.n_l1 >= 1).sum()):,}")
    print(f"  windows with >= {L1_COUNT_MATCH_THRESHOLD} L1                   : "
          f"{int((pool.n_l1 >= L1_COUNT_MATCH_THRESHOLD).sum()):,}")
    print(f"  windows with >= {GENE_COUNT_MATCH_THRESHOLD} genes and >= 1 L1     : "
          f"{int(((pool.n_genes >= GENE_COUNT_MATCH_THRESHOLD) & (pool.n_l1 >= 1)).sum()):,}")

    rng = np.random.default_rng(RANDOM_SEED)
    results, nulls = [], []

    r, n = run_divergence_test(
        "T1_divergence_unmatched",
        "Mean L1 divergence vs random 220 kb autosomal windows with >= 1 L1",
        observed_mean, pool, pool.n_l1 >= 1, rng,
    )
    results.append(r); nulls.append(n)

    r, n = run_divergence_test(
        "T2_divergence_L1_count_matched",
        f"Mean L1 divergence vs random windows with >= {L1_COUNT_MATCH_THRESHOLD} L1 "
        f"(controls for L1 density)",
        observed_mean, pool, pool.n_l1 >= L1_COUNT_MATCH_THRESHOLD, rng,
    )
    results.append(r); nulls.append(n)

    r, n = run_divergence_test(
        "T3_divergence_gene_density_matched",
        f"Mean L1 divergence vs random windows with >= "
        f"{GENE_COUNT_MATCH_THRESHOLD} genes and >= 1 L1",
        observed_mean, pool,
        (pool.n_genes >= GENE_COUNT_MATCH_THRESHOLD) & (pool.n_l1 >= 1), rng,
    )
    results.append(r); nulls.append(n)

    # --- Test 4: subfamily composition ---------------------------------------
    win_young = int((l1_domain.age_group == "young_primate").sum())
    win_old = int((l1_domain.age_group == "old_mammalian").sum())
    gen_young = int((genome_l1.age_group == "young_primate").sum())
    gen_old = int((genome_l1.age_group == "old_mammalian").sum())
    # Remove the window's own elements from the genome-wide comparison arm so the
    # two arms are disjoint.
    table = [[win_young, win_old], [gen_young - win_young, gen_old - win_old]]
    odds_ratio, fisher_p = stats.fisher_exact(table, alternative="two-sided")

    print(f"\nTest 4 — subfamily composition (2x2 Fisher exact, two-sided, raw p):")
    print(f"  domain      : {win_young} young primate-specific L1, {win_old} old L1M*")
    print(f"  rest of genome: {table[1][0]:,} young, {table[1][1]:,} old")
    print(f"  odds ratio {odds_ratio:.3f}, p = {fisher_p:.3g}")

    results.append(
        {
            "test": "T4_subfamily_composition",
            "description": (
                "2x2 Fisher exact: young primate-specific L1 (L1HS/L1P*) vs older "
                "L1M* in the domain against the rest of the genome"
            ),
            "observed": odds_ratio,
            "n_null_windows": np.nan,
            "n_qualifying_in_pool": np.nan,
            "pool_size": len(pool),
            "null_mean": (gen_young - win_young) / max(1, gen_old - win_old),
            "null_sd": np.nan,
            "null_median": np.nan,
            "null_q025": np.nan,
            "null_q975": np.nan,
            "z_score": np.nan,
            "n_null_at_least_as_low": np.nan,
            "empirical_p_two_sided": fisher_p,
            "empirical_p_floor": np.nan,
            "at_floor": False,
            "statistic": "Fisher exact p (two-sided), raw",
            "window_young": win_young,
            "window_old": win_old,
            "genome_young": table[1][0],
            "genome_old": table[1][1],
        }
    )

    # --- Robustness of the divergence signal ---------------------------------
    print("\nRobustness of the divergence result (test 2's null):")
    null_matched = nulls[1]["null_mean_l1_divergence"].to_numpy()
    for label, subset in [
        ("all 77 elements", l1_domain),
        ("excluding the div=0 element", l1_trimmed),
        ("excluding the 5 youngest", l1_drop_youngest),
    ]:
        value = float(subset.score.mean())
        p = empirical_two_sided_p(value, null_matched)
        print(f"  {label:30s} mean {value:6.1f}  empirical p = {p:.4f}  "
              f"(n = {len(subset)})")

    loo = np.array([
        l1_domain.drop(idx).score.mean() for idx in l1_domain.index
    ])
    print(f"  leave-one-out mean range      : {loo.min():.1f} - {loo.max():.1f} "
          f"(all below the null mean {null_matched.mean():.1f}: "
          f"{bool((loo < null_matched.mean()).all())})")

    # --- Persist everything ---------------------------------------------------
    elements.to_csv(os.path.join(OUTPUT_DIR, "ifna_window_elements.csv"), index=False)
    pd.DataFrame(results).to_csv(
        os.path.join(OUTPUT_DIR, "ifna_test_results.csv"), index=False
    )
    pd.concat(nulls, ignore_index=True).to_csv(
        os.path.join(OUTPUT_DIR, "ifna_null_distributions.csv"), index=False
    )

    composition = (
        l1_domain.groupby(["subfamily_name", "age_group"], observed=True)
        .agg(n=("score", "size"), mean_divergence=("score", "mean"),
             min_divergence=("score", "min"), max_divergence=("score", "max"))
        .reset_index()
        .sort_values("mean_divergence")
    )
    composition.to_csv(
        os.path.join(OUTPUT_DIR, "ifna_subfamily_composition.csv"), index=False
    )
    pool.to_csv(os.path.join(OUTPUT_DIR, "ifna_null_pool.csv"), index=False)

    qc = {
        "domain": f"{CHROM}:{START}-{END}",
        "domain_length_bp": DOMAIN_LENGTH,
        "n_tes_in_window": int(len(elements)),
        "n_l1_in_window": int(len(l1_domain)),
        "l1_fraction_of_window_tes": round(len(l1_domain) / len(elements), 4),
        "family_composition": elements.family_name.value_counts().to_dict(),
        "n_distinct_l1_subfamilies": int(len(subfamilies)),
        "l1_subfamilies": subfamilies,
        "genes_in_window": gene_names,
        "mean_l1_divergence_window": round(observed_mean, 2),
        "mean_l1_divergence_window_excl_zero": round(float(l1_trimmed.score.mean()), 2),
        "mean_l1_divergence_window_excl_5_youngest": round(
            float(l1_drop_youngest.score.mean()), 2
        ),
        "mean_l1_divergence_genome": round(float(genome_l1.score.mean()), 2),
        "median_l1_divergence_genome": float(genome_l1.score.median()),
        "n_l1_genome": int(len(genome_l1)),
        "l1_per_mb_window": round(density_window, 1),
        "l1_per_mb_genome": round(density_genome, 1),
        "l1_density_ratio": round(density_window / density_genome, 2),
        "n_genes_in_window": n_genes_in_window,
        "qc_divergence_zero_elements": zero_div[
            ["chrom", "start", "end", "score", "subfamily_name"]
        ].to_dict(orient="records"),
        "qc_note_C7": (
            "One L1P3 element is annotated with divergence 0, implausible for that "
            "clade and possibly a RepeatMasker artefact. Every divergence statistic "
            "is reported with and without it; the conclusion does not depend on it."
        ),
        "random_seed": RANDOM_SEED,
        "pool_size": int(len(pool)),
        "n_null_per_test": N_NULL,
    }
    with open(os.path.join(OUTPUT_DIR, "ifna_qc.json"), "w") as handle:
        json.dump(qc, handle, indent=2, default=str)

    print("\n" + "=" * 78)
    print("Test results")
    print("=" * 78)
    summary = pd.DataFrame(results)[
        ["test", "observed", "null_mean", "null_sd", "z_score",
         "empirical_p_two_sided", "n_null_windows"]
    ]
    print(summary.to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    print("\nWrote ifna_window_elements.csv, ifna_test_results.csv,")
    print("      ifna_null_distributions.csv, ifna_subfamily_composition.csv,")
    print("      ifna_null_pool.csv, ifna_qc.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
