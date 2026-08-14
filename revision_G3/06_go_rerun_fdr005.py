#!/usr/bin/env python
"""WP6 — re-run every GO enrichment level at FDR 0.05 (reviewer minor comment 2, D2).

Decision D2 is 0.05 everywhere, with NO "suggestive 0.05-0.1" band in the main
text, the supplementary figures or the supplementary tables. This script is the
single place that threshold is applied for the G3 revision.

What is reproduced exactly, and what changes
--------------------------------------------
The foreground gene sets are built by the SAME construction as the frozen
`Gene_ontology_analysis.ipynb`, quirks included, so that the only difference
between the published tables and these is the FDR cut:

  * classes by count (frozen cell 20): LINE / LTR / SINE / DNA take the top
    **1,436** genes by `{class}_number` — a hard-coded number in the original,
    not a recomputed 5 %, and it is kept hard-coded here for that reason.
    SVA and Helitron are NOT top-5 % sets: they are ALL genes with any
    Retroposon / RC element, which is why their gene counts differ from the rest.
    TE_top / TE_bottom are the top and bottom 1,436 genes by `TE_number`.
  * classes by divergence (frozen cell 90): sort by `Divergence_Avg_{class}`
    (or `Average_Divergence_Score` for TE_all), drop NaN, then take the top and
    bottom `len // 20`.
  * families by count (frozen cell 149): for each of the 44 curated families,
    genes with a non-zero `{family}_number`, sorted descending, truncated to
    `n_unique_genes // 20`.
  * subfamilies (frozen `GO_subfamilies.py`): the same, per subfamily.

Nothing is written outside `revision_G3/output/`. In particular the shared
`../GO_tables/` directory is left at FDR 0.1, because `GO_subfamilies.py` belongs
to the companion subfamilies manuscript whose refresh is deferred (§7.4, C3).

Retrieval at 0.1, publication at 0.05
-------------------------------------
Each level is retrieved at `fdr_threshold=0.1` — which reproduces the published
term set exactly — and then written out twice:

    GO_<level>_fdr01_reference.csv   every term with FDR < 0.1  (NOT for publication)
    GO_<level>_fdr005.csv            every term with FDR < 0.05 (this is what ships)

The 0.1 reference file exists only so the effect of tightening the threshold is
measured rather than asserted: `--report` prints the term counts at both cuts and
names every term that leaves the paper, which is what plan §3.2 tabulates and
what the response letter quotes. No figure or supplementary table is built from
the 0.1 file.

Usage
-----
    python revision_G3/06_go_rerun_fdr005.py --level classes_count
    python revision_G3/06_go_rerun_fdr005.py --level all          # excludes subfamilies
    python revision_G3/06_go_rerun_fdr005.py --level subfamilies  # ~1,143 studies, hours
    python revision_G3/06_go_rerun_fdr005.py --report
"""

from __future__ import annotations

import argparse
import os
import sys
import time
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import revision_lib as rl  # noqa: E402

OUTPUT_DIR = rl.OUTPUT_DIR
REFERENCE_FDR = 0.10  # retrieval cut, reproduces the published set
PUBLISHED_FDR = rl.FDR_THRESHOLD  # 0.05, what ships

# Hard-coded in the frozen notebook (cell 20) as the top-5 % size. Kept verbatim
# so the class-level foreground sets are identical to the published ones.
CLASS_TOP_N = 1436

MAJOR_CLASSES = ["LINE", "LTR", "SINE", "DNA"]

# Copied from the frozen Gene_ontology_analysis.ipynb cell 134/136.
FAMILY_NAMES_CURATED = [
    "I-Jockey", "TcMar-Tigger", "Kolobok", "hAT", "hAT-Tag1", "MULE-MuDR",
    "ERVL", "tRNA-RTE", "TcMar-Pogo", "hAT-Tip100", "ERVK", "Penelope",
    "L1-Tx1", "ERVL-MaLR", "Gypsy", "PiggyBac", "TcMar-Mariner", "RTE-BovB",
    "hAT-hAT19", "L2", "hAT?", "Crypton-A", "tRNA-Deu", "TcMar-Tc1", "CR1",
    "hAT-Tip100?", "hAT-Ac", "Alu", "hAT-Charlie", "5S-Deu-L2", "Crypton",
    "PIF-Harbinger", "MIR", "RTE-X", "L1", "ERV1", "Merlin", "Helitron",
    "SVA", "tRNA", "Dong-R4", "TcMar-Tc2", "TcMar", "hAT-Blackjack",
]

# Terms plan §3.2 predicts will leave the paper at 0.05. `--report` checks each.
TERMS_EXPECTED_TO_DROP = [
    "flavone metabolic process",
    "detoxification of copper ion",
    "response to cadmium ion",
    "beta-2-microglobulin binding",
    "antigen processing and presentation of endogenous peptide antigen via MHC class I",
]
TERMS_EXPECTED_TO_SURVIVE = [
    "type I interferon receptor binding",
    "termination of RNA polymerase II transcription",
    "olfactory receptor activity",
    "cellular response to zinc ion",
    "MHC class I protein complex",
]


def load_gene_table() -> pd.DataFrame:
    """The per-TSS-window table with class counts and per-class divergence."""
    path = os.path.join(rl.REPO_DIR, "TEs_on_genes.csv")
    df = pd.read_csv(path, low_memory=False)
    print(f"Loaded {os.path.basename(path)}: {len(df):,} TSS windows, "
          f"{df['Gene_name'].nunique():,} distinct genes")
    return df


def add_family_counts(df: pd.DataFrame) -> pd.DataFrame:
    """Per-window counts of each curated family, from the comma-joined column.

    Same construction as frozen cell 134: count exact matches in the
    comma-separated `TE_families` string, so `hAT` does not absorb `hAT-Ac`.
    """
    families = df["TE_families"].fillna("").astype(str)
    split_lists = families.str.split(",")
    for family_name in FAMILY_NAMES_CURATED:
        df[f"{family_name}_number"] = split_lists.apply(
            lambda names, target=family_name: sum(1 for n in names if n == target)
        )
    return df


def run_level_classes_count(df: pd.DataFrame) -> pd.DataFrame:
    """GO for the six classes plus the TE_top / TE_bottom groups (frozen cell 20)."""
    background = list(df["Gene_name"].unique())
    frames = []

    for class_name in MAJOR_CLASSES:
        genes = list(
            df.sort_values(f"{class_name}_number", ascending=False)[["Gene_name"]]
            .drop_duplicates()["Gene_name"]
        )[:CLASS_TOP_N]
        frames.append(_study(genes, background, class_name=class_name))

    # SVA and Helitron are every gene with any such element, not a top-5 % set.
    for class_name, column in [("SVA", "Retroposon_number"), ("Helitron", "RC_number")]:
        genes = list(df[df[column] != 0]["Gene_name"].unique())
        frames.append(_study(genes, background, class_name=class_name))

    for class_name, ascending in [("TE_top", False), ("TE_bottom", True)]:
        genes = list(
            df.sort_values("TE_number", ascending=ascending)[["Gene_name"]]
            .drop_duplicates()["Gene_name"]
        )[:CLASS_TOP_N]
        frames.append(_study(genes, background, class_name=class_name))

    return pd.concat(frames, axis=0, ignore_index=True)


def run_level_classes_divergence(df: pd.DataFrame) -> pd.DataFrame:
    """GO for the highest- and lowest-divergence 5 % per class (frozen cell 90)."""
    background = list(df["Gene_name"].unique())
    frames = []

    targets = [(c, f"Divergence_Avg_{c}") for c in MAJOR_CLASSES]
    targets.append(("TE_all", "Average_Divergence_Score"))

    for class_name, column in targets:
        ranked = list(
            df.sort_values(column, ascending=False)
            .dropna(subset=[column])[["Gene_name"]]
            .drop_duplicates()["Gene_name"]
        )
        take = len(ranked) // 20
        for group, genes in [("highest", ranked[:take]), ("lowest", ranked[-take:])]:
            frames.append(
                _study(genes, background, class_name=class_name, divergence_group=group)
            )

    return pd.concat(frames, axis=0, ignore_index=True)


def run_level_families(df: pd.DataFrame) -> pd.DataFrame:
    """GO for each of the 44 curated families (frozen cell 149)."""
    families_by_classes = (
        pd.read_csv(os.path.join(rl.REPO_DIR, "families_by_classes_TE.csv"))
        .dropna()
        .set_index("family_name")
    )
    family_to_class = families_by_classes["class_name"].to_dict()

    background = list(df["Gene_name"].unique())
    take = len(background) // 20
    frames = []

    for family_name in FAMILY_NAMES_CURATED:
        column = f"{family_name}_number"
        genes = list(
            df[df[column] != 0]
            .sort_values(column, ascending=False)[["Gene_name"]]
            .drop_duplicates()["Gene_name"]
        )
        genes = genes[:take] if len(genes) >= take else genes
        if not genes:
            print(f"  {family_name}: no genes, skipped")
            continue
        frames.append(
            _study(
                genes,
                background,
                family_name=family_name,
                class_name=family_to_class.get(family_name, "Other"),
            )
        )

    return pd.concat(frames, axis=0, ignore_index=True)


def run_level_subfamilies() -> pd.DataFrame:
    """GO for every annotated subfamily (frozen GO_subfamilies.py).

    This is the long one: ~1,143 independent studies. Writes one CSV per
    subfamily into output/GO_tables_fdr005/ as it goes, so an interrupted run
    resumes rather than restarting, and the shared ../GO_tables/ is untouched.
    """
    path = os.path.join(rl.REPO_DIR, "TEs_on_genes_counts_subfamilies.csv")
    subfam_counts = pd.read_csv(path, low_memory=False)
    background = list(subfam_counts["Gene_name"].unique())
    take = len(background) // 20

    mapping = pd.read_csv(os.path.join(rl.REPO_DIR, "individuals_by_classes_TE.csv"))[
        ["class_name", "individual_name"]
    ]
    mapping.columns = ["class_name", "subfamily_name"]
    subfamily_to_class = dict(zip(mapping.subfamily_name, mapping.class_name))

    per_subfamily_dir = os.path.join(OUTPUT_DIR, "GO_tables_fdr005")
    os.makedirs(per_subfamily_dir, exist_ok=True)

    frames = []
    names = list(mapping["subfamily_name"])
    for i, subfamily_name in enumerate(names, 1):
        out_path = os.path.join(
            per_subfamily_dir, f"GO_terms_by_subfamilies_{subfamily_name}.csv"
        )
        if os.path.exists(out_path):
            existing = pd.read_csv(out_path)
            if not existing.empty:
                frames.append(existing)
            continue

        column = f"{subfamily_name}_number"
        if column not in subfam_counts.columns:
            continue
        genes = list(
            subfam_counts[subfam_counts[column] != 0]
            .sort_values(column, ascending=False)[["Gene_name"]]
            .drop_duplicates()["Gene_name"]
        )
        if not genes:
            continue
        genes = genes[:take] if len(genes) >= take else genes

        result = _study(
            genes,
            background,
            subfamily_name=subfamily_name,
            class_name=subfamily_to_class.get(subfamily_name, "Other"),
        )
        result.to_csv(out_path, index=False)
        if not result.empty:
            frames.append(result)
        if i % 25 == 0:
            print(f"  [{i}/{len(names)}] {subfamily_name}", flush=True)

    return pd.concat(frames, axis=0, ignore_index=True) if frames else pd.DataFrame()


def _study(genes, background, **annotations) -> pd.DataFrame:
    """One goatools study at the reference cut, annotated with its group labels."""
    label = " / ".join(str(v) for v in annotations.values())
    start = time.time()
    result = rl.run_goatools_enrichment(
        genes, background, fdr_threshold=REFERENCE_FDR
    )
    if result.empty:
        result = pd.DataFrame(
            columns=[
                "Term ID", "Term Name", "Term Database", "P-value", "FDR",
                "Fold Enrichment", "Overlap Count", "Total Term Genes (Human)",
                "Overlapping Genes", "Full Term Gene List",
            ]
        )
    for key, value in annotations.items():
        result[key] = value
    result["n_foreground_genes"] = len(genes)

    n_005 = int((result["FDR"] < PUBLISHED_FDR).sum()) if not result.empty else 0
    print(
        f"  {label:28s} {len(genes):>6,} genes -> "
        f"{len(result):>4} terms at FDR<{REFERENCE_FDR}, "
        f"{n_005:>4} at FDR<{PUBLISHED_FDR}   ({time.time() - start:.1f} s)",
        flush=True,
    )
    return result


def write_level(name: str, table: pd.DataFrame) -> None:
    reference_path = os.path.join(OUTPUT_DIR, f"GO_{name}_fdr01_reference.csv")
    published_path = os.path.join(OUTPUT_DIR, f"GO_{name}_fdr005.csv")

    table.to_csv(reference_path, index=False)
    published = table[table["FDR"] < PUBLISHED_FDR].copy()
    published.to_csv(published_path, index=False)

    print(
        f"\n{name}: {len(table):,} terms at FDR<{REFERENCE_FDR} -> "
        f"{len(published):,} at FDR<{PUBLISHED_FDR} "
        f"({100 * (len(published) / len(table) - 1):+.1f} %)"
    )
    print(f"  wrote {os.path.basename(reference_path)} (reference, NOT for publication)")
    print(f"  wrote {os.path.basename(published_path)} (this is what ships)")


def report() -> int:
    """Print the measured effect of 0.1 -> 0.05, i.e. plan §3.2, from real files."""
    print("=" * 78)
    print("Measured effect of tightening the GO FDR threshold from 0.1 to 0.05")
    print("=" * 78)

    levels = [
        ("classes_count", "class_name"),
        ("classes_divergence", "class_name"),
        ("families", "family_name"),
        ("subfamilies", "subfamily_name"),
    ]
    any_found = False
    for name, group_col in levels:
        path = os.path.join(OUTPUT_DIR, f"GO_{name}_fdr01_reference.csv")
        if not os.path.exists(path):
            print(f"\n{name}: not run yet")
            continue
        any_found = True
        table = pd.read_csv(path, low_memory=False)
        n01 = len(table)
        n005 = int((table.FDR < PUBLISHED_FDR).sum())
        print(f"\n{name}: {n01} terms at FDR<0.1 -> {n005} at FDR<0.05 "
              f"({100 * (n005 / n01 - 1):+.1f} %)" if n01 else f"\n{name}: empty")

        if group_col in table.columns:
            lost_groups = []
            for group, grp in table.groupby(group_col):
                if len(grp) and not (grp.FDR < PUBLISHED_FDR).any():
                    lost_groups.append(f"{group} (had {len(grp)})")
            print(f"  groups losing ALL terms: "
                  f"{', '.join(lost_groups) if lost_groups else 'none'}")

    if not any_found:
        print("\nNothing to report — run the levels first.")
        return 1

    print("\n" + "-" * 78)
    print("Headline terms named in plan §3.2:")
    all_ref = []
    for name, _ in levels:
        path = os.path.join(OUTPUT_DIR, f"GO_{name}_fdr01_reference.csv")
        if os.path.exists(path):
            t = pd.read_csv(path, low_memory=False)
            t["level"] = name
            all_ref.append(t)
    if not all_ref:
        return 1
    ref = pd.concat(all_ref, ignore_index=True)

    for term, expectation in (
        [(t, "survive") for t in TERMS_EXPECTED_TO_SURVIVE]
        + [(t, "drop") for t in TERMS_EXPECTED_TO_DROP]
    ):
        rows = ref[ref["Term Name"] == term]
        if rows.empty:
            print(f"  {term[:60]:62s} not present in any level")
            continue
        best = rows.FDR.min()
        actual = "survive" if best < PUBLISHED_FDR else "drop"
        flag = "OK " if actual == expectation else "!! "
        print(
            f"  {flag}{term[:58]:60s} min FDR = {best:.4g}  -> {actual}"
            + ("" if actual == expectation else f"  (plan predicted {expectation})")
        )
    print("-" * 78)
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument(
        "--level",
        choices=[
            "classes_count",
            "classes_divergence",
            "families",
            "subfamilies",
            "all",
        ],
        help="'all' runs the three fast levels; subfamilies must be asked for",
    )
    parser.add_argument("--report", action="store_true", help="print the §3.2 table")
    args = parser.parse_args()

    if args.report:
        return report()
    if not args.level:
        parser.error("give --level ... or --report")

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    if args.level == "subfamilies":
        print("Level: subfamilies (~1,143 studies; resumable)")
        table = run_level_subfamilies()
        if not table.empty:
            write_level("subfamilies", table)
        return 0

    df = load_gene_table()
    wanted = (
        ["classes_count", "classes_divergence", "families"]
        if args.level == "all"
        else [args.level]
    )

    if "families" in wanted:
        print("Deriving per-family counts from TE_families ...")
        df = add_family_counts(df)

    for level in wanted:
        print(f"\n=== level: {level} ===")
        if level == "classes_count":
            write_level(level, run_level_classes_count(df))
        elif level == "classes_divergence":
            write_level(level, run_level_classes_divergence(df))
        elif level == "families":
            write_level(level, run_level_families(df))

    return 0


if __name__ == "__main__":
    sys.exit(main())
