#!/usr/bin/env python
"""WP5 — GO enrichment at the top/bottom 10 % gene sets vs. the published 5 %.

The reviewer's comment
----------------------
"...sensitivity analyses for key findings (e.g., the IFNA-L1 association, the
SVA-termination association) using alternative window sizes (5 kb and 20 kb) and
alternative percentiles (e.g., 10%)."

`05b_window_sensitivity.py` covers the window half. This script covers the
percentile half: it rebuilds every foreground gene set at 10 % instead of 5 %,
re-runs GO at FDR 0.05, and measures how much the published conclusions depend on
that cutoff.

How the 10 % sets are built
---------------------------
The constructions are the frozen notebook's, reused verbatim from
`06_go_rerun_fdr005.py` (which documents each one against its source cell) with a
single parameter changed. Reusing that module rather than re-deriving the sets is
deliberate: the 5 % arm of this comparison IS the published arm, so any drift
between the two constructions would show up as a fake percentile effect.

  * classes by count      top 2 x 1,436 = 2,872 genes instead of 1,436
  * classes by divergence top and bottom len // 10 instead of len // 20
  * families by count     truncated at len(background) // 10 instead of // 20

2,872 is used rather than 28,738 // 10 = 2,873 because the published cut is the
hard-coded 1,436 and doubling it keeps the arithmetic relationship exact; the
one-gene difference is immaterial. Either way the 5 % set is a strict prefix of
the 10 % set — the script asserts that nesting, because it is what makes the
term-set comparison interpretable as "what does adding the next 5 % of genes do".

Two groups have no percentile parameter at all: SVA and Helitron are every gene
carrying any Retroposon / RC element (frozen cell 20), not a ranked cut. They are
reported as percentile-invariant instead of being re-run, and the same is true of
any family with fewer genes than the truncation limit.

Outputs
    output/GO_<level>_p10_fdr01_reference.csv     retrieval cut, not for publication
    output/GO_<level>_p10_fdr005.csv              the 10 % arm, what a table would show
    output/percentile_sensitivity_summary.csv     per group: term counts, Jaccard, drift
    output/percentile_sensitivity_terms.csv       every term gained or lost, named
    output/percentile_sensitivity_headline.csv    each abstract-level claim at 5 % and 10 %

Usage
    python revision_G3/05c_percentile_sensitivity.py                 # all three levels
    python revision_G3/05c_percentile_sensitivity.py --level families
    python revision_G3/05c_percentile_sensitivity.py --summary       # report only
"""

from __future__ import annotations

import argparse
import importlib
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import revision_lib as rl  # noqa: E402

# The module name starts with a digit, so it cannot be imported with `import`.
# Reusing it is what guarantees the 5 % and 10 % arms share one construction.
go_rerun = importlib.import_module("06_go_rerun_fdr005")

OUTPUT_DIR = rl.OUTPUT_DIR
REFERENCE_FDR = go_rerun.REFERENCE_FDR  # 0.10 retrieval cut
PUBLISHED_FDR = rl.FDR_THRESHOLD  # 0.05, what ships

PERCENTILE = 10
CLASS_TOP_N_10 = 2 * go_rerun.CLASS_TOP_N  # 2,872

LEVELS = ["classes_count", "classes_divergence", "families"]

GROUP_COLUMN = {
    "classes_count": "class_name",
    "classes_divergence": "class_name",
    "families": "family_name",
}

# Groups whose construction has no percentile parameter (frozen cell 20).
PERCENTILE_INVARIANT_GROUPS = {"SVA", "Helitron"}

# The abstract- and title-level claims, as (level, group, term). Each is looked up
# at both percentiles so the response letter can state survival rather than imply it.
HEADLINE_CLAIMS = [
    ("classes_divergence", "LINE", "type I interferon receptor binding"),
    ("classes_divergence", "LINE", "olfactory receptor activity"),
    ("families", "SVA", "termination of RNA polymerase II transcription"),
    ("families", "L1", "flavone metabolic process"),
    ("families", "MIR", "cellular response to zinc ion"),
    ("families", "MIR", "intracellular zinc ion homeostasis"),
    ("families", "hAT-Charlie", "MHC class I protein complex"),
    ("classes_count", "LINE", "olfactory receptor activity"),
    ("classes_count", "SINE", "mRNA splicing, via spliceosome"),
]


def build_classes_count(df: pd.DataFrame) -> pd.DataFrame:
    """GO for the four major classes and TE_top / TE_bottom at the 10 % cut.

    SVA and Helitron are omitted: their gene sets are not ranked cuts, so they are
    identical at any percentile and are reported as invariant instead.
    """
    background = list(df["Gene_name"].unique())
    frames = []

    for class_name in go_rerun.MAJOR_CLASSES:
        genes = list(
            df.sort_values(f"{class_name}_number", ascending=False)[["Gene_name"]]
            .drop_duplicates()["Gene_name"]
        )[:CLASS_TOP_N_10]
        frames.append(go_rerun._study(genes, background, class_name=class_name))

    for class_name, ascending in [("TE_top", False), ("TE_bottom", True)]:
        genes = list(
            df.sort_values("TE_number", ascending=ascending)[["Gene_name"]]
            .drop_duplicates()["Gene_name"]
        )[:CLASS_TOP_N_10]
        frames.append(go_rerun._study(genes, background, class_name=class_name))

    return pd.concat(frames, axis=0, ignore_index=True)


def build_classes_divergence(df: pd.DataFrame) -> pd.DataFrame:
    """GO for the highest- and lowest-divergence 10 % per class (frozen cell 90)."""
    background = list(df["Gene_name"].unique())
    frames = []

    targets = [(c, f"Divergence_Avg_{c}") for c in go_rerun.MAJOR_CLASSES]
    targets.append(("TE_all", "Average_Divergence_Score"))

    for class_name, column in targets:
        ranked = list(
            df.sort_values(column, ascending=False)
            .dropna(subset=[column])[["Gene_name"]]
            .drop_duplicates()["Gene_name"]
        )
        take = len(ranked) // PERCENTILE
        for group, genes in [("highest", ranked[:take]), ("lowest", ranked[-take:])]:
            frames.append(
                go_rerun._study(
                    genes, background, class_name=class_name, divergence_group=group
                )
            )

    return pd.concat(frames, axis=0, ignore_index=True)


def build_families(df: pd.DataFrame) -> pd.DataFrame:
    """GO for each of the 44 curated families at the 10 % truncation (frozen cell 149)."""
    families_by_classes = (
        pd.read_csv(os.path.join(rl.REPO_DIR, "families_by_classes_TE.csv"))
        .dropna()
        .set_index("family_name")
    )
    family_to_class = families_by_classes["class_name"].to_dict()

    background = list(df["Gene_name"].unique())
    take = len(background) // PERCENTILE
    frames = []

    for family_name in go_rerun.FAMILY_NAMES_CURATED:
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
            go_rerun._study(
                genes,
                background,
                family_name=family_name,
                class_name=family_to_class.get(family_name, "Other"),
            )
        )

    return pd.concat(frames, axis=0, ignore_index=True)


BUILDERS = {
    "classes_count": build_classes_count,
    "classes_divergence": build_classes_divergence,
    "families": build_families,
}


def assert_nesting(df: pd.DataFrame) -> None:
    """The 5 % foreground must be a prefix of the 10 % one, for every class.

    Both arms sort on the same column in the same direction, so the smaller set is
    a prefix of the larger by construction. Checking it here turns that argument
    into a test — if it ever fails, the two arms are not measuring the same thing
    and the whole comparison is void.
    """
    for class_name in go_rerun.MAJOR_CLASSES:
        ranked = list(
            df.sort_values(f"{class_name}_number", ascending=False)[["Gene_name"]]
            .drop_duplicates()["Gene_name"]
        )
        five = ranked[: go_rerun.CLASS_TOP_N]
        ten = ranked[:CLASS_TOP_N_10]
        assert five == ten[: len(five)], f"{class_name}: 5 % set is not a prefix of 10 %"
    print(f"nesting check: 5 % sets are strict prefixes of the 10 % sets "
          f"({go_rerun.CLASS_TOP_N:,} -> {CLASS_TOP_N_10:,} genes)")


def group_key(table: pd.DataFrame, level: str) -> pd.Series:
    """One label per GO study — class, or class + divergence direction."""
    column = GROUP_COLUMN[level]
    if "divergence_group" in table.columns:
        return table[column].astype(str) + " / " + table["divergence_group"].astype(str)
    return table[column].astype(str)


def compare(level: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Term-set agreement between the published 5 % arm and the new 10 % arm."""
    path5 = os.path.join(OUTPUT_DIR, f"GO_{level}_fdr01_reference.csv")
    path10 = os.path.join(OUTPUT_DIR, f"GO_{level}_p10_fdr01_reference.csv")
    if not (os.path.exists(path5) and os.path.exists(path10)):
        return pd.DataFrame(), pd.DataFrame()

    five = pd.read_csv(path5, low_memory=False)
    ten = pd.read_csv(path10, low_memory=False)
    five["group"] = group_key(five, level)
    ten["group"] = group_key(ten, level)

    # Only terms below the publication threshold enter the comparison: the 0.1
    # files exist purely so the retrieval cut is identical in both arms.
    five = five[five["FDR"] < PUBLISHED_FDR]
    ten = ten[ten["FDR"] < PUBLISHED_FDR]

    summary_rows, term_rows = [], []
    for group in sorted(set(five["group"]) | set(ten["group"])):
        g5 = five[five["group"] == group]
        g10 = ten[ten["group"] == group]

        # SVA and Helitron at class level are not ranked cuts, so they were not
        # re-run at 10 %. Their 5 % result IS their 10 % result; carrying it over is
        # what keeps them out of the "lost at 10 %" list they would otherwise fill.
        invariant = group.split(" / ")[0] in PERCENTILE_INVARIANT_GROUPS
        if invariant and g10.empty:
            g10 = g5

        t5, t10 = set(g5["Term ID"]), set(g10["Term ID"])
        shared = t5 & t10
        union = t5 | t10

        n5_fg = int(g5["n_foreground_genes"].iloc[0]) if len(g5) else 0
        n10_fg = int(g10["n_foreground_genes"].iloc[0]) if len(g10) else 0

        summary_rows.append(
            {
                "level": level,
                "group": group,
                "percentile_dependent": not invariant,
                "n_foreground_genes_5pct": n5_fg,
                "n_foreground_genes_10pct": n10_fg,
                "foreground_changed": n5_fg != n10_fg,
                "n_terms_5pct": len(t5),
                "n_terms_10pct": len(t10),
                "n_shared": len(shared),
                "n_lost_at_10pct": len(t5 - t10),
                "n_gained_at_10pct": len(t10 - t5),
                "jaccard": len(shared) / len(union) if union else float("nan"),
                "fraction_of_5pct_preserved": len(shared) / len(t5) if t5 else float("nan"),
            }
        )

        names = dict(zip(g5["Term ID"], g5["Term Name"]))
        names.update(dict(zip(g10["Term ID"], g10["Term Name"])))
        fdr5 = dict(zip(g5["Term ID"], g5["FDR"]))
        fdr10 = dict(zip(g10["Term ID"], g10["FDR"]))
        for term_id in sorted(t5 ^ t10):
            term_rows.append(
                {
                    "level": level,
                    "group": group,
                    "term_id": term_id,
                    "term_name": names.get(term_id, ""),
                    "status": "lost_at_10pct" if term_id in t5 else "gained_at_10pct",
                    "fdr_5pct": fdr5.get(term_id),
                    "fdr_10pct": fdr10.get(term_id),
                }
            )

    return pd.DataFrame(summary_rows), pd.DataFrame(term_rows)


def headline_table() -> pd.DataFrame:
    """Each abstract-level claim, with its FDR at 5 % and at 10 %."""
    rows = []
    for level, group, term in HEADLINE_CLAIMS:
        record = {"level": level, "group": group, "term_name": term}
        for arm, suffix in [("5pct", ""), ("10pct", "_p10")]:
            path = os.path.join(OUTPUT_DIR, f"GO_{level}{suffix}_fdr01_reference.csv")
            if not os.path.exists(path):
                record[f"fdr_{arm}"] = None
                record[f"survives_{arm}"] = None
                continue
            table = pd.read_csv(path, low_memory=False)
            table["group"] = group_key(table, level)
            # 'LINE' matches both 'LINE / highest' and 'LINE / lowest'; the claim is
            # about the class, so the strongest of its divergence groups is taken.
            hit = table[
                table["group"].str.startswith(group) & (table["Term Name"] == term)
            ]
            record[f"fdr_{arm}"] = float(hit["FDR"].min()) if len(hit) else None
            record[f"survives_{arm}"] = (
                bool(hit["FDR"].min() < PUBLISHED_FDR) if len(hit) else False
            )
        rows.append(record)

    table = pd.DataFrame(rows)
    table["verdict"] = [
        "survives both"
        if r.survives_5pct and r.survives_10pct
        else "5 % only"
        if r.survives_5pct
        else "10 % only"
        if r.survives_10pct
        else "neither"
        for r in table.itertuples()
    ]
    return table


def report() -> int:
    summary_path = os.path.join(OUTPUT_DIR, "percentile_sensitivity_summary.csv")
    if not os.path.exists(summary_path):
        print("Run without --summary first.", file=sys.stderr)
        return 1

    summary = pd.read_csv(summary_path)
    print("=" * 78)
    print(f"GO term stability, published top/bottom 5 % vs. {PERCENTILE} %, FDR < "
          f"{PUBLISHED_FDR}")
    print("=" * 78)
    for level in LEVELS:
        part = summary[summary["level"] == level]
        if part.empty:
            continue
        dependent = part[part["percentile_dependent"] & part["foreground_changed"]]
        print(f"\n{level}: {len(part)} groups, {len(dependent)} percentile-dependent")
        print(f"  terms at  5 %: {int(part['n_terms_5pct'].sum()):,}")
        print(f"  terms at 10 %: {int(part['n_terms_10pct'].sum()):,}")
        if len(dependent):
            print(f"  median Jaccard (percentile-dependent groups): "
                  f"{dependent['jaccard'].median():.3f}")
            print(f"  median fraction of 5 % terms preserved:       "
                  f"{dependent['fraction_of_5pct_preserved'].median():.3f}")
            worst = dependent.nsmallest(5, "jaccard")
            print("  least stable groups:")
            for r in worst.itertuples():
                print(f"    {r.group:24s} J = {r.jaccard:.3f}  "
                      f"({r.n_terms_5pct} -> {r.n_terms_10pct} terms, "
                      f"{r.n_lost_at_10pct} lost, {r.n_gained_at_10pct} gained)")

    headline_path = os.path.join(OUTPUT_DIR, "percentile_sensitivity_headline.csv")
    if os.path.exists(headline_path):
        print("\n" + "-" * 78)
        print("Headline claims at both percentiles:")
        head = pd.read_csv(headline_path)
        for r in head.itertuples():
            f5 = "n/a" if pd.isna(r.fdr_5pct) else f"{r.fdr_5pct:.4g}"
            f10 = "n/a" if pd.isna(r.fdr_10pct) else f"{r.fdr_10pct:.4g}"
            print(f"  {r.group:14s} {str(r.term_name)[:46]:48s} "
                  f"5 % {f5:>10s}   10 % {f10:>10s}   {r.verdict}")
        print("-" * 78)
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--level", choices=LEVELS + ["all"], default="all")
    parser.add_argument("--summary", action="store_true",
                        help="print the comparison from existing outputs and exit")
    parser.add_argument("--compare-only", action="store_true",
                        help="rebuild the comparison tables from GO output already on "
                             "disk, without re-running any enrichment")
    args = parser.parse_args()

    if args.summary:
        return report()

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    wanted = LEVELS if args.level == "all" else [args.level]

    if not args.compare_only:
        df = go_rerun.load_gene_table()
        assert_nesting(df)
        if "families" in wanted:
            print("Deriving per-family counts from TE_families ...")
            df = go_rerun.add_family_counts(df)

        for level in wanted:
            print(f"\n=== level: {level} at top/bottom {PERCENTILE} % ===")
            table = BUILDERS[level](df)
            reference_path = os.path.join(
                OUTPUT_DIR, f"GO_{level}_p10_fdr01_reference.csv")
            published_path = os.path.join(OUTPUT_DIR, f"GO_{level}_p10_fdr005.csv")
            table.to_csv(reference_path, index=False)
            table[table["FDR"] < PUBLISHED_FDR].to_csv(published_path, index=False)
            print(f"  wrote {os.path.basename(reference_path)} and "
                  f"{os.path.basename(published_path)}")

    summaries, term_changes = [], []
    for level in wanted:
        summary, terms = compare(level)
        if not summary.empty:
            summaries.append(summary)
        if not terms.empty:
            term_changes.append(terms)

    if summaries:
        pd.concat(summaries, ignore_index=True).to_csv(
            os.path.join(OUTPUT_DIR, "percentile_sensitivity_summary.csv"), index=False
        )
    if term_changes:
        pd.concat(term_changes, ignore_index=True).to_csv(
            os.path.join(OUTPUT_DIR, "percentile_sensitivity_terms.csv"), index=False
        )
    headline_table().to_csv(
        os.path.join(OUTPUT_DIR, "percentile_sensitivity_headline.csv"), index=False
    )
    print("\nWrote percentile_sensitivity_{summary,terms,headline}.csv")

    return report()


if __name__ == "__main__":
    sys.exit(main())
