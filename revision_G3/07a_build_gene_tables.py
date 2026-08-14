#!/usr/bin/env python
"""Build the per-window TE table for each TSS window size (the GO grid's input).

Why this script exists
----------------------
The published `../TEs_on_genes.csv` is the per-window TE table at 10 kb, and it is
what every GO gene-set construction in `06_go_rerun_fdr005.py` and
`05c_percentile_sensitivity.py` reads. Those constructions take the table as an
argument, so running GO at 5 kb or 20 kb needs nothing more than the same table
built at those widths — and that table does not exist. This script builds all
three from the canonical inputs.

Semantics, deliberately identical to the published table
-------------------------------------------------------
1. The key is the **TSS window**, not the gene. A gene with several annotated TSS
   contributes several windows, so an element within 10 kb of two TSS of the same
   gene is counted twice. That is the published design, flagged in the
   manuscript's Limitations, and reproducing it is what makes the 5 kb and 20 kb
   arms comparable with the published one rather than a different analysis. 214
   window intervals are shared by two genes, so the key is
   (Chromosome, Start, End, Gene_name), which is unique.
2. Windows with no intersecting TE are **retained with zero counts**. This matters:
   the GO background is `df["Gene_name"].unique()`, so dropping the 343 empty
   windows would shrink the background below 28,738 genes and make the FDRs
   incomparable between window sizes.
3. The comma-joined `TE_subfamilies` / `TE_families` / `TE_classes` /
   `Divergence_scores` columns are written in the same format
   `add_family_counts()` parses, with `.` for an empty window as in the published
   file.

What is NOT reproduced
----------------------
`Signal` and `Signals` in the published table come from the epigenomic bedGraph
the legacy pipeline carried alongside the TSS windows. No GO gene-set construction
reads them, and inventing empty columns for them would be worse than omitting
them, so they are absent here and the absence is asserted to be harmless: every
column the three builders touch is checked to be present.

The regression gate
-------------------
`--verify-10kb` rebuilds the 10 kb table and compares it against the published
`TEs_on_genes.csv`. Six invariants must hold exactly (row count, gene count, empty
windows, total elements, per-class totals, and `Average_Divergence_Score` on every
window to within 1e-9), plus the per-window multiset of family, class and subfamily
names. If those hold, the same code path at 5 kb and 20 kb is trustworthy; that is
the whole argument for the new window arms, so the script exits non-zero on any
failure rather than warning.

Usage
    python revision_G3/07a_build_gene_tables.py --verify-10kb   # gate only
    python revision_G3/07a_build_gene_tables.py                 # all three windows
    python revision_G3/07a_build_gene_tables.py --window 5kb
    python revision_G3/07a_build_gene_tables.py --summary
"""

from __future__ import annotations

import argparse
import io
import json
import os
import subprocess
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import revision_lib as rl  # noqa: E402

OUTPUT_DIR = rl.OUTPUT_DIR
WINDOWS = ["5kb", "10kb", "20kb"]

WINDOW_KEY = ["Chromosome", "Start", "End", "Gene_name"]

# Column order of the published table, minus the two epigenomic columns.
OUTPUT_COLUMNS = (
    WINDOW_KEY
    + ["Divergence_scores", "TE_subfamilies", "TE_families", "TE_classes",
       "Average_Divergence_Score"]
    + [f"{c}_number" for c in rl.CLASS_NAMES]
    + ["TE_number"]
    + [f"Divergence_Avg_{c}" for c in rl.CLASS_NAMES]
)

# Every column the three gene-set builders in 06/05c actually read.
COLUMNS_REQUIRED_BY_BUILDERS = (
    ["Gene_name", "TE_families", "TE_number", "Average_Divergence_Score"]
    + [f"{c}_number" for c in rl.CLASS_NAMES]
    + [f"Divergence_Avg_{c}" for c in ["LINE", "LTR", "SINE", "DNA"]]
)

# Published values, from ../CLAUDE.md and re-measured on TEs_on_genes.csv.
PUBLISHED_10KB = {
    "rows": 38_704,
    "genes": 28_738,
    "empty_windows": 343,
    "total_elements": 582_540,
    "per_class": {"LINE": 169_930, "LTR": 51_103, "SINE": 302_480,
                  "DNA": 57_684, "Retroposon": 1_170, "RC": 173},
}

# 05a's log. Printed as the cross-check that the other two windows ran the same
# code path — the totals are not independent knowledge, they are 05a's intersect.
EXPECTED_TOTAL_ELEMENTS = {"5kb": 293_652, "10kb": 582_540, "20kb": 1_157_235}

# The empty-window placeholder the published file uses in its joined columns.
EMPTY_JOINED = "."


def table_path(window: str) -> str:
    return os.path.join(OUTPUT_DIR, f"TEs_on_genes_{window}.csv")


def windows_path(window: str) -> str:
    return os.path.join(OUTPUT_DIR, f"windows_{window}.bed")


def read_windows(window: str) -> pd.DataFrame:
    """The full window list, which the aggregate is left-joined onto."""
    frame = pd.read_csv(
        windows_path(window), sep="\t", header=None, names=WINDOW_KEY,
        dtype={"Chromosome": str, "Start": np.int64, "End": np.int64,
               "Gene_name": str},
    )
    assert not frame.duplicated(WINDOW_KEY).any(), "window+gene key is not unique"
    return frame


def intersect_hits(window: str) -> pd.DataFrame:
    """One row per (element, window) overlap, as `bedtools intersect -wa -wb`.

    `-wa -wb` and not `-u`: the published table counts an element once per window
    it falls in, which is the multiple-TSS property documented above.
    """
    command = (
        f"bedtools intersect -a {rl.REPEATS_BED} -b {windows_path(window)} -wa -wb"
    )
    print(f"  {command}")
    start = time.time()
    raw = subprocess.run(
        command, shell=True, check=True, capture_output=True,
        env=dict(os.environ, LC_ALL="C"),
    ).stdout
    hits = pd.read_csv(
        io.BytesIO(raw), sep="\t", header=None,
        names=["te_chrom", "te_start", "te_end", "score", "subfamily_name",
               "family_name", "class_name",
               "Chromosome", "Start", "End", "Gene_name"],
        usecols=["score", "subfamily_name", "family_name", "class_name",
                 "Chromosome", "Start", "End", "Gene_name"],
        # A handful of subfamilies have an empty family_name in the RepeatMasker
        # annotation; without this pandas turns them into NaN and every later
        # groupby silently drops them (same trap as read_counts_file).
        keep_default_na=False,
        dtype={"score": np.int64, "subfamily_name": str, "family_name": str,
               "class_name": str, "Chromosome": str, "Start": np.int64,
               "End": np.int64, "Gene_name": str},
    )
    print(f"  {len(hits):,} TE-window intersections in {time.time() - start:.0f} s")
    return hits


def aggregate(hits: pd.DataFrame, windows: pd.DataFrame) -> pd.DataFrame:
    """Collapse the per-overlap rows into one row per window.

    The joined string columns keep the order `bedtools intersect` emitted, which
    is TE coordinate order because `-a` is the coordinate-sorted repeat file.
    """
    grouped = hits.groupby(WINDOW_KEY, sort=False)

    joined = pd.DataFrame({
        "Divergence_scores": grouped["score"].agg(
            lambda values: ",".join(str(int(v)) for v in values)),
        "TE_subfamilies": grouped["subfamily_name"].agg(",".join),
        "TE_families": grouped["family_name"].agg(",".join),
        "TE_classes": grouped["class_name"].agg(",".join),
        "Average_Divergence_Score": grouped["score"].mean(),
        "TE_number": grouped.size(),
    })

    per_class_count = (
        hits.groupby(WINDOW_KEY + ["class_name"], sort=False).size()
        .unstack("class_name", fill_value=0)
    )
    per_class_divergence = (
        hits.groupby(WINDOW_KEY + ["class_name"], sort=False)["score"].mean()
        .unstack("class_name")
    )
    for class_name in rl.CLASS_NAMES:
        if class_name not in per_class_count.columns:
            per_class_count[class_name] = 0
            per_class_divergence[class_name] = np.nan
    per_class_count = per_class_count[rl.CLASS_NAMES].add_suffix("_number")
    per_class_divergence = (
        per_class_divergence[rl.CLASS_NAMES]
        .rename(columns=lambda c: f"Divergence_Avg_{c}")
    )

    aggregated = joined.join(per_class_count).join(per_class_divergence)

    # Left join, so the empty windows survive with zeros rather than disappearing.
    table = windows.merge(aggregated.reset_index(), how="left", on=WINDOW_KEY)
    count_columns = [f"{c}_number" for c in rl.CLASS_NAMES] + ["TE_number"]
    table[count_columns] = table[count_columns].fillna(0).astype(np.int64)
    for column in ["Divergence_scores", "TE_subfamilies", "TE_families", "TE_classes"]:
        table[column] = table[column].fillna(EMPTY_JOINED)
    return table[OUTPUT_COLUMNS]


def build(window: str) -> pd.DataFrame:
    print(f"\n=== {window} ===")
    windows = read_windows(window)
    table = aggregate(intersect_hits(window), windows)

    missing = [c for c in COLUMNS_REQUIRED_BY_BUILDERS if c not in table.columns]
    assert not missing, f"missing columns the GO builders need: {missing}"

    total = int(table["TE_number"].sum())
    expected = EXPECTED_TOTAL_ELEMENTS[window]
    flag = "OK" if total == expected else "MISMATCH"
    print(f"  {len(table):,} windows, {table['Gene_name'].nunique():,} genes, "
          f"{int((table.TE_number == 0).sum()):,} empty, "
          f"{total:,} elements (05a expected {expected:,}: {flag})")
    return table


def write(window: str, table: pd.DataFrame) -> None:
    path = table_path(window)
    table.to_csv(path, index=False)
    print(f"  wrote {os.path.basename(path)} "
          f"({os.path.getsize(path) / 1e6:.1f} MB)")


def verify_10kb() -> int:
    """The regression gate: a rebuilt 10 kb table must be the published one."""
    print("=" * 78)
    print("07a regression gate — rebuilt 10 kb vs. the published TEs_on_genes.csv")
    print("=" * 78)

    built = build("10kb")
    published = pd.read_csv(
        os.path.join(rl.REPO_DIR, "TEs_on_genes.csv"), low_memory=False, index_col=0
    )

    failures = []

    def check(name, got, want):
        ok = got == want
        print(f"  {'OK  ' if ok else 'FAIL'} {name:34s} {got:>12,}  (expected {want:,})")
        if not ok:
            failures.append(name)

    check("rows", len(built), PUBLISHED_10KB["rows"])
    check("unique Gene_name", built["Gene_name"].nunique(), PUBLISHED_10KB["genes"])
    check("windows with TE_number == 0",
          int((built["TE_number"] == 0).sum()), PUBLISHED_10KB["empty_windows"])
    check("TE_number total", int(built["TE_number"].sum()),
          PUBLISHED_10KB["total_elements"])
    for class_name, want in PUBLISHED_10KB["per_class"].items():
        check(f"{class_name}_number total",
              int(built[f"{class_name}_number"].sum()), want)

    # Per-window comparison, aligned on the (window, gene) key rather than on row
    # order, so a different sort order is not reported as a value difference.
    left = built.set_index(WINDOW_KEY).sort_index()
    right = published.set_index(WINDOW_KEY).sort_index()
    if not left.index.equals(right.index):
        print("  FAIL window keys differ between the rebuilt and published tables")
        failures.append("window keys")
    else:
        diff = (left["Average_Divergence_Score"] - right["Average_Divergence_Score"])
        # Both are NaN on the 343 empty windows, so compare where either is set.
        both_nan = (left["Average_Divergence_Score"].isna()
                    & right["Average_Divergence_Score"].isna())
        worst = float(np.nanmax(np.abs(diff[~both_nan]))) if (~both_nan).any() else 0.0
        ok = worst < 1e-9 and bool(
            (left["Average_Divergence_Score"].isna()
             == right["Average_Divergence_Score"].isna()).all()
        )
        print(f"  {'OK  ' if ok else 'FAIL'} Average_Divergence_Score       "
              f"max |difference| = {worst:.3g}  (tolerance 1e-9)")
        if not ok:
            failures.append("Average_Divergence_Score")

        # Stronger than the totals: the same elements in every window, whatever
        # the order. This is what makes add_family_counts() give identical results.
        for column in ["TE_families", "TE_classes", "TE_subfamilies"]:
            def as_multiset(series):
                # NaN -> "" and not -> ".": three windows (RPSAP52, RBFOX3, PITX2)
                # contain a single element whose RepeatMasker `family_name` is
                # empty, so their joined string is empty in both tables and pandas
                # reads it back as NaN. Filling with "." would turn that read
                # artefact into a false difference.
                return series.fillna("").map(
                    lambda text: tuple(sorted(str(text).split(",")))
                )

            mismatched = int((as_multiset(left[column])
                              != as_multiset(right[column])).sum())
            ok = mismatched == 0
            print(f"  {'OK  ' if ok else 'FAIL'} {column:30s} "
                  f"{mismatched:,} windows differ as a multiset")
            if not ok:
                failures.append(column)

    print("-" * 78)
    if failures:
        print(f"GATE FAILED: {', '.join(failures)}")
        return 1

    write("10kb", built)
    print("GATE PASSED — the 10 kb rebuild reproduces the published table, so the "
          "same code path at 5 kb and 20 kb is comparable with it.")
    print(f"\nCross-check, from 05a's log: 5 kb {EXPECTED_TOTAL_ELEMENTS['5kb']:,} "
          f"and 20 kb {EXPECTED_TOTAL_ELEMENTS['20kb']:,} elements expected.")
    return 0


def summary() -> int:
    """Print what is on disk without recomputing anything."""
    print(f"{'window':8s} {'rows':>8s} {'genes':>8s} {'empty':>7s} {'elements':>12s}")
    found = 0
    for window in WINDOWS:
        path = table_path(window)
        if not os.path.exists(path):
            print(f"{window:8s} not built")
            continue
        found += 1
        table = pd.read_csv(path, low_memory=False)
        print(f"{window:8s} {len(table):>8,} {table.Gene_name.nunique():>8,} "
              f"{int((table.TE_number == 0).sum()):>7,} "
              f"{int(table.TE_number.sum()):>12,}")
    return 0 if found else 1


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--window", choices=WINDOWS,
                        help="build one window size instead of all three")
    parser.add_argument("--verify-10kb", action="store_true",
                        help="run the regression gate against the published table")
    parser.add_argument("--summary", action="store_true",
                        help="report the tables already on disk and exit")
    args = parser.parse_args()

    if args.summary:
        return summary()
    if args.verify_10kb:
        return verify_10kb()

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    wanted = [args.window] if args.window else WINDOWS

    # 10 kb always goes through the gate, so the other two are never written on
    # the strength of a code path that has not been shown to reproduce the paper.
    if "10kb" in wanted:
        if verify_10kb() != 0:
            return 1
        wanted = [w for w in wanted if w != "10kb"]
    elif not os.path.exists(table_path("10kb")):
        print("Run --verify-10kb first: the 5 kb / 20 kb tables are only "
              "trustworthy once the 10 kb rebuild has reproduced the published one.",
              file=sys.stderr)
        return 1

    manifest = {}
    for window in wanted:
        table = build(window)
        write(window, table)
        manifest[window] = {
            "rows": int(len(table)),
            "genes": int(table["Gene_name"].nunique()),
            "empty_windows": int((table["TE_number"] == 0).sum()),
            "total_elements": int(table["TE_number"].sum()),
            "per_class": {c: int(table[f"{c}_number"].sum()) for c in rl.CLASS_NAMES},
        }

    if manifest:
        path = os.path.join(OUTPUT_DIR, "gene_tables_manifest.json")
        existing = json.load(open(path)) if os.path.exists(path) else {}
        existing.update(manifest)
        with open(path, "w") as handle:
            json.dump(existing, handle, indent=2)
        print(f"\nwrote {os.path.basename(path)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
