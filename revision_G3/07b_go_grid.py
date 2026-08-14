#!/usr/bin/env python
"""GO enrichment across the full 3 windows x 2 percentiles grid, at FDR 0.05.

What was missing
----------------
GO existed at 10 kb / 5 % (published, `06_go_rerun_fdr005.py`) and 10 kb / 10 %
(`05c_percentile_sensitivity.py`). The other four cells of the 3 x 2 grid — 5 kb
and 20 kb at both percentiles — had never been run, so the robustness statement in
the Results was narrower than it read and the headline-claim table honestly said
`not evaluated` in its window columns. This script fills the grid.

Why it calls the existing builders instead of rebuilding the gene sets
---------------------------------------------------------------------
`06_go_rerun_fdr005.py`'s three `run_level_*(df)` functions and
`05c_percentile_sensitivity.py`'s three `build_*(df)` functions already take the
per-window gene table as an argument. So the new work only has to supply a
different `df` — one per window, from `07a_build_gene_tables.py` — and call them.
That is deliberate: a third copy of the gene-set construction would make a
construction difference indistinguishable from a window effect, which is precisely
the thing this grid exists to measure. Neither module is edited, and neither
module's own outputs are touched.

The two 10 kb cells are reused, not recomputed
----------------------------------------------
`--reuse-10kb` (the default) copies the existing published tables into the grid
naming rather than re-running them, so the published arm stays byte-stable and the
grid does not silently become a second, subtly different version of the paper.
`--check-reuse` proves the copies are row-identical to their sources and, unless
`--no-spot-check` is given, re-runs one level at 10 kb from the rebuilt gene table
and checks that it reproduces the published table term for term. That spot check is
what licenses reusing the other five.

Two things do not vary across the grid, and are asserted rather than assumed
--------------------------------------------------------------------------
* The gene count is 28,738 at every window (`07a` builds all three tables against
  the same 38,704-window list, empty windows retained), so the 1,436 / 2,872 gene
  cuts mean the same fraction in every cell.
* SVA and Helitron are every gene carrying any Retroposon / RC element, not a
  ranked cut, so they are percentile-invariant by construction at every window.
  They are carried through from the 5 % arm and flagged, never re-run as if the
  percentile changed them.

Outputs (all under output/GO_grid/)
    GO_<level>_<window>_p<pct>_fdr005.csv           18 files, what ships
    GO_<level>_<window>_p<pct>_fdr01_reference.csv  18 files, retrieval cut only
    INDEX.csv       one row per cell: level, window, percentile, counts, path
    MANIFEST.json   builders used, gene-table md5s, run time, reuse decisions

Usage
    python revision_G3/07b_go_grid.py                    # the four new cells
    python revision_G3/07b_go_grid.py --window 5kb
    python revision_G3/07b_go_grid.py --check-reuse
    python revision_G3/07b_go_grid.py --summary
"""

from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import os
import shutil
import sys
import time

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import revision_lib as rl  # noqa: E402

# Both module names start with a digit, so `import` cannot reach them. Calling
# their builders is what keeps all six grid cells built by one construction.
go_rerun = importlib.import_module("06_go_rerun_fdr005")
percentile = importlib.import_module("05c_percentile_sensitivity")

OUTPUT_DIR = rl.OUTPUT_DIR
GRID_DIR = os.path.join(OUTPUT_DIR, "GO_grid")

WINDOWS = ["5kb", "10kb", "20kb"]
PERCENTILES = [5, 10]
LEVELS = ["classes_count", "classes_divergence", "families"]
PUBLISHED_WINDOW = "10kb"
PUBLISHED_PERCENTILE = 5

REFERENCE_FDR = go_rerun.REFERENCE_FDR  # 0.10 retrieval cut
PUBLISHED_FDR = rl.FDR_THRESHOLD  # 0.05, what ships

# The 5 % arm is 06's builders; the 10 % arm is 05c's. Same constructions, one
# parameter apart — see each module's docstring for the frozen source cell.
BUILDERS = {
    (5, "classes_count"): go_rerun.run_level_classes_count,
    (5, "classes_divergence"): go_rerun.run_level_classes_divergence,
    (5, "families"): go_rerun.run_level_families,
    (10, "classes_count"): percentile.build_classes_count,
    (10, "classes_divergence"): percentile.build_classes_divergence,
    (10, "families"): percentile.build_families,
}

# Where the two published 10 kb cells are copied from.
REUSE_SOURCES = {
    (10, "10kb", "classes_count"): "GO_classes_count_p10",
    (10, "10kb", "classes_divergence"): "GO_classes_divergence_p10",
    (10, "10kb", "families"): "GO_families_p10",
    (5, "10kb", "classes_count"): "GO_classes_count",
    (5, "10kb", "classes_divergence"): "GO_classes_divergence",
    (5, "10kb", "families"): "GO_families",
}

EXPECTED_GENES = 28_738

# Groups with no percentile parameter (frozen GO cell 20), same set 05c uses.
PERCENTILE_INVARIANT_GROUPS = percentile.PERCENTILE_INVARIANT_GROUPS

GROUP_COLUMN = {
    "classes_count": "class_name",
    "classes_divergence": "class_name",
    "families": "family_name",
}


def cell_stem(level: str, window: str, pct: int) -> str:
    return f"GO_{level}_{window}_p{pct}"


def published_path(level: str, window: str, pct: int) -> str:
    return os.path.join(GRID_DIR, f"{cell_stem(level, window, pct)}_fdr005.csv")


def reference_path(level: str, window: str, pct: int) -> str:
    return os.path.join(GRID_DIR, f"{cell_stem(level, window, pct)}_fdr01_reference.csv")


def gene_table_path(window: str) -> str:
    return os.path.join(OUTPUT_DIR, f"TEs_on_genes_{window}.csv")


def md5(path: str) -> str:
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_window_table(window: str, with_families: bool) -> pd.DataFrame:
    """The per-window gene table from 07a, with the family counts if needed."""
    path = gene_table_path(window)
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"{os.path.basename(path)} is missing. Run "
            f"07a_build_gene_tables.py --verify-10kb first, then the other windows."
        )
    df = pd.read_csv(path, low_memory=False)
    n_genes = df["Gene_name"].nunique()
    assert n_genes == EXPECTED_GENES, (
        f"{window}: {n_genes:,} genes, expected {EXPECTED_GENES:,} — the 1,436 / "
        f"2,872 gene cuts would not be the same fraction across the grid"
    )
    print(f"  {os.path.basename(path)}: {len(df):,} windows, {n_genes:,} genes")
    if with_families:
        print("  deriving per-family counts from TE_families ...")
        df = go_rerun.add_family_counts(df)
    return df


def write_cell(level: str, window: str, pct: int, table: pd.DataFrame) -> dict:
    table.to_csv(reference_path(level, window, pct), index=False)
    ships = table[table["FDR"] < PUBLISHED_FDR].copy()
    ships.to_csv(published_path(level, window, pct), index=False)
    return summarise_cell(level, window, pct, table, ships, source="computed")


def summarise_cell(level, window, pct, reference, ships, source) -> dict:
    group_column = GROUP_COLUMN[level]
    n_foreground = (
        int(reference["n_foreground_genes"].max())
        if "n_foreground_genes" in reference.columns and len(reference) else 0
    )
    record = {
        "level": level,
        "window": window,
        "percentile": pct,
        "is_published_cell": window == PUBLISHED_WINDOW and pct == PUBLISHED_PERCENTILE,
        "source": source,
        "n_groups": int(reference[group_column].nunique()) if len(reference) else 0,
        "foreground_cut_size": n_foreground,
        "n_terms_01": int(len(reference)),
        "n_terms_005": int(len(ships)),
        "n_groups_with_terms": int(ships[group_column].nunique()) if len(ships) else 0,
        "path": os.path.relpath(published_path(level, window, pct), rl.REVISION_DIR),
    }
    print(f"  {level:20s} {window:5s} {pct:>3}%  {record['n_terms_01']:>5} terms at "
          f"FDR<{REFERENCE_FDR} -> {record['n_terms_005']:>5} at FDR<{PUBLISHED_FDR} "
          f"({record['n_groups_with_terms']}/{record['n_groups']} groups, {source})")
    return record


def reuse_10kb(percentiles: list[int], levels: list[str]) -> list[dict]:
    """Copy the two existing 10 kb cells into the grid naming, unchanged."""
    print("\n=== reusing the two published 10 kb cells (not recomputed) ===")
    records = []
    for (pct, window, level), stem in sorted(REUSE_SOURCES.items()):
        if pct not in percentiles or level not in levels:
            continue
        for suffix in ["fdr005", "fdr01_reference"]:
            source = os.path.join(OUTPUT_DIR, f"{stem}_{suffix}.csv")
            if not os.path.exists(source):
                raise FileNotFoundError(
                    f"{os.path.basename(source)} is missing — run "
                    f"06_go_rerun_fdr005.py --level all and "
                    f"05c_percentile_sensitivity.py first."
                )
            target = os.path.join(
                GRID_DIR, f"{cell_stem(level, window, pct)}_{suffix}.csv"
            )
            shutil.copyfile(source, target)
        reference = pd.read_csv(reference_path(level, window, pct), low_memory=False)
        ships = pd.read_csv(published_path(level, window, pct), low_memory=False)
        records.append(
            summarise_cell(level, window, pct, reference, ships,
                           source=f"reused from {stem}_*.csv")
        )
    return records


def compute_cells(windows: list[str], percentiles: list[int],
                  levels: list[str]) -> list[dict]:
    """Run the grid cells that are not the reused 10 kb ones."""
    records = []
    for window in windows:
        print(f"\n=== window {window} ===")
        df = load_window_table(window, with_families="families" in levels)
        for pct in percentiles:
            for level in levels:
                print(f"\n--- {level} at {window} / top-bottom {pct} % ---")
                start = time.time()
                table = BUILDERS[(pct, level)](df)
                records.append(write_cell(level, window, pct, table))
                print(f"  cell finished in {time.time() - start:.0f} s")
    return records


def write_index(records: list[dict]) -> None:
    path = os.path.join(GRID_DIR, "INDEX.csv")
    index = pd.DataFrame(records)
    if os.path.exists(path):
        previous = pd.read_csv(path)
        index = (
            pd.concat([previous, index], ignore_index=True)
            .drop_duplicates(["level", "window", "percentile"], keep="last")
        )
    index = index.sort_values(["level", "window", "percentile"])
    index.to_csv(path, index=False)
    print(f"\nwrote {os.path.relpath(path, rl.REVISION_DIR)} ({len(index)} cells)")


def write_manifest(records: list[dict]) -> None:
    path = os.path.join(GRID_DIR, "MANIFEST.json")
    manifest = json.load(open(path)) if os.path.exists(path) else {}
    manifest.update(
        {
            "script": os.path.basename(__file__),
            "created": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "retrieval_fdr": REFERENCE_FDR,
            "published_fdr": PUBLISHED_FDR,
            "builders": {
                f"{pct}pct/{level}": f"{fn.__module__}.{fn.__name__}"
                for (pct, level), fn in BUILDERS.items()
            },
            "class_top_n": {"5pct": go_rerun.CLASS_TOP_N,
                            "10pct": percentile.CLASS_TOP_N_10},
            "percentile_invariant_groups": sorted(PERCENTILE_INVARIANT_GROUPS),
            "gene_table_md5": {
                window: md5(gene_table_path(window))
                for window in WINDOWS if os.path.exists(gene_table_path(window))
            },
        }
    )
    cells = manifest.setdefault("cells", {})
    for record in records:
        cells[f"{record['level']}/{record['window']}/p{record['percentile']}"] = {
            "source": record["source"],
            "n_terms_005": record["n_terms_005"],
            "n_terms_01": record["n_terms_01"],
        }
    with open(path, "w") as handle:
        json.dump(manifest, handle, indent=2)
    print(f"wrote {os.path.relpath(path, rl.REVISION_DIR)}")


def check_reuse(spot_check_level: str | None) -> int:
    """Prove the copied 10 kb cells are the published ones, term for term."""
    print("=" * 78)
    print("07b reuse check — the six 10 kb grid files against their sources")
    print("=" * 78)
    failures = []
    for (pct, window, level), stem in sorted(REUSE_SOURCES.items()):
        target = published_path(level, window, pct)
        source = os.path.join(OUTPUT_DIR, f"{stem}_fdr005.csv")
        if not os.path.exists(target):
            print(f"  FAIL {level:20s} p{pct:<3} grid file not written yet")
            failures.append(f"{level}/p{pct}")
            continue
        a = pd.read_csv(source, low_memory=False)
        b = pd.read_csv(target, low_memory=False)
        identical = a.equals(b)
        print(f"  {'OK  ' if identical else 'FAIL'} {level:20s} p{pct:<3} "
              f"{len(b):>5} rows  {'row-identical' if identical else 'DIFFERS'} "
              f"to {stem}_fdr005.csv")
        if not identical:
            failures.append(f"{level}/p{pct}")

    if spot_check_level:
        print(f"\nspot check: re-running {spot_check_level} at 10 kb / 5 % from the "
              f"rebuilt gene table")
        df = load_window_table(PUBLISHED_WINDOW,
                               with_families=spot_check_level == "families")
        rerun = BUILDERS[(5, spot_check_level)](df)
        rerun = rerun[rerun["FDR"] < PUBLISHED_FDR]
        published = pd.read_csv(
            os.path.join(OUTPUT_DIR, f"GO_{spot_check_level}_fdr005.csv"),
            low_memory=False,
        )
        group_column = GROUP_COLUMN[spot_check_level]
        key_columns = [group_column, "Term ID"]
        if "divergence_group" in published.columns:
            key_columns.append("divergence_group")
        left = set(map(tuple, rerun[key_columns].astype(str).to_numpy()))
        right = set(map(tuple, published[key_columns].astype(str).to_numpy()))
        identical = left == right
        print(f"  {'OK  ' if identical else 'FAIL'} {len(left)} (group, term) pairs "
              f"re-run vs {len(right)} published; "
              f"{len(left ^ right)} differ")
        if not identical:
            for item in sorted(left ^ right)[:10]:
                print(f"    only in {'re-run' if item in left else 'published'}: {item}")
            failures.append(f"spot check {spot_check_level}")

    print("-" * 78)
    if failures:
        print(f"REUSE CHECK FAILED: {', '.join(failures)}")
        return 1
    print("REUSE CHECK PASSED — the 10 kb grid cells are the published tables.")
    return 0


def summary() -> int:
    path = os.path.join(GRID_DIR, "INDEX.csv")
    if not os.path.exists(path):
        print("No grid yet — run without --summary first.", file=sys.stderr)
        return 1
    index = pd.read_csv(path)
    print("=" * 78)
    print(f"GO grid: {len(index)} of {len(WINDOWS) * len(PERCENTILES) * len(LEVELS)} "
          f"cells, terms at FDR < {PUBLISHED_FDR}")
    print("=" * 78)
    for level in LEVELS:
        part = index[index.level == level]
        if part.empty:
            continue
        print(f"\n{level}")
        header = "  " + "".join(f"{w:>12s}" for w in WINDOWS)
        print(f"  {'':8s}{header}")
        for pct in PERCENTILES:
            cells = []
            for window in WINDOWS:
                hit = part[(part.window == window) & (part.percentile == pct)]
                cells.append(f"{int(hit.n_terms_005.iloc[0]):,}" if len(hit) else "-")
            marker = "  (published)" if pct == PUBLISHED_PERCENTILE else ""
            print(f"  {pct:>3}% " + "".join(f"{c:>12s}" for c in cells) + marker)
    missing = [
        (level, window, pct)
        for level in LEVELS for window in WINDOWS for pct in PERCENTILES
        if index[(index.level == level) & (index.window == window)
                 & (index.percentile == pct)].empty
    ]
    print(f"\nmissing cells: {missing if missing else 'none'}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--window", choices=WINDOWS, action="append",
                        help="restrict to one window (repeatable)")
    parser.add_argument("--percentile", type=int, choices=PERCENTILES, action="append",
                        help="restrict to one percentile arm (repeatable)")
    parser.add_argument("--level", choices=LEVELS, action="append",
                        help="restrict to one GO level (repeatable)")
    parser.add_argument("--no-reuse-10kb", action="store_true",
                        help="recompute the 10 kb cells instead of copying them "
                             "(this changes the published arm — do not use for the "
                             "submission grid)")
    parser.add_argument("--check-reuse", action="store_true",
                        help="prove the copied 10 kb cells are the published tables")
    parser.add_argument("--no-spot-check", action="store_true",
                        help="with --check-reuse, skip the re-run comparison")
    parser.add_argument("--spot-check-level", choices=LEVELS, default="classes_count")
    parser.add_argument("--summary", action="store_true")
    args = parser.parse_args()

    if args.summary:
        return summary()

    os.makedirs(GRID_DIR, exist_ok=True)

    if args.check_reuse:
        return check_reuse(None if args.no_spot_check else args.spot_check_level)

    windows = args.window or WINDOWS
    percentiles = args.percentile or PERCENTILES
    levels = args.level or LEVELS

    records = []
    if PUBLISHED_WINDOW in windows and not args.no_reuse_10kb:
        records += reuse_10kb(percentiles, levels)
    records += compute_cells(
        [w for w in windows if not (w == PUBLISHED_WINDOW and not args.no_reuse_10kb)],
        percentiles, levels,
    )

    if records:
        write_index(records)
        write_manifest(records)
    return summary()


if __name__ == "__main__":
    sys.exit(main())
