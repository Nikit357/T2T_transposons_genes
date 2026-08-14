#!/usr/bin/env python
"""WP1b — reconstruct a legacy permutation BED from the compacted counts store.

This is what makes the compaction defensible: the legacy format is never truly
lost. Anyone who needs to re-run the frozen `TEs_mapped_on_TSS_analysis.ipynb`,
which reads `consolidated_random_data.csv`, can rebuild exactly what it expects.

What "exact" means here
-----------------------
The source BEDs were an *unordered multiset* — the row order bedtools happened to
emit. That order was never recorded and carries no information, so the
reconstruction is exact **as a multiset**, which is what every downstream
statistic sees. Compare with `sort` on both sides:

    python revision_G3/01c_expand_counts.py --seed 1 | sort > /tmp/rebuilt.bed
    sort ../epigenomic_files/permutation_results/\\
repeats_intersected_with_TSS_random_1.bed > /tmp/orig.bed
    diff -q /tmp/orig.bed /tmp/rebuilt.bed

`--check` does exactly that comparison for you, when the source still exists.

Usage
-----
    python revision_G3/01c_expand_counts.py --seed 1               # to stdout
    python revision_G3/01c_expand_counts.py --seed 1 --check       # prove it
    python revision_G3/01c_expand_counts.py --seed 1 -o out.bed
    python revision_G3/01c_expand_counts.py --all-seeds --consolidated out.csv
"""

from __future__ import annotations

import argparse
import os
import shlex
import subprocess
import sys

REPO_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LEGACY_DIR = os.path.abspath(
    os.path.join(REPO_DIR, "..", "epigenomic_files", "permutation_results")
)

# `sorted` is required for round-trip proofs; `n` is a weight, so each tuple is
# emitted `n` times. Printing $2..$5 drops the seed column, giving the exact
# 4-column legacy BED layout (score, subfamily_name, family_name, class_name).
AWK_EXPAND = r"{ for (i = 0; i < $6; i++) print $2, $3, $4, $5 }"

# The consolidated CSV is tab-separated despite its extension, and carries a
# header plus permutation_id as the first column.
CONSOLIDATED_HEADER = "permutation_id\tscore\tsubfamily_name\tfamily_name\tclass_name"
AWK_EXPAND_WITH_SEED = r"{ for (i = 0; i < $6; i++) print $1, $2, $3, $4, $5 }"


def store_dir(window: str) -> str:
    return os.path.join(REPO_DIR, "revision_G3", "output", f"permutation_counts_{window}")


def counts_path(seed: int, window: str) -> str:
    return os.path.join(store_dir(window), f"counts_seed_{seed}.tsv.zst")


def run(cmd: str, stdout=None) -> None:
    env = dict(os.environ, LC_ALL="C")
    subprocess.run(cmd, shell=True, check=True, env=env, stdout=stdout)


def capture(cmd: str) -> str:
    env = dict(os.environ, LC_ALL="C")
    return subprocess.run(
        cmd, shell=True, check=True, capture_output=True, text=True, env=env
    ).stdout


def expand_command(seed: int, window: str, with_seed: bool = False) -> str:
    path = counts_path(seed, window)
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    prog = AWK_EXPAND_WITH_SEED if with_seed else AWK_EXPAND
    return (
        f"zstd -dc {shlex.quote(path)} "
        f"| awk -F'\\t' -v OFS='\\t' {shlex.quote(prog)}"
    )


def check_seed(seed: int, window: str) -> bool:
    """Byte-compare the sorted reconstruction against the sorted source BED."""
    src = os.path.join(LEGACY_DIR, f"repeats_intersected_with_TSS_random_{seed}.bed")
    if not os.path.exists(src):
        print(
            f"Cannot --check seed {seed}: the source BED is gone "
            f"(that is expected after compaction). The MANIFEST.json entry holds "
            f"`source_sorted_md5`, which 01b verified before deleting it.",
            file=sys.stderr,
        )
        manifest_md5 = manifest_source_md5(seed, window)
        if manifest_md5 is None:
            return False
        rebuilt = capture(f"{expand_command(seed, window)} | sort | md5sum").split()[0]
        ok = rebuilt == manifest_md5
        print(
            f"seed {seed}: reconstruction md5 {rebuilt} "
            f"{'==' if ok else '!='} manifest {manifest_md5}"
        )
        return ok

    src_md5 = capture(f"sort {shlex.quote(src)} | md5sum").split()[0]
    rebuilt_md5 = capture(f"{expand_command(seed, window)} | sort | md5sum").split()[0]
    ok = src_md5 == rebuilt_md5
    n_src = int(capture(f"wc -l < {shlex.quote(src)}").strip())
    n_new = int(capture(f"{expand_command(seed, window)} | wc -l").strip())
    print(f"seed {seed}: {n_src:,} rows in source, {n_new:,} reconstructed")
    print(f"seed {seed}: sorted md5 {src_md5} {'==' if ok else '!='} {rebuilt_md5}")
    print("seed %d reconstructs EXACTLY" % seed if ok else "MISMATCH")
    return ok


def manifest_source_md5(seed: int, window: str) -> str | None:
    import json

    path = os.path.join(store_dir(window), "MANIFEST.json")
    if not os.path.exists(path):
        return None
    with open(path) as handle:
        manifest = json.load(handle)
    record = manifest.get("seeds", {}).get(str(seed))
    return record.get("source_sorted_md5") if record else None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--seed", type=int, help="permutation index to expand")
    parser.add_argument(
        "--all-seeds", action="store_true", help="expand every seed in the store"
    )
    parser.add_argument("--window", default="10kb", choices=["5kb", "10kb", "20kb"])
    parser.add_argument("-o", "--output", help="write here instead of stdout")
    parser.add_argument(
        "--consolidated",
        help="rebuild a `consolidated_random_data.csv` equivalent at this path",
    )
    parser.add_argument(
        "--sorted",
        action="store_true",
        help="emit in sorted order (required for md5 comparison against a source)",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="prove the reconstruction matches the source (or the manifest md5)",
    )
    args = parser.parse_args()

    if args.check:
        if not args.seed:
            parser.error("--check needs --seed")
        return 0 if check_seed(args.seed, args.window) else 1

    if args.consolidated:
        import glob

        paths = sorted(
            glob.glob(os.path.join(store_dir(args.window), "counts_seed_*.tsv.zst")),
            key=lambda p: int(p.rsplit("_", 1)[1].split(".")[0]),
        )
        if not paths:
            parser.error(f"no counts files in {store_dir(args.window)}")
        with open(args.consolidated, "w") as handle:
            handle.write(CONSOLIDATED_HEADER + "\n")
            handle.flush()
            for path in paths:
                seed = int(path.rsplit("_", 1)[1].split(".")[0])
                run(expand_command(seed, args.window, with_seed=True), stdout=handle)
        print(
            f"Rebuilt {args.consolidated} "
            f"({os.path.getsize(args.consolidated):,} B) from {len(paths)} seeds",
            file=sys.stderr,
        )
        return 0

    if args.all_seeds:
        parser.error("--all-seeds is only meaningful with --consolidated")
    if not args.seed:
        parser.error("give --seed N (or --consolidated PATH)")

    cmd = expand_command(args.seed, args.window)
    if args.sorted:
        cmd += " | sort"
    if args.output:
        with open(args.output, "w") as handle:
            run(cmd, stdout=handle)
        print(
            f"Wrote {args.output} ({os.path.getsize(args.output):,} B)", file=sys.stderr
        )
    else:
        run(cmd)
    return 0


if __name__ == "__main__":
    sys.exit(main())
