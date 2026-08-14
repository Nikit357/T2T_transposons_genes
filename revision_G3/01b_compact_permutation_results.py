#!/usr/bin/env python
"""WP1b — losslessly compact `permutation_results/` from 11 GB to ~110 MB.

Why this is lossless
--------------------
Each per-seed file `repeats_intersected_with_TSS_random_${SEED}.bed` is an
*unordered multiset* of 4-tuples `(score, subfamily_name, family_name,
class_name)` — one row per intersecting TE, in the arbitrary order bedtools
emitted. Row order carries no information, so run-length encoding by tuple
loses nothing that any downstream statistic (count, mean, SD, quantile,
distribution) can detect. Measured on seed 1: 545,099 rows collapse to 70,040
unique tuples, 10,729,853 B -> 217,587 B with `zstd -19` (49x). Across 500 seeds
that is ~109 MB, and `consolidated_random_data.csv` (6.37 GB) becomes redundant
because it is only the concatenation of the per-seed files with a
`permutation_id` column, which the compact store carries as `seed`.

Output schema (identical to what `01_permutations_stream.sh` produces for the
new 5 kb / 20 kb windows, so compaction and the new runs share one reader):

    seed <TAB> score <TAB> subfamily_name <TAB> family_name <TAB> class_name <TAB> n

Two data facts that dictate the implementation
----------------------------------------------
1. Rows are TAB-separated with exactly 4 fields, but 2,730 rows per seed carry
   an EMPTY `family_name`. Every awk invocation therefore sets `-F'\\t'` and
   `OFS='\\t'`; a default-FS awk would silently see those rows as 3 fields and
   shift `class_name` into `family_name`.
2. All sorting is done under `LC_ALL=C` so the byte-exactness comparison between
   the source and the reconstruction is not at the mercy of the locale.

Safety
------
Nothing is deleted unless `--delete-verified` is passed, and then only after that
seed's own verification has passed: total row count, per-class totals, and a
byte-exact md5 comparison of `sort <source>` against `sort <reconstruction>`.
Peak extra disk use is one seed (~12 MB), not the whole store. If any check
fails the script stops and the source is kept (caveat C2 — never delete a source
whose verification did not pass).

Usage
-----
    # smoke test on 5 seeds, nothing deleted
    python revision_G3/01b_compact_permutation_results.py --limit 5

    # full run, freeing disk as it goes
    python revision_G3/01b_compact_permutation_results.py --delete-verified

    # re-verify an existing store against whatever sources remain
    python revision_G3/01b_compact_permutation_results.py --verify-only

    # aggregate check against the 6.37 GB CSV, then allow its deletion
    python revision_G3/01b_compact_permutation_results.py --aggregate-check
    python revision_G3/01b_compact_permutation_results.py --delete-consolidated
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import shlex
import subprocess
import sys

SCRIPT_VERSION = "1.0.0"

REPO_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LEGACY_DIR = os.path.abspath(
    os.path.join(REPO_DIR, "..", "epigenomic_files", "permutation_results")
)
STORE_DIR = os.path.join(REPO_DIR, "revision_G3", "output", "permutation_counts_10kb")
CONSOLIDATED_CSV = os.path.join(LEGACY_DIR, "consolidated_random_data.csv")
MANIFEST_PATH = os.path.join(STORE_DIR, "MANIFEST.json")

ZSTD_LEVEL = 19
N_SEEDS = 500

# One awk program, reused, that run-length-encodes a per-seed BED. `$0` is
# printed verbatim so an empty family_name survives untouched.
AWK_ENCODE = (
    r"{ c[$0]++ } END { for (k in c) print seed, k, c[k] }"
)
# The exact inverse: expand `n` back into `n` identical 4-column rows.
AWK_EXPAND = r"{ for (i = 0; i < $6; i++) print $2, $3, $4, $5 }"


def sh(cmd: str, **kwargs) -> str:
    """Run a shell pipeline under LC_ALL=C and return stdout, raising on error."""
    env = dict(os.environ, LC_ALL="C")
    result = subprocess.run(
        cmd, shell=True, check=True, capture_output=True, text=True, env=env, **kwargs
    )
    return result.stdout


def source_bed(seed: int) -> str:
    return os.path.join(LEGACY_DIR, f"repeats_intersected_with_TSS_random_{seed}.bed")


def counts_zst(seed: int) -> str:
    return os.path.join(STORE_DIR, f"counts_seed_{seed}.tsv.zst")


def md5_of_sorted(path: str) -> str:
    """md5 of the file's lines in LC_ALL=C sorted order."""
    out = sh(f"sort {shlex.quote(path)} | md5sum")
    return out.split()[0]


def md5_of_sorted_expansion(zst_path: str) -> str:
    """md5 of the sorted 4-column BED reconstructed from a counts file."""
    out = sh(
        f"zstd -dc {shlex.quote(zst_path)} "
        f"| awk -F'\\t' -v OFS='\\t' {shlex.quote(AWK_EXPAND)} "
        f"| sort | md5sum"
    )
    return out.split()[0]


def per_class_totals(path: str, is_counts: bool) -> dict[str, int]:
    """Rows per class_name, from either a source BED or a counts file."""
    if is_counts:
        prog = r"{ t[$5] += $6 } END { for (k in t) print k, t[k] }"
        cmd = (
            f"zstd -dc {shlex.quote(path)} "
            f"| awk -F'\\t' -v OFS='\\t' {shlex.quote(prog)}"
        )
    else:
        prog = r"{ t[$4] += 1 } END { for (k in t) print k, t[k] }"
        cmd = f"awk -F'\\t' -v OFS='\\t' {shlex.quote(prog)} {shlex.quote(path)}"
    totals = {}
    for line in sh(cmd).splitlines():
        name, count = line.rsplit("\t", 1)
        totals[name] = int(count)
    return totals


def file_md5(path: str, chunk: int = 1 << 20) -> str:
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(chunk), b""):
            digest.update(block)
    return digest.hexdigest()


def free_bytes(path: str) -> int:
    stat = os.statvfs(path)
    return stat.f_bavail * stat.f_frsize


def encode_seed(seed: int) -> tuple[str, int, int]:
    """Write `counts_seed_${seed}.tsv.zst`; return (path, n_tuples, size)."""
    src = source_bed(seed)
    tmp = os.path.join(STORE_DIR, f".counts_seed_{seed}.tsv.tmp")
    out = counts_zst(seed)

    # Sorting the encoded output keeps the store deterministic, which is what
    # makes the per-file md5s in MANIFEST.json meaningful as a provenance record.
    # A plain whole-line LC_ALL=C sort is used rather than `sort -t$'\t' -k...`:
    # it is equally deterministic and avoids passing a literal tab through a
    # shell that may be dash rather than bash.
    sh(
        f"awk -F'\\t' -v OFS='\\t' -v seed={seed} {shlex.quote(AWK_ENCODE)} "
        f"{shlex.quote(src)} | sort > {shlex.quote(tmp)}"
    )
    n_tuples = int(sh(f"wc -l < {shlex.quote(tmp)}").strip())
    sh(f"zstd -q -{ZSTD_LEVEL} -f {shlex.quote(tmp)} -o {shlex.quote(out)}")

    # The compressed file must decompress to exactly what we just wrote.
    plain_md5 = sh(f"md5sum < {shlex.quote(tmp)}").split()[0]
    round_trip = sh(f"zstd -dc {shlex.quote(out)} | md5sum").split()[0]
    if plain_md5 != round_trip:
        raise RuntimeError(f"seed {seed}: zstd round-trip mismatch")
    os.remove(tmp)
    return out, n_tuples, os.path.getsize(out)


def verify_seed(seed: int) -> dict:
    """Prove the counts file is a lossless encoding of the source BED.

    Three independent checks, all of which must pass before the source may be
    deleted: the total row count, the per-class row totals, and a byte-exact md5
    of the sorted source against the sorted reconstruction. The third subsumes
    the first two, but a mismatch in one of the cheap checks localises the fault,
    which is worth the few seconds.
    """
    src = source_bed(seed)
    zst = counts_zst(seed)

    n_src = int(sh(f"wc -l < {shlex.quote(src)}").strip())
    n_encoded = int(
        sh(
            f"zstd -dc {shlex.quote(zst)} | awk -F'\\t' "
            f"'{{ s += $6 }} END {{ print s+0 }}'"
        ).strip()
    )
    if n_src != n_encoded:
        raise RuntimeError(
            f"seed {seed}: row count {n_encoded} from counts != {n_src} in source"
        )

    totals_src = per_class_totals(src, is_counts=False)
    totals_enc = per_class_totals(zst, is_counts=True)
    if totals_src != totals_enc:
        raise RuntimeError(
            f"seed {seed}: per-class totals differ\n"
            f"  source: {sorted(totals_src.items())}\n"
            f"  counts: {sorted(totals_enc.items())}"
        )

    md5_src = md5_of_sorted(src)
    md5_rebuilt = md5_of_sorted_expansion(zst)
    if md5_src != md5_rebuilt:
        raise RuntimeError(
            f"seed {seed}: multiset md5 mismatch ({md5_src} vs {md5_rebuilt})"
        )

    return {
        "seed": seed,
        "source_file": os.path.basename(src),
        "source_bytes": os.path.getsize(src),
        "source_rows": n_src,
        "source_sorted_md5": md5_src,
        "counts_file": os.path.basename(zst),
        "counts_bytes": os.path.getsize(zst),
        "counts_tuples": int(sh(f"zstd -dc {shlex.quote(zst)} | wc -l").strip()),
        "counts_md5": file_md5(zst),
        "per_class_rows": totals_src,
        "verified": True,
    }


def load_manifest() -> dict:
    if os.path.exists(MANIFEST_PATH):
        with open(MANIFEST_PATH) as handle:
            return json.load(handle)
    return {
        "script": os.path.basename(__file__),
        "script_version": SCRIPT_VERSION,
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "purpose": (
            "Lossless run-length encoding of the N=500 permutation background for "
            "the 10 kb TSS window (plan WP1b). Row order in the source BEDs is "
            "meaningless, so encoding by unique (score, subfamily, family, class) "
            "tuple preserves every downstream statistic exactly."
        ),
        "source_directory": LEGACY_DIR,
        "window": "10kb",
        "n_permutations": N_SEEDS,
        "schema": [
            "seed",
            "score",
            "subfamily_name",
            "family_name",
            "class_name",
            "n",
        ],
        "schema_notes": (
            "Tab-separated, no header. `n` is a WEIGHT: the number of intersecting "
            "elements with that exact tuple in that permutation. `family_name` is "
            "empty for a small number of subfamilies (2,730 rows in seed 1); this "
            "is faithful to the source and must be read with explicit tab "
            "separation and NA-keeping off."
        ),
        "compression": f"zstd -{ZSTD_LEVEL}",
        "reconstruct_with": "revision_G3/01c_expand_counts.py --seed N",
        "seeds": {},
        "aggregate_check": None,
        "consolidated_csv_deleted": False,
    }


def save_manifest(manifest: dict) -> None:
    manifest["updated_utc"] = dt.datetime.now(dt.timezone.utc).isoformat(
        timespec="seconds"
    )
    done = manifest["seeds"]
    manifest["n_seeds_compacted"] = len(done)
    manifest["total_source_bytes"] = sum(v["source_bytes"] for v in done.values())
    manifest["total_counts_bytes"] = sum(v["counts_bytes"] for v in done.values())
    if manifest["total_counts_bytes"]:
        manifest["compression_ratio"] = round(
            manifest["total_source_bytes"] / manifest["total_counts_bytes"], 2
        )
    tmp = MANIFEST_PATH + ".tmp"
    with open(tmp, "w") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=False)
    os.replace(tmp, MANIFEST_PATH)


def aggregate_check(manifest: dict) -> dict:
    """Compare the compacted store against `consolidated_random_data.csv`.

    Streams the 6.37 GB CSV once, accumulating per (permutation_id, class_name)
    and per (permutation_id, family_name) row totals plus an exact per-class
    score histogram, then recomputes the same three things from the store and
    compares. The histogram makes the divergence deciles exact rather than
    approximate, because scores are small integers.
    """
    if not os.path.exists(CONSOLIDATED_CSV):
        raise FileNotFoundError(f"{CONSOLIDATED_CSV} is already gone")

    # The CSV is one row per intersecting element, with a header; the store is one
    # row per unique tuple with a weight in $6. The two awk programs below
    # therefore differ only in `1` vs `$6` as the increment, and in the columns
    # they read (the CSV carries permutation_id in $1 and score in $2; the store
    # carries seed in $1 and score in $2, so the indices happen to line up).
    csv_prog = r"""
        NR == 1 { next }
        {
            pc[$1 "\t" $5] += 1
            pf[$1 "\t" $4] += 1
            hist[$5 "\t" $2] += 1
            n += 1
        }
        END {
            for (k in pc) print "PC\t" k "\t" pc[k]
            for (k in pf) print "PF\t" k "\t" pf[k]
            for (k in hist) print "H\t" k "\t" hist[k]
            print "N\t" n
        }
    """
    store_prog = r"""
        {
            pc[$1 "\t" $5] += $6
            pf[$1 "\t" $4] += $6
            hist[$5 "\t" $2] += $6
            n += $6
        }
        END {
            for (k in pc) print "PC\t" k "\t" pc[k]
            for (k in pf) print "PF\t" k "\t" pf[k]
            for (k in hist) print "H\t" k "\t" hist[k]
            print "N\t" n
        }
    """

    print(f"Streaming {CONSOLIDATED_CSV} ({os.path.getsize(CONSOLIDATED_CSV):,} B) ...")
    csv_out = sh(f"awk -F'\\t' {shlex.quote(csv_prog)} {shlex.quote(CONSOLIDATED_CSV)}")

    print("Recomputing the same statistics from the compacted store ...")
    store_out = sh(
        f"zstd -dc {shlex.quote(STORE_DIR)}/counts_seed_*.tsv.zst "
        f"| awk -F'\\t' {shlex.quote(store_prog)}"
    )

    def parse(blob: str):
        pc, pf, hist, total = {}, {}, {}, 0
        for line in blob.splitlines():
            parts = line.split("\t")
            kind = parts[0]
            if kind == "N":
                total = int(parts[1])
            elif kind == "PC":
                pc[(parts[1], parts[2])] = int(parts[3])
            elif kind == "PF":
                pf[(parts[1], parts[2])] = int(parts[3])
            elif kind == "H":
                hist[(parts[1], int(parts[2]))] = int(parts[3])
        return pc, pf, hist, total

    pc_a, pf_a, hist_a, n_a = parse(csv_out)
    pc_b, pf_b, hist_b, n_b = parse(store_out)

    problems = []
    if n_a != n_b:
        problems.append(f"total rows {n_a:,} (csv) vs {n_b:,} (store)")
    if pc_a != pc_b:
        diff = {k for k in set(pc_a) | set(pc_b) if pc_a.get(k) != pc_b.get(k)}
        problems.append(f"{len(diff)} (permutation, class) totals differ, e.g. {sorted(diff)[:3]}")
    if pf_a != pf_b:
        diff = {k for k in set(pf_a) | set(pf_b) if pf_a.get(k) != pf_b.get(k)}
        problems.append(f"{len(diff)} (permutation, family) totals differ, e.g. {sorted(diff)[:3]}")
    if hist_a != hist_b:
        diff = {k for k in set(hist_a) | set(hist_b) if hist_a.get(k) != hist_b.get(k)}
        problems.append(f"{len(diff)} (class, score) histogram bins differ, e.g. {sorted(diff)[:3]}")

    if problems:
        raise RuntimeError("aggregate check FAILED:\n  " + "\n  ".join(problems))

    summary = divergence_summary(hist_a)
    print("\nPer-class divergence distribution (identical in both sources):")
    for cls, stats in sorted(summary.items()):
        print(
            f"  {cls:12s} n={stats['n']:>12,}  mean={stats['mean']:7.2f}  "
            f"sd={stats['sd']:7.2f}  deciles={[int(d) for d in stats['deciles']]}"
        )

    return {
        "checked_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "consolidated_csv": os.path.basename(CONSOLIDATED_CSV),
        "consolidated_csv_bytes": os.path.getsize(CONSOLIDATED_CSV),
        "total_rows": n_a,
        "n_permutation_class_cells": len(pc_a),
        "n_permutation_family_cells": len(pf_a),
        "n_class_score_bins": len(hist_a),
        "all_checks_passed": True,
        "per_class_divergence": summary,
    }


def divergence_summary(hist: dict) -> dict:
    """Exact mean, SD and deciles per class from a (class, score) histogram."""
    by_class: dict[str, dict[int, int]] = {}
    for (cls, score), count in hist.items():
        by_class.setdefault(cls, {})[score] = by_class.setdefault(cls, {}).get(
            score, 0
        ) + count

    out = {}
    for cls, counts in by_class.items():
        total = sum(counts.values())
        mean = sum(s * c for s, c in counts.items()) / total
        var = sum(c * (s - mean) ** 2 for s, c in counts.items()) / total
        deciles = []
        target_idx = [total * q / 10 for q in range(1, 10)]
        cum, ti = 0, 0
        for score in sorted(counts):
            cum += counts[score]
            while ti < 9 and cum >= target_idx[ti]:
                deciles.append(score)
                ti += 1
        while len(deciles) < 9:
            deciles.append(max(counts))
        out[cls] = {
            "n": total,
            "mean": round(mean, 4),
            "sd": round(var**0.5, 4),
            "deciles": deciles,
        }
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--seeds", help="comma-separated seeds instead of 1..500")
    parser.add_argument("--limit", type=int, help="process only the first N seeds")
    parser.add_argument(
        "--delete-verified",
        action="store_true",
        help="delete each source BED once its own verification has passed",
    )
    parser.add_argument(
        "--verify-only",
        action="store_true",
        help="re-verify existing counts files against whatever sources remain",
    )
    parser.add_argument(
        "--aggregate-check",
        action="store_true",
        help="compare the whole store against consolidated_random_data.csv",
    )
    parser.add_argument(
        "--delete-consolidated",
        action="store_true",
        help="delete the 6.37 GB CSV (requires a passed aggregate check first)",
    )
    args = parser.parse_args()

    os.makedirs(STORE_DIR, exist_ok=True)
    manifest = load_manifest()

    if args.seeds:
        seeds = [int(s) for s in args.seeds.split(",")]
    else:
        seeds = list(range(1, N_SEEDS + 1))
    if args.limit:
        seeds = seeds[: args.limit]

    if args.delete_consolidated:
        check = manifest.get("aggregate_check")
        if not (check and check.get("all_checks_passed")):
            print(
                "REFUSING to delete: no passed aggregate check in MANIFEST.json. "
                "Run --aggregate-check first.",
                file=sys.stderr,
            )
            return 2
        if len(manifest["seeds"]) != N_SEEDS:
            print(
                f"REFUSING to delete: only {len(manifest['seeds'])}/{N_SEEDS} seeds "
                f"are compacted and verified.",
                file=sys.stderr,
            )
            return 2
        size = os.path.getsize(CONSOLIDATED_CSV)
        os.remove(CONSOLIDATED_CSV)
        manifest["consolidated_csv_deleted"] = True
        manifest["consolidated_csv_deleted_utc"] = dt.datetime.now(
            dt.timezone.utc
        ).isoformat(timespec="seconds")
        save_manifest(manifest)
        print(f"Deleted {CONSOLIDATED_CSV} ({size:,} B)")
        print(f"Free disk now: {free_bytes(REPO_DIR) / 1e9:.1f} GB")
        return 0

    if args.aggregate_check:
        manifest["aggregate_check"] = aggregate_check(manifest)
        save_manifest(manifest)
        print("\nAggregate check PASSED — the CSV is now provably redundant.")
        print("Delete it with --delete-consolidated.")
        return 0

    if args.verify_only:
        failures = []
        for seed in seeds:
            if not os.path.exists(counts_zst(seed)):
                failures.append(f"seed {seed}: no counts file")
                continue
            if not os.path.exists(source_bed(seed)):
                # Source already deleted: fall back to the manifest's record.
                record = manifest["seeds"].get(str(seed))
                if not record:
                    failures.append(f"seed {seed}: source gone and no manifest record")
                elif file_md5(counts_zst(seed)) != record["counts_md5"]:
                    failures.append(f"seed {seed}: counts md5 drifted from manifest")
                continue
            try:
                verify_seed(seed)
            except RuntimeError as exc:
                failures.append(str(exc))
        if failures:
            print("VERIFY FAILED:", file=sys.stderr)
            for line in failures:
                print("  " + line, file=sys.stderr)
            return 1
        print(f"Verified {len(seeds)} seed(s): all lossless.")
        return 0

    total_src = total_dst = 0
    for i, seed in enumerate(seeds, 1):
        src = source_bed(seed)
        if not os.path.exists(src):
            if str(seed) in manifest["seeds"]:
                print(f"[{i}/{len(seeds)}] seed {seed}: already compacted, source gone")
                continue
            raise FileNotFoundError(src)

        zst, n_tuples, size = encode_seed(seed)
        record = verify_seed(seed)
        manifest["seeds"][str(seed)] = record
        total_src += record["source_bytes"]
        total_dst += record["counts_bytes"]

        msg = (
            f"[{i}/{len(seeds)}] seed {seed}: {record['source_rows']:,} rows -> "
            f"{n_tuples:,} tuples, {record['source_bytes']:,} -> {size:,} B "
            f"({record['source_bytes'] / size:.1f}x) verified"
        )
        if args.delete_verified:
            os.remove(src)
            msg += ", source deleted"
        print(msg, flush=True)

        if i % 25 == 0 or i == len(seeds):
            save_manifest(manifest)

    save_manifest(manifest)
    if total_dst:
        print(
            f"\n{len(seeds)} seed(s): {total_src:,} B -> {total_dst:,} B "
            f"({total_src / total_dst:.1f}x)"
        )
    print(f"Store: {STORE_DIR}")
    print(f"Free disk: {free_bytes(REPO_DIR) / 1e9:.1f} GB")
    if not args.delete_verified:
        print("\nSources kept (no --delete-verified). Re-run with that flag to free disk.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
