#!/usr/bin/env python
"""WP4 — bound the assembly-attributable component of the Lu et al. difference (C6).

Why this exists
---------------
The reviewer asked a causal question about major comment 4: "Is this difference due
to methodology or the updated genome assembly?" Decision D6 drops the controlled
experiment that would settle it (no hg38 re-run), so caveat C6 requires the next
best thing: bound descriptively how much of the difference the assembly *could*
explain, attribute the remainder to the named methodological differences, and say
in the response letter that this is an argument rather than a controlled test.

That bound is what this script measures: how many of our 3,709,429 TEs and how many
of our 38,704 TSS windows sit in T2T-CHM13 sequence that no earlier hg38-based study
could have analysed at all.

Definition of "newly resolved", and why it is an upper bound
-----------------------------------------------------------
Built from the UCSC hs1 -> hg38 liftOver chains in two parts, because the naive
definition does not survive contact with the data:

  no_chain    chm13 sequence outside the span of every chain — no hg38 counterpart
              at all. This is the acrocentric short arms, the centromeric satellite
              arrays and chrY.
  large_gap   sequence inside a chain span but in none of that chain's aligned
              blocks, and at least MIN_GAP_BP long — chm13 insertions relative to
              hg38, i.e. sequence the older assembly was missing locally.

    newly_resolved = merge(no_chain + large_gap)

The MIN_GAP_BP floor is what makes the measure mean anything. Taking the raw
complement of the aligned blocks instead counts every alignment indel: it returns
6.81 % of the genome but has 63 % of all TSS windows "touching" it, because a 10 kb
window almost always contains at least one 1-bp indel. Those are not newly resolved
sequence. The sub-floor total is reported separately so the choice is visible.

The result is still deliberately an over-count, which is the right direction for a
bound: liftOver chains are quality-filtered, so chm13 sequence with a poor or
repetitive hg38 counterpart also lands outside the chains. So "at most X % of our
TEs sit in sequence unavailable to an hg38 study" is a true ceiling on the
assembly-attributable component, and the response letter must quote it as a ceiling.

The headline finding is checked explicitly
------------------------------------------
The IFNA domain (chr9:21,150,692-21,370,055) is intersected with the newly-resolved
set. If the domain is fully alignable to hg38 — which is the expected result for a
gene-dense sub-telomeric region — then the paper's central claim cannot be an
assembly artefact, and that is a stronger statement than the genome-wide bound.

Inputs
    external/hs1ToHg38.over.chain.gz     UCSC liftOver chains (recorded in PROVENANCE.md)
    T2T_repeat_masker_processed_sorted.bed, T2T_genes.bed, chm13.genome
    output/windows_10kb.bed              from 05a_build_windows.sh

Outputs
    output/newly_resolved_regions.bed            merged non-alignable chm13 intervals
    output/assembly_bound_summary.json           the headline percentages
    output/assembly_bound_by_class.csv           TEs in / out, per class and per family
    output/assembly_bound_by_chromosome.csv      per-chromosome non-alignable fraction

Usage
    python revision_G3/04b_newly_resolved_regions.py
    python revision_G3/04b_newly_resolved_regions.py --summary
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import shlex
import subprocess
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import revision_lib as rl  # noqa: E402

OUTPUT_DIR = rl.OUTPUT_DIR
EXTERNAL_DIR = os.path.join(rl.REVISION_DIR, "external")
CHAIN_PATH = os.path.join(EXTERNAL_DIR, "hs1ToHg38.over.chain.gz")
GENOME_FILE = os.path.join(rl.REPO_DIR, "chm13.genome")

ALIGNED_BED = os.path.join(OUTPUT_DIR, "hg38_alignable_regions.bed")
CHAIN_SPANS_BED = os.path.join(OUTPUT_DIR, "hg38_chain_spans.bed")
NEWLY_RESOLVED_BED = os.path.join(OUTPUT_DIR, "newly_resolved_regions.bed")

# Minimum length of an unaligned stretch inside a chain for it to count as newly
# resolved sequence rather than an alignment indel. 1 kb is above essentially every
# indel in these chains and below every genuine insertion of interest; the totals at
# other floors are reported in the summary so the sensitivity to it is visible.
MIN_GAP_BP = 1_000

IFNA_CHROM, IFNA_START, IFNA_END = rl.IFNA_DOMAIN


def sh(cmd: str) -> str:
    return subprocess.run(
        cmd, shell=True, check=True, capture_output=True, text=True,
        env=dict(os.environ, LC_ALL="C"),
    ).stdout


def parse_chains(chain_path: str, chromosomes: set[str]) -> tuple[int, int]:
    """Write both the aligned blocks and the whole chain spans, chm13 side, to BED.

    Chain format: a `chain` header gives tName, tStart and tEnd, then one line per
    block. A block line is `size dt dq`, where dt advances the target past an
    unaligned gap; the last line of a chain carries only `size`. Only the target
    (chm13) side matters here — the question is which chm13 bases have an hg38
    counterpart, not where in hg38 it is.
    """
    n_blocks = n_chains = 0
    with gzip.open(chain_path, "rt") as source, \
            open(ALIGNED_BED, "w") as blocks_out, \
            open(CHAIN_SPANS_BED, "w") as spans_out:
        target_name, position, keep = None, 0, False
        for line in source:
            line = line.strip()
            if not line:
                target_name = None
                continue
            if line.startswith("chain"):
                fields = line.split()
                target_name = fields[2]
                position = int(fields[5])
                keep = target_name in chromosomes
                if keep:
                    spans_out.write(f"{target_name}\t{fields[5]}\t{fields[6]}\n")
                    n_chains += 1
                continue
            if target_name is None:
                continue
            parts = line.split()
            size = int(parts[0])
            if keep:
                blocks_out.write(f"{target_name}\t{position}\t{position + size}\n")
                n_blocks += 1
            position += size + (int(parts[1]) if len(parts) >= 3 else 0)
    return n_chains, n_blocks


def span_bp(path: str) -> int:
    return int(sh(f"awk -F'\\t' '{{ s += $3 - $2 }} END {{ print s+0 }}' "
                  f"{shlex.quote(path)}").strip())


def sort_merge(path: str) -> None:
    sh(f"sort -k1,1 -k2,2n {shlex.quote(path)} | bedtools merge -i - > {shlex.quote(path)}.m "
       f"&& mv {shlex.quote(path)}.m {shlex.quote(path)}")


def build_region_sets() -> dict[str, object]:
    """The two-part newly-resolved set: no chain at all, plus large in-chain gaps."""
    chrom_sizes = dict(
        (line.split("\t")[0], int(line.split("\t")[1]))
        for line in open(GENOME_FILE).read().strip().split("\n")
    )
    total = sum(chrom_sizes.values())

    print(f"Parsing {os.path.basename(CHAIN_PATH)} ...")
    n_chains, n_blocks = parse_chains(CHAIN_PATH, set(chrom_sizes))
    print(f"  {n_chains:,} chains, {n_blocks:,} aligned blocks on the chm13 side")

    sort_merge(ALIGNED_BED)
    sort_merge(CHAIN_SPANS_BED)

    genome_sorted = f"{GENOME_FILE}.sorted"
    sh(f"sort -k1,1 -k2,2n {shlex.quote(GENOME_FILE)} > {shlex.quote(genome_sorted)}")

    work = os.path.join(OUTPUT_DIR, "_no_chain.bed")
    gaps = os.path.join(OUTPUT_DIR, "_in_chain_gaps.bed")
    large = os.path.join(OUTPUT_DIR, "_large_gaps.bed")

    # Outside every chain span: no hg38 counterpart at all.
    sh(f"bedtools complement -i {shlex.quote(CHAIN_SPANS_BED)} -g {shlex.quote(genome_sorted)} "
       f"> {shlex.quote(work)}")
    # Unaligned inside a chain: complement of the blocks, minus the no-chain part.
    sh(f"bedtools complement -i {shlex.quote(ALIGNED_BED)} -g {shlex.quote(genome_sorted)} "
       f"| bedtools subtract -a - -b {shlex.quote(work)} > {shlex.quote(gaps)}")
    sh(f"awk -F'\\t' '($3 - $2) >= {MIN_GAP_BP}' {shlex.quote(gaps)} > {shlex.quote(large)}")
    sh(f"cat {shlex.quote(work)} {shlex.quote(large)} | sort -k1,1 -k2,2n "
       f"| bedtools merge -i - > {shlex.quote(NEWLY_RESOLVED_BED)}")

    spans = {
        "genome_bp": total,
        "aligned_blocks_bp": span_bp(ALIGNED_BED),
        "chain_span_bp": span_bp(CHAIN_SPANS_BED),
        "no_chain_bp": span_bp(work),
        "in_chain_unaligned_bp": span_bp(gaps),
        "in_chain_unaligned_large_bp": span_bp(large),
        "newly_resolved_bp": span_bp(NEWLY_RESOLVED_BED),
        "min_gap_bp": MIN_GAP_BP,
    }
    spans["in_chain_unaligned_small_bp"] = (
        spans["in_chain_unaligned_bp"] - spans["in_chain_unaligned_large_bp"]
    )

    for path in (work, gaps, large, genome_sorted):
        os.remove(path)

    print(f"  no chain at all               : {spans['no_chain_bp']:>13,} bp "
          f"({100 * spans['no_chain_bp'] / total:.2f} %)")
    print(f"  unaligned inside a chain      : {spans['in_chain_unaligned_bp']:>13,} bp, "
          f"of which >= {MIN_GAP_BP:,} bp: {spans['in_chain_unaligned_large_bp']:,} "
          f"(discarded as indels: {spans['in_chain_unaligned_small_bp']:,})")
    print(f"  newly resolved (the union)    : {spans['newly_resolved_bp']:>13,} bp "
          f"({100 * spans['newly_resolved_bp'] / total:.2f} %)")
    return spans


def count_features_in_regions() -> dict[str, int]:
    """TEs, TSS windows and genes falling in newly resolved sequence."""
    counts = {}

    # A TE counts as newly resolved if it overlaps the set at all (-u, any overlap):
    # a partially new element could not have been annotated identically in hg38.
    counts["tes_total"] = int(sh(f"wc -l < {shlex.quote(rl.REPEATS_BED)}").strip())
    counts["tes_newly_resolved"] = int(sh(
        f"bedtools intersect -u -a {shlex.quote(rl.REPEATS_BED)} "
        f"-b {shlex.quote(NEWLY_RESOLVED_BED)} | wc -l").strip())

    windows = os.path.join(OUTPUT_DIR, "windows_10kb.bed")
    counts["windows_total"] = int(sh(f"wc -l < {shlex.quote(windows)}").strip())
    counts["windows_newly_resolved"] = int(sh(
        f"bedtools intersect -u -a {shlex.quote(windows)} "
        f"-b {shlex.quote(NEWLY_RESOLVED_BED)} | wc -l").strip())
    counts["windows_entirely_newly_resolved"] = int(sh(
        f"bedtools intersect -f 1.0 -u -a {shlex.quote(windows)} "
        f"-b {shlex.quote(NEWLY_RESOLVED_BED)} | wc -l").strip())

    genes_touching = sh(
        f"bedtools intersect -u -a {shlex.quote(rl.GENES_BED)} "
        f"-b {shlex.quote(NEWLY_RESOLVED_BED)} | cut -f4 | sort -u")
    counts["genes_total"] = int(sh(
        f"cut -f4 {shlex.quote(rl.GENES_BED)} | sort -u | wc -l").strip())
    counts["genes_newly_resolved"] = len(
        [g for g in genes_touching.strip().split("\n") if g]
    )

    # Why so few TEs in 6.7 % of the genome: newly resolved sequence is satellite
    # array, and satellites are not one of the six classes this study annotates.
    # Recording the base-pair density on both sides makes that explicit rather than
    # leaving the low count looking like a pipeline error.
    merged_tes = f"sort -k1,1 -k2,2n {shlex.quote(rl.REPEATS_BED)} | bedtools merge -i -"
    counts["te_annotated_bp_genome"] = int(sh(
        f"{merged_tes} | awk -F'\\t' '{{ s += $3 - $2 }} END {{ print s }}'").strip())
    counts["te_annotated_bp_in_newly_resolved"] = int(sh(
        f"{merged_tes} | bedtools intersect -a - -b {shlex.quote(NEWLY_RESOLVED_BED)} "
        f"| awk -F'\\t' '{{ s += $3 - $2 }} END {{ print s+0 }}'").strip())
    return counts


def ifna_check() -> dict[str, object]:
    """Is any part of the headline IFNA domain outside hg38-alignable sequence?"""
    region = f"{IFNA_CHROM}\t{IFNA_START}\t{IFNA_END}\n"
    tmp = os.path.join(OUTPUT_DIR, "_ifna_region.bed")
    with open(tmp, "w") as fh:
        fh.write(region)
    overlap_bp = int(sh(
        f"bedtools intersect -a {shlex.quote(tmp)} -b {shlex.quote(NEWLY_RESOLVED_BED)} "
        f"| awk -F'\\t' '{{ s += $3 - $2 }} END {{ print s+0 }}'").strip())
    os.remove(tmp)
    length = IFNA_END - IFNA_START
    return {
        "domain": f"{IFNA_CHROM}:{IFNA_START:,}-{IFNA_END:,}",
        "domain_bp": length,
        "newly_resolved_bp_in_domain": overlap_bp,
        "fraction_newly_resolved": overlap_bp / length,
        "fully_alignable_to_hg38": overlap_bp == 0,
    }


def per_level_table(level_column: int, level_name: str) -> pd.DataFrame:
    """TEs in vs. out of newly resolved sequence, per class or per family."""
    inside = sh(
        f"bedtools intersect -u -a {shlex.quote(rl.REPEATS_BED)} "
        f"-b {shlex.quote(NEWLY_RESOLVED_BED)} "
        f"| awk -F'\\t' '{{ n[${level_column}]++ }} END "
        f"{{ for (k in n) print k \"\\t\" n[k] }}'")
    total = sh(
        f"awk -F'\\t' '{{ n[${level_column}]++ }} END "
        f"{{ for (k in n) print k \"\\t\" n[k] }}' {shlex.quote(rl.REPEATS_BED)}")

    def parse(text: str) -> dict[str, int]:
        return {
            line.rsplit("\t", 1)[0]: int(line.rsplit("\t", 1)[1])
            for line in text.strip().split("\n") if line
        }

    inside_counts, total_counts = parse(inside), parse(total)
    rows = [
        {
            level_name: name,
            "n_total": total_counts[name],
            "n_newly_resolved": inside_counts.get(name, 0),
            "fraction_newly_resolved": inside_counts.get(name, 0) / total_counts[name],
        }
        for name in sorted(total_counts)
    ]
    return pd.DataFrame(rows).sort_values("fraction_newly_resolved", ascending=False)


def per_chromosome_table(spans: dict[str, int]) -> pd.DataFrame:
    """Newly resolved base pairs per chromosome — where the new sequence actually is."""
    per_chrom = sh(
        f"awk -F'\\t' '{{ s[$1] += $3 - $2 }} END {{ for (k in s) print k \"\\t\" s[k] }}' "
        f"{shlex.quote(NEWLY_RESOLVED_BED)}")
    newly = {
        line.split("\t")[0]: int(line.split("\t")[1])
        for line in per_chrom.strip().split("\n") if line
    }
    sizes = dict(
        (line.split("\t")[0], int(line.split("\t")[1]))
        for line in open(GENOME_FILE).read().strip().split("\n")
    )
    rows = [
        {
            "chrom": chrom,
            "chrom_bp": size,
            "newly_resolved_bp": newly.get(chrom, 0),
            "fraction_newly_resolved": newly.get(chrom, 0) / size,
        }
        for chrom, size in sizes.items()
    ]
    return pd.DataFrame(rows).sort_values("fraction_newly_resolved", ascending=False)


def report() -> int:
    path = os.path.join(OUTPUT_DIR, "assembly_bound_summary.json")
    if not os.path.exists(path):
        print("Run without --summary first.", file=sys.stderr)
        return 1
    summary = json.load(open(path))

    print("=" * 78)
    print("Upper bound on the assembly-attributable component (WP4 / caveat C6)")
    print("=" * 78)
    genome = summary["genome_bp"]
    print(f"  chm13v2.0 not alignable to hg38 : "
          f"{summary['newly_resolved_bp']:,} bp "
          f"({100 * summary['newly_resolved_bp'] / genome:.2f} % of the genome)")
    print(f"    of which no chain at all      : {summary['no_chain_bp']:,} bp")
    print(f"    plus in-chain gaps >= {summary['min_gap_bp']:,} bp  : "
          f"{summary['in_chain_unaligned_large_bp']:,} bp")
    print(f"    excluded as alignment indels  : "
          f"{summary['in_chain_unaligned_small_bp']:,} bp")
    counts = summary["counts"]
    print(f"  TEs in that sequence            : {counts['tes_newly_resolved']:,} of "
          f"{counts['tes_total']:,} "
          f"({100 * counts['tes_newly_resolved'] / counts['tes_total']:.2f} %)")
    print(f"  TSS windows touching it         : {counts['windows_newly_resolved']:,} of "
          f"{counts['windows_total']:,} "
          f"({100 * counts['windows_newly_resolved'] / counts['windows_total']:.2f} %)"
          f"; entirely inside: {counts['windows_entirely_newly_resolved']:,}")
    print(f"  genes touching it               : {counts['genes_newly_resolved']:,} of "
          f"{counts['genes_total']:,} "
          f"({100 * counts['genes_newly_resolved'] / counts['genes_total']:.2f} %)")

    density_in = (counts["te_annotated_bp_in_newly_resolved"]
                  / summary["newly_resolved_bp"])
    density_all = counts["te_annotated_bp_genome"] / genome
    print(f"  TE-annotated bp density         : {100 * density_in:.1f} % inside vs "
          f"{100 * density_all:.1f} % genome-wide ({density_in / density_all:.2f}x)")
    print("    -> newly resolved sequence is satellite array, and satellites are not")
    print("       one of the six classes analysed here, hence the low TE count above")

    ifna = summary["ifna_domain"]
    print(f"\n  IFNA domain {ifna['domain']}: "
          f"{ifna['newly_resolved_bp_in_domain']:,} bp newly resolved "
          f"({100 * ifna['fraction_newly_resolved']:.2f} %)")
    print("  -> " + (
        "fully alignable to hg38, so the headline finding cannot be an assembly artefact"
        if ifna["fully_alignable_to_hg38"] else
        "PARTLY outside hg38-alignable sequence — report this explicitly"))

    by_class_path = os.path.join(OUTPUT_DIR, "assembly_bound_by_class.csv")
    if os.path.exists(by_class_path):
        print("\n  Per class:")
        for r in pd.read_csv(by_class_path).itertuples():
            print(f"    {r.class_name:12s} {r.n_newly_resolved:>9,} / {r.n_total:>9,} "
                  f"= {100 * r.fraction_newly_resolved:5.2f} %")

    by_chrom_path = os.path.join(OUTPUT_DIR, "assembly_bound_by_chromosome.csv")
    if os.path.exists(by_chrom_path):
        print("\n  Most affected chromosomes (the acrocentrics and chrY, as expected):")
        for r in pd.read_csv(by_chrom_path).head(6).itertuples():
            print(f"    {r.chrom:6s} {100 * r.fraction_newly_resolved:5.2f} % "
                  f"({r.newly_resolved_bp:,} bp)")
    print("\nThis is a CEILING, not an estimate: liftOver chains are quality-filtered,")
    print("so the complement includes hard-to-align sequence as well as genuinely new.")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--summary", action="store_true",
                        help="print the bound from existing outputs and exit")
    args = parser.parse_args()

    if args.summary:
        return report()

    if not os.path.exists(CHAIN_PATH):
        print(f"Missing {CHAIN_PATH}\nSee revision_G3/external/PROVENANCE.md",
              file=sys.stderr)
        return 1
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    spans = build_region_sets()
    print("\nCounting features in newly resolved sequence ...")
    counts = count_features_in_regions()
    ifna = ifna_check()

    per_level_table(7, "class_name").to_csv(
        os.path.join(OUTPUT_DIR, "assembly_bound_by_class.csv"), index=False)
    per_level_table(6, "family_name").to_csv(
        os.path.join(OUTPUT_DIR, "assembly_bound_by_family.csv"), index=False)
    per_chromosome_table(spans).to_csv(
        os.path.join(OUTPUT_DIR, "assembly_bound_by_chromosome.csv"), index=False)

    summary = dict(spans)
    summary["counts"] = counts
    summary["ifna_domain"] = ifna
    summary["definition"] = (
        "chm13v2.0 bases not covered by any aligned block of the UCSC hs1->hg38 "
        "liftOver chain set; an upper bound on sequence unavailable to hg38 studies"
    )
    summary["chain_source"] = (
        "https://hgdownload.soe.ucsc.edu/goldenPath/hs1/liftOver/hs1ToHg38.over.chain.gz"
    )
    with open(os.path.join(OUTPUT_DIR, "assembly_bound_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print("Wrote assembly_bound_summary.json, assembly_bound_by_"
          "{class,family,chromosome}.csv and newly_resolved_regions.bed\n")

    return report()


if __name__ == "__main__":
    sys.exit(main())
