#!/bin/bash
# WP5 / 05a — build the 5 kb and 20 kb TSS neighbourhoods and map TEs onto them.
#
# Reviewer major comment 5 asks for sensitivity analyses at alternative window
# sizes (5 kb and 20 kb). This script derives those windows from EXACTLY the same
# TSS definition that produced the published 10 kb set, so the three window sizes
# differ in width and in nothing else.
#
# How the TSS definition is recovered
# -----------------------------------
# `T2T_genes.bed` is the published 10 kb neighbourhood file (its interval set is
# provably identical to the `mapped_on_TSS` bedGraph the legacy permutations used
# as their -b file). 38,700 of its 38,704 intervals are exactly 10,000 bp wide;
# the four exceptions all start at 0, i.e. they were clipped where TSS - 5000
# fell off the start of the chromosome. So for every row, without exception:
#
#     TSS = end - 5000
#
# and the published window is [max(0, TSS-5000), TSS+5000]. This script applies
# the same construction with half-widths of 2,500 and 10,000, clamping to
# [0, chromosome length] at both ends — the 20 kb windows are wide enough that
# the upper clamp actually bites, which the 10 kb set never needed.
#
# Outputs (all in revision_G3/output/)
#   windows_5kb.bed  / windows_20kb.bed        chrom, start, end, gene
#   windows_10kb.bed                           a copy of T2T_genes.bed, for symmetry
#   observed_TEs_in_windows_{5,10,20}kb.tsv.zst    score, subfamily, family, class
#   observed_TEs_by_gene_{5,10,20}kb.tsv.zst       + the gene each hit belongs to
#
# The 4-column projection is the observed counterpart of the permutation
# background (same `cut -f4,5,6,7` projection), and feeds 05b_window_sensitivity.py.
# The gene-annotated version feeds 05c_percentile_sensitivity.py.

set -eo pipefail
export LC_ALL=C

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(dirname "${SCRIPT_DIR}")"
OUT_DIR="${SCRIPT_DIR}/output"

REPEATS="${REPO_DIR}/T2T_repeat_masker_processed_sorted.bed"
GENES_10KB="${REPO_DIR}/T2T_genes.bed"
GENOME_FILE="${REPO_DIR}/chm13.genome"

PUBLISHED_HALF_WIDTH=5000

mkdir -p "${OUT_DIR}"

for f in "${REPEATS}" "${GENES_10KB}" "${GENOME_FILE}"; do
  [ -s "$f" ] || { echo "ERROR: missing or empty input: $f" >&2; exit 1; }
done

# --------------------------------------------------------------------------- #
# 1. Build the window files
# --------------------------------------------------------------------------- #

build_windows () {
  local half_width=$1
  local label=$2
  local out="${OUT_DIR}/windows_${label}.bed"

  awk -v FS='\t' -v OFS='\t' \
      -v half="${half_width}" \
      -v published_half="${PUBLISHED_HALF_WIDTH}" \
      'NR == FNR { len[$1] = $2; next }
       {
         tss = $3 - published_half
         start = tss - half
         if (start < 0) start = 0
         end = tss + half
         if ($1 in len && end > len[$1]) end = len[$1]
         if (end > start) print $1, start, end, $4
         else dropped++
       }
       END { if (dropped) print "WARNING: dropped " dropped " empty windows" > "/dev/stderr" }' \
      "${GENOME_FILE}" "${GENES_10KB}" \
    | sort -k1,1 -k2,2n > "${out}"

  printf '  %-22s %8d intervals, widths: ' "$(basename "${out}")" "$(wc -l < "${out}")"
  awk -v FS='\t' '{print $3-$2}' "${out}" | sort -n | uniq -c \
    | sort -rn | head -3 | awk '{printf "%sx%s ", $1, $2}'
  printf '\n'
}

echo "Building TSS neighbourhoods from ${GENES_10KB} (TSS = end - ${PUBLISHED_HALF_WIDTH}):"
build_windows 2500  "5kb"
build_windows 10000 "20kb"

# Keep a 10 kb copy alongside the new ones so every downstream loop can treat the
# three window sizes uniformly. It must remain identical to the published file.
sort -k1,1 -k2,2n "${GENES_10KB}" > "${OUT_DIR}/windows_10kb.bed"
printf '  %-22s %8d intervals (sorted copy of the published set)\n' \
  "windows_10kb.bed" "$(wc -l < "${OUT_DIR}/windows_10kb.bed")"

# Sanity check: the 10 kb copy must have the same interval set as the source, and
# every 5 kb window must sit inside its 20 kb counterpart.
if ! diff -q <(cut -f1-3 "${GENES_10KB}" | sort) \
             <(cut -f1-3 "${OUT_DIR}/windows_10kb.bed" | sort) > /dev/null; then
  echo "ERROR: the 10 kb copy is not interval-identical to ${GENES_10KB}" >&2
  exit 1
fi
echo "  10 kb copy is interval-identical to the published set"

# --------------------------------------------------------------------------- #
# 2. Map TEs onto each window set
# --------------------------------------------------------------------------- #

map_tes () {
  local label=$1
  local windows="${OUT_DIR}/windows_${label}.bed"
  local proj="${OUT_DIR}/observed_TEs_in_windows_${label}.tsv.zst"
  local by_gene="${OUT_DIR}/observed_TEs_by_gene_${label}.tsv.zst"

  # -wa -wb keeps both sides, so one intersect serves both outputs. Columns:
  #   1-3 TE locus, 4 divergence score, 5 subfamily, 6 family, 7 class,
  #   8-10 window locus, 11 gene
  bedtools intersect -a "${REPEATS}" -b "${windows}" -wa -wb \
    | tee >(cut -f4,5,6,7 | zstd -q -19 -o "${proj}" -f) \
    | cut -f4,5,6,7,11 | zstd -q -19 -o "${by_gene}" -f

  local n_hits n_genes
  n_hits=$(zstd -dc "${proj}" | wc -l)
  n_genes=$(zstd -dc "${by_gene}" | cut -f5 | sort -u | wc -l)
  printf '  %-4s window: %10d TE-window intersections over %6d distinct genes  (%s, %s)\n' \
    "${label}" "${n_hits}" "${n_genes}" \
    "$(du -h "${proj}" | cut -f1)" "$(du -h "${by_gene}" | cut -f1)"
}

echo
echo "Mapping ${REPEATS##*/} onto each window set:"
map_tes "5kb"
map_tes "10kb"
map_tes "20kb"

echo
echo "Done. Outputs in ${OUT_DIR}/"
