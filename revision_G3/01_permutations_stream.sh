#!/bin/bash
# WP1 / WP5 — streaming permutation background, N = 500, one window size per run.
#
# Purpose changed under decision D1
# ---------------------------------
# This script does NOT re-run the published 10 kb background. That background
# stays exactly as published at N = 500 (D1), which is why no enrichment CSV,
# table or figure is regenerated for permutation reasons. What this script is for
# is the two NEW window sizes the reviewer asked for in major comment 5: 5 kb and
# 20 kb, also at N = 500 so they are directly comparable with the retained 10 kb
# background.
#
# It can also reproduce the 10 kb background on demand (`--window 10kb`), which
# is how `01a_consolidate_counts.py --check-legacy` proves this pipeline is
# equivalent to the legacy one before anything relies on it.
#
# Why streaming
# -------------
# The legacy pipeline wrote a full shuffled BED (162 MB) and a full intersect BED
# (11 MB) per seed, and then a 6.37 GB concatenation of all 500. This version
# never lands an intermediate: it shuffles, intersects, projects, and run-length
# encodes in one pipe, writing ~215 KB per seed. Output schema is identical to
# what 01b_compact_permutation_results.py produces, so compaction and new runs
# converge on one format and one reader (revision_lib.load_permutation_counts).
#
#     seed <TAB> score <TAB> subfamily_name <TAB> family_name <TAB> class_name <TAB> n
#
# Semantics preserved from the legacy run
# ---------------------------------------
#   * same -i file: T2T_repeat_masker_processed_sorted.bed is md5-identical to
#     the epigenomic_files/repeats_all.bed the legacy run shuffled
#   * same -g file: chm13.genome, rebuilt by the same two columns of
#     sequence_report.tsv the frozen download notebook used
#   * same seeds 1..500, same `bedtools shuffle | bedtools intersect -wa` order
#   * same `cut -f4,5,6,7` projection
# `bedtools shuffle -seed N` is deterministic, so identical inputs plus identical
# seeds give byte-identical results.
#
# Note on awk: every invocation sets -F'\t' and OFS='\t' because a small number of
# subfamilies have an EMPTY family_name. Under awk's default field splitting those
# rows look like 3 fields and class_name silently shifts into family_name.
#
# Usage
#   bash revision_G3/01_permutations_stream.sh --window 5kb
#   bash revision_G3/01_permutations_stream.sh --window 20kb
#   bash revision_G3/01_permutations_stream.sh --window 10kb --seeds 1  # regression
#   bash revision_G3/01_permutations_stream.sh --window 5kb --jobs 8 --resume

set -eo pipefail
export LC_ALL=C

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(dirname "${SCRIPT_DIR}")"
OUT_ROOT="${SCRIPT_DIR}/output"

WINDOW=""
NUM_PERMUTATIONS=500
N_JOBS="$(nproc)"
SEEDS=""
RESUME=0
ZSTD_LEVEL=19

while [ $# -gt 0 ]; do
  case "$1" in
    --window)   WINDOW="$2"; shift 2 ;;
    --jobs)     N_JOBS="$2"; shift 2 ;;
    --seeds)    SEEDS="$2"; shift 2 ;;
    --n)        NUM_PERMUTATIONS="$2"; shift 2 ;;
    --resume)   RESUME=1; shift ;;
    -h|--help)  sed -n '2,50p' "$0"; exit 0 ;;
    *)          echo "ERROR: unknown argument $1" >&2; exit 1 ;;
  esac
done

case "${WINDOW}" in
  5kb|10kb|20kb) ;;
  *) echo "ERROR: --window must be one of 5kb, 10kb, 20kb (got '${WINDOW}')" >&2; exit 1 ;;
esac

INPUT_REPEATS="${REPO_DIR}/T2T_repeat_masker_processed_sorted.bed"
GENOME_FILE="${REPO_DIR}/chm13.genome"
OUTPUT_DIR="${OUT_ROOT}/permutation_counts_${WINDOW}"

# 10 kb uses the published interval file itself; the other two use 05a's output.
if [ "${WINDOW}" = "10kb" ]; then
  TSS_FILE="${REPO_DIR}/T2T_genes.bed"
else
  TSS_FILE="${OUT_ROOT}/windows_${WINDOW}.bed"
fi

for f in "${INPUT_REPEATS}" "${GENOME_FILE}" "${TSS_FILE}"; do
  [ -s "$f" ] || { echo "ERROR: missing or empty input: $f" >&2; exit 1; }
done

# Writing the 10 kb store would overwrite the verified compaction of the
# published background, so that target is protected: regression runs go to a
# scratch directory instead.
if [ "${WINDOW}" = "10kb" ]; then
  OUTPUT_DIR="${OUT_ROOT}/permutation_counts_10kb_regression"
  echo "NOTE: --window 10kb writes to ${OUTPUT_DIR##*/} so the verified"
  echo "      compacted store of the published background is not touched."
fi

mkdir -p "${OUTPUT_DIR}"

export INPUT_REPEATS GENOME_FILE TSS_FILE OUTPUT_DIR ZSTD_LEVEL RESUME LC_ALL

run_one () {
  SEED=$1
  OUT="${OUTPUT_DIR}/counts_seed_${SEED}.tsv.zst"

  if [ "${RESUME}" = "1" ] && [ -s "${OUT}" ]; then
    return 0
  fi

  TMP="${OUTPUT_DIR}/.counts_seed_${SEED}.tsv.part"

  bedtools shuffle -i "${INPUT_REPEATS}" -g "${GENOME_FILE}" -seed "${SEED}" \
    | bedtools intersect -a - -b "${TSS_FILE}" -wa \
    | cut -f4,5,6,7 \
    | awk -F'\t' -v OFS='\t' -v seed="${SEED}" \
          '{ c[$0]++ } END { for (k in c) print seed, k, c[k] }' \
    | sort > "${TMP}"

  # An empty result means the pipeline failed rather than that nothing
  # intersected — with 3.7 M elements and 38,704 windows, zero is impossible.
  if [ ! -s "${TMP}" ]; then
    echo "ERROR: seed ${SEED} produced no output" >&2
    rm -f "${TMP}"
    return 1
  fi

  zstd -q -"${ZSTD_LEVEL}" -f "${TMP}" -o "${OUT}"
  rm -f "${TMP}"
}
export -f run_one

if [ -n "${SEEDS}" ]; then
  SEED_LIST="$(echo "${SEEDS}" | tr ',' '\n')"
else
  SEED_LIST="$(seq 1 "${NUM_PERMUTATIONS}")"
fi
N_TOTAL="$(echo "${SEED_LIST}" | wc -l)"

echo "Window        : ${WINDOW}  (${TSS_FILE##*/}, $(wc -l < "${TSS_FILE}") intervals)"
echo "Repeats       : ${INPUT_REPEATS##*/}  ($(wc -l < "${INPUT_REPEATS}") intervals)"
echo "Permutations  : ${N_TOTAL}"
echo "Parallelism   : ${N_JOBS} jobs"
echo "Output        : ${OUTPUT_DIR}"
echo

START=$(date +%s)
echo "${SEED_LIST}" | xargs -P "${N_JOBS}" -I{} bash -c 'run_one {}'
END=$(date +%s)

N_DONE="$(find "${OUTPUT_DIR}" -name 'counts_seed_*.tsv.zst' -size +0 | wc -l)"
echo
echo "Finished ${N_DONE} seed files in $(( END - START )) s"
du -sh "${OUTPUT_DIR}"

if [ "${N_DONE}" -lt "${N_TOTAL}" ]; then
  echo "WARNING: ${N_DONE} of ${N_TOTAL} expected files present; re-run with --resume" >&2
  exit 1
fi
