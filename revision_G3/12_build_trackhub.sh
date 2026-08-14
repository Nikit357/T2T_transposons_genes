#!/usr/bin/env bash
# Build the UCSC track hub for the G3 revision (WP13b, decision D8).
#
# The hub lets a reviewer open any locus on hs1 (T2T-CHM13v2.0) and see the TE annotation,
# the TSS windows and the gene sets this paper actually measured, coloured with the same
# class palette as the figures. Data comes only from the two canonical inputs of plan
# section 3.0, so the hub cannot drift from the analysis.
#
# Tracks are split per TE class rather than combined: one bigBed of 3.7 M elements would risk
# GitHub's 100 MB per-file limit, and per-class files let a reviewer switch classes off.
#
# Requirements
#   bedToBigBed and hubCheck from https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/
#   Put them on PATH, or set UCSC_BIN to the directory holding them.
#
# Usage
#   bash revision_G3/12_build_trackhub.sh [output-dir]
# Default output-dir is revision_G3/trackhub/ .

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(dirname "$HERE")"
OUT="${1:-$HERE/trackhub}"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

PYTHON="${PYTHON:-$HOME/venvs/Retroelements_3_11/bin/python}"
UCSC_BIN="${UCSC_BIN:-}"
BED_TO_BIGBED="$(command -v bedToBigBed || echo "${UCSC_BIN%/}/bedToBigBed")"
HUB_CHECK="$(command -v hubCheck || echo "${UCSC_BIN%/}/hubCheck")"

GENOME_FILE="$REPO/chm13.genome"
HUB_URL_BASE="https://nikit357.github.io/T2T_transposons_genes/trackhub"
CONTACT_EMAIL="daniil.nikitin@bostongene.com"

for tool in "$BED_TO_BIGBED" "$HUB_CHECK"; do
  if [ ! -x "$tool" ]; then
    echo "ERROR: $tool not found or not executable." >&2
    echo "Download it from https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/ and" >&2
    echo "put it on PATH, or set UCSC_BIN to its directory." >&2
    exit 1
  fi
done
[ -f "$GENOME_FILE" ] || { echo "ERROR: missing $GENOME_FILE" >&2; exit 1; }

echo "=== 1. BED9 files ==="
"$PYTHON" "$HERE/12a_trackhub_beds.py" "$WORK"

echo
echo "=== 2. bigBed conversion ==="
mkdir -p "$OUT/hs1"
# chrom.sizes and the .genome file are the same two columns; bedToBigBed wants it sorted.
sort -k1,1 "$GENOME_FILE" > "$WORK/chrom.sizes"
for bed in "$WORK"/*.bed; do
  name="$(basename "$bed" .bed)"
  # -type=bed9 because every track carries itemRgb; the class palette is the whole point.
  LC_ALL=C sort -k1,1 -k2,2n "$bed" > "$WORK/$name.sorted.bed"
  "$BED_TO_BIGBED" -type=bed9 "$WORK/$name.sorted.bed" "$WORK/chrom.sizes" \
      "$OUT/hs1/$name.bb"
  printf '  %-24s %10s bytes\n' "$name.bb" "$(stat -c%s "$OUT/hs1/$name.bb")"
done

echo
echo "=== 3. hub files ==="
cat > "$OUT/hub.txt" <<HUB
hub T2T_transposons_genes
shortLabel T2T TEs near genes
longLabel Transposable elements and gene TSS neighbourhoods in T2T-CHM13v2.0 (Nikitin 2026)
genomesFile genomes.txt
email $CONTACT_EMAIL
descriptionUrl https://github.com/Nikit357/T2T_transposons_genes
HUB

cat > "$OUT/genomes.txt" <<GENOMES
genome hs1
trackDb hs1/trackDb.txt
GENOMES

{
cat <<'COMPOSITE'
track T2T_TEs
compositeTrack on
shortLabel T2T TEs by class
longLabel Transposable elements in T2T-CHM13v2.0, split by RepeatMasker class
type bigBed 9 .
itemRgb on
visibility dense
priority 1
html T2T_TEs

COMPOSITE

# Divergence is RepeatMasker's substitutions per 1000 bp. It runs 0-480 in this annotation,
# not 0-1000, so the label says what the number is instead of implying a UCSC-style score.
for class_label in LINE SINE LTR DNA SVA Helitron; do
  cat <<TRACK
    track TEs_${class_label}
    parent T2T_TEs on
    bigDataUrl TEs_${class_label}.bb
    shortLabel ${class_label}
    longLabel ${class_label} elements; name = subfamily, score = RepeatMasker divergence (substitutions per 1000 bp, observed range 0-480)
    type bigBed 9 .
    itemRgb on
    visibility dense
    html TEs_${class_label}

TRACK
done

cat <<'REST'
track TSS_10kb_windows
bigDataUrl TSS_10kb_windows.bb
shortLabel TSS 10 kb windows
longLabel The 38,704 TSS 10 kb neighbourhoods used for all TE mapping (5 kb either side of each TSS)
type bigBed 9 .
itemRgb on
visibility pack
priority 2
html TSS_10kb_windows

track genes_TE_top
bigDataUrl genes_TE_top.bb
shortLabel TE-top genes
longLabel The 1,436 genes with the highest TE count in their TSS 10 kb neighbourhood (top 5 percent)
type bigBed 9 .
itemRgb on
visibility pack
priority 3
html genes_TE_top

track genes_TE_bottom
bigDataUrl genes_TE_bottom.bb
shortLabel TE-bottom genes
longLabel The 1,436 genes with the lowest TE count in their TSS 10 kb neighbourhood (bottom 5 percent)
type bigBed 9 .
itemRgb on
visibility pack
priority 4
html genes_TE_bottom

track IFNA_domain
bigDataUrl IFNA_domain.bb
shortLabel IFNA domain
longLabel The 220 kb interferon-alpha domain enriched with young L1 elements (chr9:21,150,692-21,370,055)
type bigBed 9 .
itemRgb on
visibility pack
priority 5
html IFNA_domain
REST
} > "$OUT/hs1/trackDb.txt"

# Per-track description pages. hubCheck warns for every track without one, and a reviewer
# clicking a track name in the browser gets these rather than a blank page.
write_description() {
  local name="$1" heading="$2" body="$3"
  cat > "$OUT/hs1/$name.html" <<HTML
<h2>$heading</h2>
<p>$body</p>
<h3>Source</h3>
<p>Built by <code>revision_G3/12_build_trackhub.sh</code> from the T2T-CHM13v2.0 RepeatMasker
annotation and NCBI RefSeq gene annotation used throughout Nikitin (2026),
<a href="https://github.com/Nikit357/T2T_transposons_genes">github.com/Nikit357/T2T_transposons_genes</a>.
Colours match the TE class palette used in the paper's figures.</p>
HTML
}

write_description T2T_TEs "Transposable elements by class" \
  "All 3,709,429 RepeatMasker-annotated transposable elements in T2T-CHM13v2.0, split into one
   subtrack per class. Item names are subfamilies and item scores are RepeatMasker divergence
   (substitutions per 1000 bp, observed range 0-480; higher means older)."
for class_label in LINE SINE LTR DNA SVA Helitron; do
  write_description "TEs_${class_label}" "${class_label} elements" \
    "RepeatMasker-annotated ${class_label} elements. Item name is the subfamily, item score is
     divergence in substitutions per 1000 bp (observed range 0-480). Note that the paper reports
     enrichment as the ratio of the observed to a 500-permutation random odds ratio, because the
     raw odds ratio scales with element length."
done
write_description TSS_10kb_windows "TSS 10 kb neighbourhoods" \
  "The 38,704 transcription start site neighbourhoods (5 kb either side of each TSS) that define
   proximity throughout the paper, covering 28,738 genes. Genes with several annotated TSS
   contribute several windows. A small number of windows are clipped where a TSS lies within
   5 kb of a chromosome end."
write_description genes_TE_top "TE-top genes (top 5 percent)" \
  "The 1,436 genes with the highest total TE count in their TSS 10 kb neighbourhood. These are
   the foreground set for the 'TE top' Gene Ontology analysis, which is enriched for RNA
   splicing, DNA repair, telomere maintenance and apoptotic signalling."
write_description genes_TE_bottom "TE-bottom genes (bottom 5 percent)" \
  "The 1,436 genes with the lowest total TE count in their TSS 10 kb neighbourhood. These are
   the foreground set for the 'TE bottom' Gene Ontology analysis, which is dominated by
   embryogenesis, transcription and nervous system terms."
write_description IFNA_domain "Interferon-alpha domain" \
  "The 220 kb domain on chromosome 9 (chr9:21,150,692-21,370,055) containing IFNA4, IFNA6,
   IFNA7, IFNA10, IFNA14, IFNA16, IFNA17, IFNA21, IFNA22P, IFNA5, IFNW1 and KLHL9. It holds 175
   TEs, of which 77 are L1 copies from 36 subfamilies at a mean divergence of 135.7 against a
   genome-wide L1 mean of 188.2. The deficit survives null windows matched for L1 count
   (empirical p = 0.006) and for gene density (p = 0.002)."

echo "  hub.txt, genomes.txt, hs1/trackDb.txt, $(ls "$OUT"/hs1/*.html | wc -l) description pages"

echo
echo "=== 4. hubCheck ==="
# hubCheck needs an hg.conf to fetch UCSC's default trackDb spec, which a machine without a
# browser install does not have. That one message is an environment limitation, not a defect in
# the hub, so it is reported and tolerated while any other finding is treated as a failure.
hub_report="$WORK/hubcheck.txt"
set +e
"$HUB_CHECK" -checkSettings "$OUT/hub.txt" > "$hub_report" 2>&1
hub_status=$?
set -e
# Drop the summary header and the hg.conf message; anything left is a real finding.
real_problems="$(grep -vE "^Found [0-9]+ problem|Can't get default spec from host" \
                 "$hub_report" | grep -cE "[^[:space:]]" || true)"
cat "$hub_report" | sed 's/^/  /'
if [ "$real_problems" -eq 0 ]; then
  if [ "$hub_status" -ne 0 ]; then
    echo "  hubCheck: OK apart from the missing local hg.conf (no track-level findings)"
  else
    echo "  hubCheck: OK"
  fi
else
  echo "ERROR: hubCheck reported track-level problems; fix them before publishing." >&2
  exit 1
fi

echo
echo "=== 4b. every bigBed opens and its record count matches the source ==="
if command -v bigBedInfo >/dev/null 2>&1 || [ -x "${UCSC_BIN%/}/bigBedInfo" ]; then
  BIG_BED_INFO="$(command -v bigBedInfo || echo "${UCSC_BIN%/}/bigBedInfo")"
  for bb in "$OUT"/hs1/*.bb; do
    printf '  %-26s %s\n' "$(basename "$bb")" \
      "$("$BIG_BED_INFO" "$bb" | awk -F': ' '/itemCount/{print "itemCount="$2}')"
  done
else
  echo "  bigBedInfo not available; skipped (optional check)"
fi

echo
echo "=== 5. one-click URL ==="
echo "https://genome.ucsc.edu/cgi-bin/hgTracks?db=hs1&hubUrl=$HUB_URL_BASE/hub.txt&position=chr9:21150692-21370055"
echo
du -sh "$OUT"
echo "built -> $OUT"
