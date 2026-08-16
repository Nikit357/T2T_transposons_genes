#!/usr/bin/env bash
# Verify the PUBLISHED track hub over HTTP (WP13b).
#
# Reachability is not enough for a UCSC hub: the host must also honour byte-range requests and
# must not transcode the payload. So for every bigBed this script asks for bytes 0-3 and asserts
# they are the bigBed magic number eb f2 89 87. That can only hold if the byte offsets the
# browser asks for are the byte offsets it receives, which is the property UCSC actually needs.
#
# USAGE
#   bash revision_G3/12c_verify_trackhub_live.sh                       # local manifest
#   bash revision_G3/12c_verify_trackhub_live.sh --remote-manifest     # CI: fetch the manifest
#   bash revision_G3/12c_verify_trackhub_live.sh --base-url https://…/trackhub
#   bash revision_G3/12c_verify_trackhub_live.sh --json out.json
#
# Exit status is non-zero if any check fails, so it works as a CI gate and as a cron canary.

set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_URL="https://nikit357.github.io/T2T_transposons_genes/trackhub"
LOCAL_HUB="$HERE/trackhub"
MANIFEST_MODE="local"
JSON_OUT="$HERE/output/trackhub_live_check.json"
CURL=(curl --silent --show-error --location --max-time 60 --retry 2 --retry-delay 5)

while [ $# -gt 0 ]; do
  case "$1" in
    --base-url)        BASE_URL="${2%/}"; shift ;;
    --remote-manifest) MANIFEST_MODE="remote" ;;
    --json)            JSON_OUT="$2"; shift ;;
    -h|--help)         sed -n '2,20p' "$0"; exit 0 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
  shift
done

pass=0; fail=0; results=()
record() {  # record <status> <name> <detail>
  results+=("{\"check\": \"$2\", \"status\": \"$1\", \"detail\": \"$3\"}")
  if [ "$1" = "pass" ]; then pass=$((pass + 1)); printf '  PASS  %-38s %s\n' "$2" "$3"
  else fail=$((fail + 1)); printf '  FAIL  %-38s %s\n' "$2" "$3"; fi
}

echo "=== base URL: $BASE_URL ==="

echo
echo "=== 1. manifest ==="
work="$(mktemp -d)"; trap 'rm -rf "$work"' EXIT
if [ "$MANIFEST_MODE" = "remote" ]; then
  if "${CURL[@]}" -o "$work/MANIFEST.json" "$BASE_URL/MANIFEST.json"; then
    record pass manifest "fetched from the live site"
  else
    record fail manifest "cannot fetch $BASE_URL/MANIFEST.json"
  fi
  manifest="$work/MANIFEST.json"
else
  manifest="$LOCAL_HUB/MANIFEST.json"
  [ -s "$manifest" ] && record pass manifest "local $manifest" \
                     || record fail manifest "missing $manifest — run 12b first"
fi
# path<TAB>bytes, from the manifest, without a JSON dependency
python3 - "$manifest" > "$work/files.tsv" <<'PY' || true
import json, sys
for entry in json.load(open(sys.argv[1]))["files"]:
    print(entry["path"], entry["bytes"], sep="\t")
PY

echo
echo "=== 2. text files return 200 and match ==="
for rel in hub.txt genomes.txt hs1/trackDb.txt; do
  code=$("${CURL[@]}" -o "$work/$(basename "$rel")" -w '%{http_code}' "$BASE_URL/$rel")
  if [ "$code" = "200" ]; then
    if [ "$MANIFEST_MODE" = "local" ] && [ -f "$LOCAL_HUB/$rel" ] \
       && ! diff -q "$LOCAL_HUB/$rel" "$work/$(basename "$rel")" >/dev/null; then
      record fail "$rel" "200 but content differs from the local build"
    else
      record pass "$rel" "200, $(wc -c < "$work/$(basename "$rel")") bytes"
    fi
  else
    record fail "$rel" "HTTP $code"
  fi
done

echo
echo "=== 3. every bigBed: size, 206, and magic number at offset 0 ==="
while IFS=$'\t' read -r rel bytes; do
  case "$rel" in *.bb) ;; *) continue ;; esac
  headers="$("${CURL[@]}" -o /dev/null -D - -r 0-3 "$BASE_URL/$rel")"
  code=$(printf '%s' "$headers" | awk '/^HTTP/{c=$2} END{print c}')
  crange=$(printf '%s' "$headers" | tr -d '\r' | awk -F': ' 'tolower($1)=="content-range"{print $2}')
  magic=$("${CURL[@]}" -r 0-3 "$BASE_URL/$rel" | od -An -tx1 | tr -d ' \n')
  total="${crange##*/}"
  if   [ "$code" != "206" ];       then record fail "$rel range"  "expected 206, got $code"
  elif [ "$magic" != "ebf28987" ]; then record fail "$rel magic"  "bytes 0-3 = $magic, expected ebf28987 (payload transcoded?)"
  elif [ "$total" != "$bytes" ];   then record fail "$rel size"   "server says $total bytes, manifest says $bytes"
  else record pass "$rel" "206, $bytes bytes, magic ok"
  fi
done < "$work/files.tsv"

echo
echo "=== 4. description pages and landing page ==="
while IFS=$'\t' read -r rel bytes; do
  case "$rel" in *.html) ;; *) continue ;; esac
  code=$("${CURL[@]}" -o /dev/null -w '%{http_code}' "$BASE_URL/$rel")
  [ "$code" = "200" ] && record pass "$rel" "200" || record fail "$rel" "HTTP $code"
done < "$work/files.tsv"
root="${BASE_URL%/trackhub}"
code=$("${CURL[@]}" -o /dev/null -w '%{http_code}' "$root/")
[ "$code" = "200" ] && record pass "site root" "200 $root/" || record fail "site root" "HTTP $code"

echo
echo "=== 5. hubCheck against the live hub (optional) ==="
if command -v hubCheck >/dev/null 2>&1; then
  if hubCheck -checkSettings "$BASE_URL/hub.txt" > "$work/hubcheck.txt" 2>&1; then
    record pass hubCheck "clean"
  else
    # The missing local hg.conf is an environment limitation, not a hub defect (see 12_build).
    leftover=$(grep -vE "^Found [0-9]+ problem|Can't get default spec from host" "$work/hubcheck.txt" \
               | grep -cE "[^[:space:]]")
    [ "$leftover" -eq 0 ] && record pass hubCheck "clean apart from the missing local hg.conf" \
                          || record fail hubCheck "$(head -3 "$work/hubcheck.txt" | tr '\n' ' ')"
  fi
else
  echo "  hubCheck not installed; skipped (install from https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/)"
fi

mkdir -p "$(dirname "$JSON_OUT")"
{
  printf '{\n  "base_url": "%s",\n  "checked_utc": "%s",\n' "$BASE_URL" "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf '  "passed": %d,\n  "failed": %d,\n  "checks": [\n    ' "$pass" "$fail"
  ( IFS=$',\n    '; printf '%s' "${results[*]}" )
  printf '\n  ]\n}\n'
} > "$JSON_OUT"

echo
echo "$pass passed, $fail failed -> $JSON_OUT"
if [ "$fail" -gt 0 ]; then
  echo
  echo "If every check failed with 404, GitHub Pages is not enabled yet, or the gh-pages branch"
  echo "was not pushed. If only the range/magic checks failed, GitHub Pages is transcoding the"
  echo "payload and the raw.githubusercontent fallback in trackhub_ghpages_plan_260814.md §7 applies."
  exit 1
fi
echo
echo "one-click: https://genome.ucsc.edu/cgi-bin/hgTracks?db=hs1&hubUrl=$BASE_URL/hub.txt&position=chr9:21150692-21370055"
