#!/usr/bin/env bash
# Publish the built UCSC track hub to the gh-pages branch (WP13b, decision D8).
#
# WHY THIS EXISTS
#   12_build_trackhub.sh builds revision_G3/trackhub/ (105 MB), which is gitignored on the
#   analysis branch. The manuscript prints
#       https://nikit357.github.io/T2T_transposons_genes/trackhub/hub.txt
#   so the published tree must carry a TOP-LEVEL trackhub/ directory and that path may not
#   move. On 2026-08-14 the URL returned 404 for one reason only: the gh-pages branch had
#   never been pushed and GitHub Pages had never been enabled (repos API has_pages=false,
#   remote branches = [main]).
#
# WHAT IT DOES
#   1. preflight  — refuses to publish a hub that would not load (file inventory, bigBed
#                   magic number, per-file and total size, relative bigDataUrl values)
#   2. manifest   — writes trackhub/MANIFEST.json + trackhub/CHECKSUMS.md5, which
#                   12c_verify_trackhub_live.sh and the CI workflow check the live site against
#   3. stage      — a throwaway git worktree on an orphan branch, containing exactly
#                   .nojekyll, index.html and trackhub/
#   4. push       — ONLY with --push, and as an explicit force-push of a single commit
#
# USAGE
#   bash revision_G3/12b_publish_trackhub.sh                 # prepare + print the push command
#   bash revision_G3/12b_publish_trackhub.sh --push          # prepare and push (outward-facing)
#   bash revision_G3/12b_publish_trackhub.sh --preflight-only
#   bash revision_G3/12b_publish_trackhub.sh --append --push # keep gh-pages history
#
# NOTE ON HISTORY
#   The default is a single-commit orphan branch, force-pushed. The hub is a build product,
#   not a manuscript, and 105 MB per revision of history would be paid by every future clone
#   of the branch. --append keeps history if a versioned hub is ever wanted.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(dirname "$HERE")"
HUB_SRC="$HERE/trackhub"
BRANCH="gh-pages"
SITE_SUBDIR="trackhub"            # do not change: the manuscript prints .../trackhub/hub.txt
PAGES_BASE="https://nikit357.github.io/T2T_transposons_genes"
EXPECTED_BB=10
EXPECTED_HTML=11
MAX_FILE_MB=95                    # GitHub's hard per-file limit is 100 MB
MAX_TOTAL_MB=900                  # GitHub Pages' recommended site limit is 1 GB

DO_PUSH=0; APPEND=0; PREFLIGHT_ONLY=0; KEEP_WORKTREE=0
while [ $# -gt 0 ]; do
  case "$1" in
    --push)           DO_PUSH=1 ;;
    --append)         APPEND=1 ;;
    --preflight-only) PREFLIGHT_ONLY=1 ;;
    --keep-worktree)  KEEP_WORKTREE=1 ;;
    --hub-src)        HUB_SRC="$2"; shift ;;
    --branch)         BRANCH="$2"; shift ;;
    -h|--help)        sed -n '2,40p' "$0"; exit 0 ;;
    *) echo "unknown argument: $1 (try --help)" >&2; exit 2 ;;
  esac
  shift
done

fail() { echo "ERROR: $*" >&2; exit 1; }

echo "=== 1. preflight ==="
[ -d "$HUB_SRC" ] || fail "no hub at $HUB_SRC — run 12_build_trackhub.sh first"
for f in hub.txt genomes.txt hs1/trackDb.txt; do
  [ -s "$HUB_SRC/$f" ] || fail "missing or empty $HUB_SRC/$f"
done

n_bb=$(find "$HUB_SRC/hs1" -maxdepth 1 -name '*.bb' | wc -l)
n_html=$(find "$HUB_SRC/hs1" -maxdepth 1 -name '*.html' | wc -l)
[ "$n_bb" -eq "$EXPECTED_BB" ]     || fail "expected $EXPECTED_BB bigBed files, found $n_bb"
[ "$n_html" -eq "$EXPECTED_HTML" ] || fail "expected $EXPECTED_HTML description pages, found $n_html"
echo "  inventory: $n_bb bigBed, $n_html description pages, hub.txt, genomes.txt, trackDb.txt"

# Every bigBed must start with the little-endian bigBed magic number 0x8789F2EB. This catches a
# truncated copy, a Git-LFS pointer, or an HTML error page saved under a .bb name.
for bb in "$HUB_SRC"/hs1/*.bb; do
  magic="$(head -c 4 "$bb" | od -An -tx1 | tr -d ' \n')"
  [ "$magic" = "ebf28987" ] || fail "$(basename "$bb") is not a bigBed (magic $magic)"
done
echo "  every bigBed carries the bigBed magic number"

total=0
while IFS= read -r -d '' f; do
  bytes=$(stat -c%s "$f"); total=$((total + bytes))
  if [ "$bytes" -gt $((MAX_FILE_MB * 1000 * 1000)) ]; then
    fail "$(basename "$f") is $((bytes / 1000000)) MB, over the ${MAX_FILE_MB} MB guard"
  fi
done < <(find "$HUB_SRC" -type f -print0)
[ "$total" -le $((MAX_TOTAL_MB * 1000 * 1000)) ] \
  || fail "hub totals $((total / 1000000)) MB, over the ${MAX_TOTAL_MB} MB guard"
echo "  sizes: $((total / 1000000)) MB total, largest file $(du -m "$HUB_SRC"/hs1/*.bb | sort -rn | head -1 | cut -f1) MB"

# bigDataUrl must stay relative to hs1/trackDb.txt. An absolute URL here is how a hub silently
# keeps pointing at the wrong host after a move.
if grep -nE '^[[:space:]]*bigDataUrl[[:space:]]+https?://' "$HUB_SRC/hs1/trackDb.txt"; then
  fail "trackDb.txt has absolute bigDataUrl values; Pages publication needs relative ones"
fi
n_tracks=$(grep -cE '^[[:space:]]*bigDataUrl[[:space:]]' "$HUB_SRC/hs1/trackDb.txt")
[ "$n_tracks" -eq "$EXPECTED_BB" ] || fail "trackDb.txt declares $n_tracks bigDataUrl, expected $EXPECTED_BB"
echo "  trackDb.txt: $n_tracks relative bigDataUrl values"

echo
echo "=== 2. manifest ==="
# The manifest is what the live check compares Content-Length against, so a truncated upload or a
# gzip-mangled transfer is caught without downloading 105 MB.
manifest="$HUB_SRC/MANIFEST.json"
checks="$HUB_SRC/CHECKSUMS.md5"
( cd "$HUB_SRC" && find . -type f \
    ! -name MANIFEST.json ! -name CHECKSUMS.md5 -printf '%P\n' | LC_ALL=C sort \
  | xargs md5sum > "$checks" )
{
  printf '{\n  "hub": "T2T_transposons_genes",\n'
  printf '  "built_by": "revision_G3/12_build_trackhub.sh",\n'
  printf '  "published_by": "revision_G3/12b_publish_trackhub.sh",\n'
  printf '  "source_commit": "%s",\n' "$(git -C "$REPO" rev-parse HEAD)"
  printf '  "published_utc": "%s",\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf '  "base_url": "%s/%s",\n' "$PAGES_BASE" "$SITE_SUBDIR"
  printf '  "files": [\n'
  first=1
  while IFS= read -r rel; do
    [ $first -eq 1 ] || printf ',\n'; first=0
    printf '    {"path": "%s", "bytes": %s, "md5": "%s"}' \
      "$rel" "$(stat -c%s "$HUB_SRC/$rel")" "$(md5sum "$HUB_SRC/$rel" | cut -d' ' -f1)"
  done < <(awk '{print $2}' "$checks")
  printf '\n  ]\n}\n'
} > "$manifest"
echo "  MANIFEST.json and CHECKSUMS.md5 written ($(wc -l < "$checks") files)"

if [ "$PREFLIGHT_ONLY" -eq 1 ]; then
  echo; echo "preflight only — nothing staged."; exit 0
fi

echo
echo "=== 3. stage the gh-pages tree ==="
WORK="$(mktemp -d)"
WT="$WORK/worktree"
TMPBRANCH="publish-${BRANCH}-$(date +%y%m%d-%H%M%S)"
cleanup() {
  if [ "$KEEP_WORKTREE" -eq 1 ]; then
    echo "  worktree kept at $WT"
    return
  fi
  git -C "$REPO" worktree remove --force "$WT" >/dev/null 2>&1 || true
  git -C "$REPO" branch -D "$TMPBRANCH" >/dev/null 2>&1 || true
  rm -rf "$WORK"
}
trap cleanup EXIT

if [ "$APPEND" -eq 1 ]; then
  git -C "$REPO" fetch --quiet origin "$BRANCH" \
    || fail "--append needs an existing origin/$BRANCH"
  git -C "$REPO" worktree add --quiet -b "$TMPBRANCH" "$WT" "origin/$BRANCH"
  rm -rf "$WT/$SITE_SUBDIR"
else
  # A single-commit orphan branch. Verified working on git 2.43.0.
  git -C "$REPO" worktree add --quiet --orphan -b "$TMPBRANCH" "$WT"
fi

mkdir -p "$WT/$SITE_SUBDIR"
cp -a "$HUB_SRC/." "$WT/$SITE_SUBDIR/"

# .nojekyll: skip Jekyll entirely. Without it Pages runs a Jekyll build over 105 MB of binaries
# for no benefit, and Jekyll silently drops paths it does not like.
: > "$WT/.nojekyll"

cat > "$WT/index.html" <<'LANDING'
<!doctype html>
<meta charset="utf-8">
<title>T2T transposable elements and human genes — UCSC track hub</title>
<meta name="viewport" content="width=device-width, initial-scale=1">
<style>
  body { font: 16px/1.55 -apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif;
         max-width: 46rem; margin: 3rem auto; padding: 0 1.25rem; color: #1c1c1c; }
  code { background: #f2f2f2; padding: .1rem .3rem; border-radius: 3px; word-break: break-all; }
  table { border-collapse: collapse; margin: 1rem 0; }
  th, td { text-align: left; padding: .3rem .8rem .3rem 0; border-bottom: 1px solid #e6e6e6; }
  .cta { display: inline-block; padding: .55rem 1rem; background: #cc660b; color: #fff;
         text-decoration: none; border-radius: 4px; }
</style>
<h1>Transposable elements and human genes in T2T-CHM13v2.0</h1>
<p>UCSC Genome Browser track hub accompanying Nikitin (2026), <i>G3: Genes|Genomes|Genetics</i>.
Tracks are coloured with the same transposable-element class palette as the paper's figures.</p>

<p><a class="cta" href="https://genome.ucsc.edu/cgi-bin/hgTracks?db=hs1&amp;hubUrl=https://nikit357.github.io/T2T_transposons_genes/trackhub/hub.txt&amp;position=chr9:21150692-21370055">Open
the interferon-alpha domain in the UCSC Genome Browser</a></p>

<p>Hub URL, to paste into <i>My Data &rarr; Track Hubs &rarr; Connected Hubs</i>:<br>
<code>https://nikit357.github.io/T2T_transposons_genes/trackhub/hub.txt</code></p>

<h2>Tracks</h2>
<table>
<tr><th>Track</th><th>Contents</th></tr>
<tr><td>T2T TEs by class</td><td>All 3,709,429 RepeatMasker elements, one subtrack per class
(LINE, SINE, LTR, DNA, SVA, Helitron). Item name is the subfamily; item score is divergence in
substitutions per 1000 bp (observed range 0&ndash;480; higher means older).</td></tr>
<tr><td>TSS 10 kb windows</td><td>The 38,704 transcription start site neighbourhoods (5 kb either
side of each TSS) that define proximity throughout the paper.</td></tr>
<tr><td>TE-top / TE-bottom genes</td><td>The 1,436 genes with the highest and the lowest TE count in
their TSS neighbourhood, shown as their TSS windows.</td></tr>
<tr><td>IFNA domain</td><td>The 220 kb interferon-alpha domain, chr9:21,150,692&ndash;21,370,055.</td></tr>
</table>

<h2>Provenance</h2>
<p>Built by <code>revision_G3/12_build_trackhub.sh</code> and published by
<code>revision_G3/12b_publish_trackhub.sh</code> from
<a href="https://github.com/Nikit357/T2T_transposons_genes">github.com/Nikit357/T2T_transposons_genes</a>.
Per-file sizes and md5 checksums are in
<a href="trackhub/MANIFEST.json">trackhub/MANIFEST.json</a>.</p>
LANDING

git -C "$WT" add --all --force            # --force: trackhub/ is gitignored on the analysis branch
git -C "$WT" -c user.name="$(git -C "$REPO" config user.name)" \
             -c user.email="$(git -C "$REPO" config user.email)" \
  commit --quiet -m "Publish T2T TE track hub ($(date -u +%Y-%m-%d), source $(git -C "$REPO" rev-parse --short HEAD))"

staged_bb=$(git -C "$WT" ls-files "$SITE_SUBDIR/hs1/*.bb" | wc -l)
[ "$staged_bb" -eq "$EXPECTED_BB" ] || fail "staged $staged_bb bigBed files, expected $EXPECTED_BB"
# du of the worktree, not `git count-objects`: run from a worktree the latter reports the whole
# repository's pack size, which for this repo is ~600 MB and says nothing about what was staged.
echo "  staged $(git -C "$WT" ls-files | wc -l) files on $TMPBRANCH ($(du -sh --exclude=.git "$WT" | cut -f1))"
git -C "$WT" ls-files | sed 's/^/    /' | head -8
echo "    ..."

echo
echo "=== 4. push ==="
if [ "$DO_PUSH" -eq 1 ]; then
  if [ "$APPEND" -eq 1 ]; then
    git -C "$REPO" push origin "$TMPBRANCH:refs/heads/$BRANCH"
  else
    git -C "$REPO" push --force origin "$TMPBRANCH:refs/heads/$BRANCH"
  fi
  echo "  pushed to origin/$BRANCH"
  echo
  echo "NEXT (one time only, and only Daniil can do it): enable GitHub Pages"
  echo "  Settings -> Pages -> Build and deployment -> Source: Deploy from a branch"
  echo "  Branch: $BRANCH   Folder: / (root)   -> Save"
  echo "Then verify:"
  echo "  bash revision_G3/12c_verify_trackhub_live.sh"
else
  echo "  not pushed (no --push). To publish, either re-run with --push or run:"
  echo "    git -C $REPO push --force origin $TMPBRANCH:refs/heads/$BRANCH"
  echo "  (the staging worktree is removed when this script exits; use --keep-worktree to keep it)"
fi

echo
echo "one-click URL once Pages is live:"
echo "  https://genome.ucsc.edu/cgi-bin/hgTracks?db=hs1&hubUrl=$PAGES_BASE/$SITE_SUBDIR/hub.txt&position=chr9:21150692-21370055"
