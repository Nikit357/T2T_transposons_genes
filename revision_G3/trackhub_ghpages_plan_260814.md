# Track hub publication plan — fixing the 404 on `nikit357.github.io/T2T_transposons_genes/trackhub/hub.txt`

**Date:** 2026-08-14 · **Author:** Claude (plan only — nothing implemented) · **Target:** G3-2026-406828
revision, submission window closes 2026-08-28.

---

## 1. Overview

The Associate Editor asked for "a browser/track instance for the TE annotations made available with
the GitHub repository" (reviewer report, editor comment 3), and the revision built exactly that:
`revision_G3/12_build_trackhub.sh` produced a complete, `hubCheck`-clean hs1 track hub in
`revision_G3/trackhub/`. The manuscript, the README, `REPRODUCE.md`, the guidelines file and
`project_overview.md` all print the URL

```
https://nikit357.github.io/T2T_transposons_genes/trackhub/hub.txt
```

**That URL 404s because the hub was never published.** This is not a build defect and not a UCSC
problem: the `gh-pages` branch does not exist on the remote, and GitHub Pages has never been enabled
for the repository. The one item left open in `revision_G3/CLAUDE.md` — *"`gh-pages` published;
range-request check returns 206 — (needs Daniil)"* — is the whole cause.

The plan therefore does three things:

1. **Automates the publish step** that was previously five hand-typed commands in `REPRODUCE.md` §6,
   as `12b_publish_trackhub.sh`, with a preflight that refuses to publish a hub that would not load.
2. **Automates verification** of the live hub as `12c_verify_trackhub_live.sh`, and wires it into two
   GitHub Actions workflows (one manual-dispatch + weekly link-rot guard, one optional
   Actions-based Pages deployment).
3. **Lists the manual steps only Daniil can take** — the push, the one-time Pages enablement, and the
   browser confirmation — with the exact clicks and the API equivalents.

**Design decision:** publish to a `gh-pages` branch and serve it with GitHub Pages *from that branch*
(the classic "Deploy from a branch" source). This is the only option that keeps the URL already
printed in the manuscript byte-for-byte identical, which means **no manuscript edit is needed** if it
works. Every alternative changes the URL and costs a tracked edit to the Data availability statement.

---

## 2. Diagnosis — measured, not assumed

All of this was measured on 2026-08-14 from this machine.

| Check | Command | Result |
|---|---|---|
| Pages enabled for the repo? | `GET api.github.com/repos/Nikit357/T2T_transposons_genes` | **`has_pages: false`**, `private: false`, `default_branch: main` |
| Pages API record exists? | `GET …/repos/…/pages` | **404** (no Pages site configured) |
| Does `gh-pages` exist on the remote? | `git ls-remote --heads origin` | **only `refs/heads/main`** |
| Branches via API | `GET …/branches` | **`['main']`** |
| Live hub URL | `GET …/trackhub/hub.txt` | **404** |
| Live Pages root | `GET https://nikit357.github.io/T2T_transposons_genes/` | **404** |
| Local hub built? | `ls revision_G3/trackhub/` | **yes**, built 2026-08-04 |
| Local hub intact? | bigBed magic + `itemCount` read from each file's header | **10/10 valid** (`magic 0x8789F2EB`, version 4); LINE 1,005,214 · SINE 1,706,485 · LTR 531,410 · DNA 458,177 · SVA 6,274 · Helitron 1,869 · TSS windows 38,704 — every count matches the documented value |
| Sizes | `du`, `stat` | **109.0 MB total** (10 `.bb` + 11 `.html` + 3 text files); largest `TEs_SINE.bb` **42.5 MB**, well under GitHub's 100 MB hard limit |

**Conclusion: the artefact is fine, only the publication step is missing.** Nothing needs rebuilding
(and rebuilding would in fact be harder: `bedToBigBed`/`hubCheck`/`bigBedInfo` are no longer on
`PATH` on this machine, and the 155 MB `T2T_repeat_masker_processed_sorted.bed` input is gitignored).

### 2.1 Does GitHub Pages satisfy UCSC's byte-range requirement?

UCSC states hub files "must be located in web-accessible locations that support byte-range requests"
and lists GitHub among the acceptable hosts. Measured on GitHub-served static assets:

```
GET https://pages.github.com/images/slideshow/bootstrap.png   Range: bytes=100-199
  → HTTP/2 206 · accept-ranges: bytes · content-range: bytes 100-199/156915
GET https://raw.githubusercontent.com/Nikit357/T2T_transposons_genes/main/README.md  Range: 100-199
  → HTTP/2 206 · content-range: bytes 100-199/9696
```

So both GitHub Pages (binary asset) and `raw.githubusercontent.com` honour `Range` and return 206.
The one thing that cannot be tested before publishing is how Pages labels the unfamiliar `.bb`
extension: if it were served as a gzip-compressed text type, ranges would be meaningless. That is
exactly why `12c_verify_trackhub_live.sh` below does not just check for 206 — **it requests bytes
0–3 of every bigBed and asserts they equal the bigBed magic number `eb f2 89 87`**, which can only
be true if the byte offsets the browser asks for are the byte offsets it gets.

### 2.2 The editor's other complaint

The letter also says "the GitHub page link provided in the data availability statement is broken".
That is a separate, already-addressed item (edit K replaced the stale `T2T_genes_evolution` URL); the
repository URL `https://github.com/Nikit357/T2T_transposons_genes` resolves and the repo is public.
Item **V7** in §8 re-verifies it in the final docx so both halves of the editor's complaint are
answered in the response letter.

---

## 3. Constraints that fix the design

1. **The URL is already in the manuscript.** `13_manuscript_tracked_edits.py:1761` writes
   `https://nikit357.github.io/T2T_transposons_genes/trackhub/hub.txt` and the one-click
   `hgTracks?db=hs1&hubUrl=…&position=chr9:21150692-21370055` into the Data availability statement,
   and the same string is in `README.md:16-17`, `REPRODUCE.md:293,300`,
   `G3_article_guidelines.md:309`, `project_overview.md:389` and
   `12_build_trackhub.sh:34`. Therefore the published tree **must** have a top-level `trackhub/`
   directory: `gh-pages:/trackhub/hub.txt`, not `gh-pages:/hub.txt`.
2. **`trackDb.txt` uses relative `bigDataUrl` values** (`TEs_LINE.bb`, resolved against the location
   of `trackDb.txt`). That is correct for Pages and must stay relative — it is also what makes the
   `raw.githubusercontent` fallback a one-command rewrite rather than a rebuild.
3. **CI cannot rebuild the hub.** `12a_trackhub_beds.py` needs the 155 MB gitignored RepeatMasker BED
   and UCSC's `bedToBigBed`. So the Actions workflows **verify and deploy, they never build**. This
   is a deliberate limitation, stated in the workflow header comment so nobody "fixes" it later.
4. **The working tree on `main` is dirty** — `.gitignore` modified and 16 root-level
   `Supplementary File/Figure *` files deleted but not committed. The publish path must not touch it,
   so `12b` works in a **separate `git worktree`** under `$TMPDIR` and never checks out a branch in
   the main working directory.
5. **`revision_G3/trackhub/` stays gitignored on `main`** (`.gitignore:24`). Publishing copies out of
   an ignored directory, so the staging step uses `git add --all --force`.
6. **No frozen file is touched.** The four md5-frozen files are irrelevant here; `12_build_trackhub.sh`
   is not frozen but gets only an additive final message plus one optional label wording change.

---

## 4. Files to create

### 4.1 NEW — `revision_G3/12b_publish_trackhub.sh`

Everything before the network: preflight, staging, `.nojekyll`, landing page, manifest, commit.
**It does not push unless `--push` is given**, per the Human-in-the-Loop rule in `CLAUDE.md`.

```bash
#!/usr/bin/env bash
# Publish the built UCSC track hub to the gh-pages branch (WP13b, decision D8).
#
# WHY THIS EXISTS
#   12_build_trackhub.sh builds revision_G3/trackhub/ (109 MB), which is gitignored on the
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
#   not a manuscript, and 109 MB per revision of history would be paid by every future clone
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
# gzip-mangled transfer is caught without downloading 109 MB.
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

# .nojekyll: skip Jekyll entirely. Without it Pages runs a Jekyll build over 109 MB of binaries
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
echo "  staged $(git -C "$WT" ls-files | wc -l) files on $TMPBRANCH ($(git -C "$WT" count-objects -vH | awk '/size-pack|size:/{print $2 $3}' | head -1))"
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
```

### 4.2 NEW — `revision_G3/12c_verify_trackhub_live.sh`

The check that proves the hub is usable, not merely reachable. Runs locally and unchanged in CI.

```bash
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
```

### 4.3 NEW — `.github/workflows/verify-trackhub.yml`

The repository has **no** `.github/` directory today; this creates it. The workflow is the link-rot
guard the journal cares about: the URL is printed in a paper, so it should be checked without anyone
remembering to.

```yaml
# Verify the published UCSC track hub is reachable AND range-readable.
#
# This workflow deliberately does NOT build the hub. Building needs the 155 MB gitignored
# T2T_repeat_masker_processed_sorted.bed and UCSC's bedToBigBed, neither of which belongs in CI.
# Build locally with revision_G3/12_build_trackhub.sh, publish with 12b, verify here.
name: verify-trackhub

on:
  workflow_dispatch:
    inputs:
      base_url:
        description: "Hub base URL to check"
        required: false
        default: "https://nikit357.github.io/T2T_transposons_genes/trackhub"
  schedule:
    # Weekly, Monday 06:17 UTC. The URL is cited in a published paper; a silent 404 is a defect.
    - cron: "17 6 * * 1"

permissions:
  contents: read

jobs:
  verify:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4

      - name: Install UCSC hubCheck
        run: |
          curl -sSfLO https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/hubCheck
          chmod +x hubCheck
          echo "$PWD" >> "$GITHUB_PATH"

      - name: Verify the live hub
        run: |
          bash revision_G3/12c_verify_trackhub_live.sh \
            --base-url "${{ inputs.base_url || 'https://nikit357.github.io/T2T_transposons_genes/trackhub' }}" \
            --remote-manifest \
            --json trackhub_live_check.json

      - name: Upload the report
        if: always()
        uses: actions/upload-artifact@v4
        with:
          name: trackhub-live-check
          path: trackhub_live_check.json
```

`--remote-manifest` is why CI needs no hub in the repository: the manifest `12b` published alongside
the hub carries the expected byte counts, so a truncated or transcoded file is still caught.

### 4.4 NEW — `.github/workflows/deploy-trackhub.yml`

**Only needed if Pages is set to the "GitHub Actions" source instead of "Deploy from a branch".**
Recommended path is a branch source (§5), which needs no workflow at all — GitHub runs its own
`pages-build-deployment`. This file is the alternative, kept dormant behind `workflow_dispatch` so it
cannot fight the branch deployment.

```yaml
# Deploy the gh-pages branch contents to GitHub Pages via Actions.
#
# Use this ONLY if Settings -> Pages -> Source is set to "GitHub Actions". If the source is
# "Deploy from a branch" (the recommended setup), GitHub deploys gh-pages by itself and this
# workflow must stay unused — two deployment paths for one site is how a site goes stale.
#
# It does not build anything: it publishes what 12b_publish_trackhub.sh pushed to gh-pages.
name: deploy-trackhub

on:
  workflow_dispatch:
    inputs:
      ref:
        description: "Branch holding the hub"
        required: false
        default: "gh-pages"

permissions:
  contents: read
  pages: write
  id-token: write

concurrency:
  group: pages
  cancel-in-progress: false

jobs:
  deploy:
    runs-on: ubuntu-latest
    environment:
      name: github-pages
      url: ${{ steps.deployment.outputs.page_url }}
    steps:
      - uses: actions/checkout@v4
        with:
          ref: ${{ inputs.ref || 'gh-pages' }}

      - name: Sanity-check the tree before deploying
        run: |
          test -s trackhub/hub.txt
          test -s trackhub/genomes.txt
          test -s trackhub/hs1/trackDb.txt
          n=$(find trackhub/hs1 -name '*.bb' | wc -l)
          [ "$n" -eq 10 ] || { echo "expected 10 bigBed files, found $n"; exit 1; }
          for bb in trackhub/hs1/*.bb; do
            magic=$(head -c 4 "$bb" | od -An -tx1 | tr -d ' \n')
            [ "$magic" = "ebf28987" ] || { echo "$bb is not a bigBed"; exit 1; }
          done
          du -sh trackhub

      - uses: actions/configure-pages@v5
      - uses: actions/upload-pages-artifact@v3
        with:
          path: .
      - id: deployment
        uses: actions/deploy-pages@v4

      - name: Verify the deployment
        run: |
          curl -sSfL -o /dev/null "${{ steps.deployment.outputs.page_url }}trackhub/hub.txt"
          echo "hub.txt reachable at ${{ steps.deployment.outputs.page_url }}trackhub/hub.txt"
```

---

## 5. Files to change

### 5.1 `revision_G3/12_build_trackhub.sh` — one additive message

The build script currently ends by printing the one-click URL, which is misleading while nothing is
published. **Before** (lines 240–244):

```bash
echo "=== 5. one-click URL ==="
echo "https://genome.ucsc.edu/cgi-bin/hgTracks?db=hs1&hubUrl=$HUB_URL_BASE/hub.txt&position=chr9:21150692-21370055"
echo
du -sh "$OUT"
echo "built -> $OUT"
```

**After:**

```bash
echo "=== 5. one-click URL (works only once the hub is published) ==="
echo "https://genome.ucsc.edu/cgi-bin/hgTracks?db=hs1&hubUrl=$HUB_URL_BASE/hub.txt&position=chr9:21150692-21370055"
echo
du -sh "$OUT"
echo "built -> $OUT"
echo
echo "=== 6. publish ==="
echo "  This directory is gitignored on the analysis branch. Nothing is public until it is"
echo "  pushed to gh-pages and GitHub Pages is enabled:"
echo "    bash revision_G3/12b_publish_trackhub.sh            # preflight + stage, no push"
echo "    bash revision_G3/12b_publish_trackhub.sh --push      # publish"
echo "    bash revision_G3/12c_verify_trackhub_live.sh         # 200 + 206 + magic-number check"
```

### 5.2 `revision_G3/12_build_trackhub.sh` — optional accuracy fix to two `longLabel`s

Measured: `genes_TE_top.bb` holds **2,075** intervals and `genes_TE_bottom.bb` **2,271**, because
`12a_trackhub_beds.py` draws each gene set as its TSS windows and a gene with several annotated TSS
contributes several windows (the per-TSS property documented throughout the project). The labels say
"The 1,436 genes with the highest TE count", so a reviewer who counts features sees a number that
does not match the label. Optional, and independent of the 404 fix.

**Before** (lines 122–140, two `longLabel` lines):

```
longLabel The 1,436 genes with the highest TE count in their TSS 10 kb neighbourhood (top 5 percent)
...
longLabel The 1,436 genes with the lowest TE count in their TSS 10 kb neighbourhood (bottom 5 percent)
```

**After:**

```
longLabel The 1,436 genes with the highest TE count in their TSS 10 kb neighbourhood (top 5 percent), shown as their 2,075 TSS windows
...
longLabel The 1,436 genes with the lowest TE count in their TSS 10 kb neighbourhood (bottom 5 percent), shown as their 2,271 TSS windows
```

`trackDb.txt` is regenerated only by a full `12_build_trackhub.sh` run, which needs `bedToBigBed`.
To take this fix without a rebuild, apply the same wording to the built file before publishing:

```bash
sed -i 's/(top 5 percent)$/(top 5 percent), shown as their 2,075 TSS windows/;
        s/(bottom 5 percent)$/(bottom 5 percent), shown as their 2,271 TSS windows/' \
  revision_G3/trackhub/hs1/trackDb.txt
```

### 5.3 `REPRODUCE.md` §6 — replace the hand-typed publish recipe

**Before** (lines 277–295): the `git switch --orphan gh-pages` / `cp -r` / `git push` block and the
lone `curl -sI -H 'Range: …'` check.

**After:**

```markdown
The hub is **not** tracked on the main branch — it is published to `gh-pages`, so that a
clone does not carry 109 MB twice. Publishing and verifying are scripted:

```bash
bash revision_G3/12b_publish_trackhub.sh            # preflight + stage; prints the push command
bash revision_G3/12b_publish_trackhub.sh --push     # publish to origin/gh-pages
```

`12b` refuses to publish a hub that would not load: it checks the file inventory, reads the
bigBed magic number out of every `.bb`, enforces the 100 MB per-file and 1 GB per-site limits,
insists that `bigDataUrl` values stay relative, and writes `MANIFEST.json` + `CHECKSUMS.md5`
next to the hub.

**One manual step, once, and only the repository owner can do it.** GitHub Pages has to be
switched on: *Settings → Pages → Build and deployment → Source: Deploy from a branch →
Branch `gh-pages`, folder `/ (root)` → Save*. Equivalent API call, with a token carrying
`administration:write` and `pages:write`:

```bash
curl -X POST -H "Authorization: Bearer $GH_TOKEN" \
     -H "Accept: application/vnd.github+json" \
     https://api.github.com/repos/Nikit357/T2T_transposons_genes/pages \
     -d '{"source":{"branch":"gh-pages","path":"/"}}'
```

The first build takes a minute or two; watch it under the repository's Actions tab as
`pages-build-deployment`. Then verify:

```bash
bash revision_G3/12c_verify_trackhub_live.sh
```

This is stronger than a status check. For every bigBed it requests bytes 0–3 and asserts they
are the bigBed magic number `eb f2 89 87`, so it proves the host honours byte ranges without
transcoding the payload — which is the property UCSC actually requires, and the only thing that
cannot be tested before publication. The same script runs weekly in
`.github/workflows/verify-trackhub.yml`, so a future 404 is caught rather than reported by a
reader.
```

(Keep the existing one-click URL block that follows.)

### 5.4 `revision_G3/CLAUDE.md`

- In **"The UCSC track hub"**: after the existing sentence, add — *"`12b_publish_trackhub.sh`
  publishes it to `gh-pages` (preflight, `.nojekyll`, landing page, `MANIFEST.json`, single-commit
  force-push; it does not push without `--push`) and `12c_verify_trackhub_live.sh` verifies the live
  site with a 200 + 206 + bigBed-magic-number check that also runs weekly in
  `.github/workflows/verify-trackhub.yml`. **Publishing is what was missing until 2026-08-14: the
  branch had never been pushed and Pages had never been enabled, so the URL printed in the
  manuscript 404'd.** The published tree must keep `trackhub/` as a top-level directory, because the
  manuscript prints `.../trackhub/hub.txt`."*
- In **"Directory layout"**, after the `12_build_trackhub.sh + 12a_…` line:
  ```
  12b_publish_trackhub.sh             publish trackhub/ to gh-pages; --push gates the network step
  12c_verify_trackhub_live.sh         live hub check: 200, 206, bigBed magic, sizes vs MANIFEST.json
  .github/workflows/                  verify-trackhub.yml (dispatch + weekly), deploy-trackhub.yml (dormant)
  ```
- In **"Run order"**, Phase 6:
  ```bash
  bash revision_G3/12_build_trackhub.sh                   # bigBeds + hub + hubCheck
  bash revision_G3/12b_publish_trackhub.sh --push         # -> origin/gh-pages (needs credentials)
  bash revision_G3/12c_verify_trackhub_live.sh            # after Pages is enabled
  ```
- In **"What is still open"**: replace *"publishing `trackhub/` to `gh-pages`"* with *"run
  `12b_publish_trackhub.sh --push`, enable Pages on `gh-pages` (one click, owner only), then
  `12c_verify_trackhub_live.sh`"*.

### 5.5 `CLAUDE.md` (parent, `T2T_transposons_genes/`)

In **"The UCSC track hub"**, replace *"It is gitignored on this branch and published to `gh-pages`"*
with *"It is gitignored on this branch and published to `gh-pages` by
`revision_G3/12b_publish_trackhub.sh`, then checked live by `12c_verify_trackhub_live.sh`. GitHub
Pages must be enabled once by hand on the `gh-pages` branch; until that is done the URL the
manuscript prints returns 404."*

### 5.6 `revision_G3/README.md` and `revision_G3/project_overview.md`

- `README.md`: add `12b`/`12c` and `.github/workflows/` to the file inventory; in the deliverables
  section note that the hub is live only after Pages is enabled.
- `project_overview.md` §8 (open items): replace the "publish `trackhub/` to `gh-pages`" bullet with
  the three-step sequence, and record the 2026-08-14 diagnosis (`has_pages: false`, remote branches
  `[main]`) so the cause is not re-diagnosed.

### 5.7 `README.md` (repository root)

Under **Genome browser**, add one line after the hub URL:

```markdown
  Landing page: <https://nikit357.github.io/T2T_transposons_genes/>. If the hub fails to load,
  it is verified weekly by [`.github/workflows/verify-trackhub.yml`](.github/workflows/verify-trackhub.yml);
  the last report is the workflow's `trackhub-live-check` artifact.
```

### 5.8 Files that do **not** change, and why

| File | Why not |
|---|---|
| `revision_G3/12a_trackhub_beds.py` | The hub data is correct and verified; nothing about publication changes how the BEDs are built |
| The four md5-frozen files | Untouched, as always |
| `Revised_manuscript/*.docx` | **No manuscript edit is needed if Pages works**, because the plan preserves the exact URL already in the Data availability statement. Only the fallbacks in §7 would require one |
| `.gitignore` | `revision_G3/trackhub/` stays ignored on the analysis branch — that is the point of publishing on a separate branch. `12b` stages with `git add --force` inside its own worktree |
| `13_manuscript_tracked_edits.py` | Its URL string is correct; changing it would rewrite the docx for no reason |
| `11_results_numbers.py`, `08`/`17` supplementary builders | No number changes |

---

## 6. Manual steps — what only Daniil can do

| # | Step | Command / clicks | Why it cannot be automated here |
|---|---|---|---|
| **M1** | Approve this plan | — | Human-in-the-Loop rule |
| **M2** | Check push credentials | `git push --dry-run origin HEAD` | A global `credential.helper store` is configured on this machine, but whether it holds a token for `github.com` with write scope is not something I should read. If it fails, create a fine-grained PAT with *Contents: read and write* (and *Pages: read and write* + *Administration: read and write* if you want to enable Pages via API) |
| **M3** | Commit and push the new scripts + workflows on `main` | `git add .github revision_G3/12b_publish_trackhub.sh revision_G3/12c_verify_trackhub_live.sh revision_G3/12_build_trackhub.sh REPRODUCE.md README.md CLAUDE.md revision_G3/CLAUDE.md revision_G3/README.md revision_G3/project_overview.md revision_G3/trackhub_ghpages_plan_260814.md && git commit && git push origin main` | **Do not use `git add -A`.** The working tree also has `.gitignore` modified and 16 root-level `Supplementary File/Figure *` files deleted but not committed — unrelated pending changes that need their own decision. Scheduled and dispatchable workflows are only picked up from the default branch, so this push is what makes the Actions tab show them |
| **M4** | Publish the hub | `bash revision_G3/12b_publish_trackhub.sh` then, after reading the staging report, `bash revision_G3/12b_publish_trackhub.sh --push` | Pushing 109 MB to a public branch is outward-facing; the script refuses to do it without the explicit flag |
| **M5** | **Enable GitHub Pages (the actual fix for the 404)** | *Settings → Pages → Build and deployment → Source: **Deploy from a branch** → Branch **`gh-pages`**, Folder **`/ (root)`** → Save*. API equivalent in §5.3 | Requires repository-admin rights. Note the error text you saw — *"There isn't a GitHub Pages site here"* — is exactly what an un-enabled Pages site returns; no amount of pushing fixes it alone |
| **M6** | Watch the first build | Repository → **Actions** → `pages-build-deployment` → green | — |
| **M7** | Verify | `bash revision_G3/12c_verify_trackhub_live.sh` — expect `pass` on every line and exit 0. Also run it from the Actions tab: **Actions → verify-trackhub → Run workflow** | — |
| **M8** | Confirm in the browser | Open the one-click URL. Check that the six class subtracks, the TSS windows, the two gene-set tracks and the IFNA domain all draw at chr9:21,150,692–21,370,055, and that clicking a track name shows its description page. Screenshot it for the response letter | Visual confirmation |
| **M9** | If the hub was ever loaded while broken | UCSC caches remote hub files for a few minutes. If a stale error persists, disconnect the hub under *My Data → Track Hubs → Connected Hubs* and reconnect | — |
| **M10** | Response letter | Add to the editor-comment reply: the hub is live at the URL in the Data availability statement, the one-click link opens at the IFNA domain, it is built by `12_build_trackhub.sh`, published by `12b`, and re-verified weekly by CI. Also confirm the repository link in the Data availability statement resolves (the letter flagged it as broken; the current URL is `https://github.com/Nikit357/T2T_transposons_genes`, which is public and resolves) | Author's text |
| **M11** | Optional | Register the hub with UCSC's public hub list (`genome.ucsc.edu` → Track Hubs → *Contact us* to propose a public hub) after publication | Editorial decision, not needed for acceptance |

---

## 7. If GitHub Pages does not work — fallbacks, in order

Each fallback is cheap, and each one after the first costs a manuscript edit.

**P1 (primary, §4–§6): Pages from the `gh-pages` branch.** Keeps the manuscript URL; no docx edit.

**P2: Pages via Actions** (`deploy-trackhub.yml`, §4.4). Same URL, same result; use only if the
branch source is unavailable (for example if the account enforces the Actions source). Set
*Settings → Pages → Source: GitHub Actions*, then **Actions → deploy-trackhub → Run workflow**.

**P3: `raw.githubusercontent.com`.** UCSC documents this explicitly — *"GitHub supports byte-range
access to files when they are accessed via the `raw.githubusercontent.com` style URLs"* — and it
measured 206 above. Requires rewriting `bigDataUrl` from relative to absolute and one tracked
manuscript edit:

```bash
RAW=https://raw.githubusercontent.com/Nikit357/T2T_transposons_genes/gh-pages/trackhub/hs1
sed -i "s#^\([[:space:]]*\)bigDataUrl \(.*\.bb\)#\1bigDataUrl $RAW/\2#" \
  revision_G3/trackhub/hs1/trackDb.txt
# hub URL becomes:
#   https://raw.githubusercontent.com/Nikit357/T2T_transposons_genes/gh-pages/trackhub/hub.txt
```

Caveats: the per-track `html` description pages do not render as HTML from `raw`, and the URL
embeds a branch name, so it breaks if the branch is renamed. Only take this if P1 and P2 both fail.

**P4: UCSC Hub Space** — UCSC's own hosting, 10 GB free per user, at
`genome.ucsc.edu/cgi-bin/hgHubConnect#hubUpload`. The most robust option for a published paper
because the host is the browser itself, but the URL is UCSC-issued, so it needs the same one-line
docx edit as P3.

**P5 (not recommended): a separate repository** for the hub, to keep the analysis repository small.
Changes the URL, splits the deliverable from the code, and adds a second thing to maintain.

**If a manuscript edit becomes necessary**, do it as a new `20_trackhub_url_edit.py` in the
established shape: read `T2T_genes_article_G3_revision_260810.docx`, write a new dated file, operate
on whole `<w:sdt>` elements, assert the content-control count before and after, and set
`wrt._rid[0]` above the highest existing revision id (the four-pass duplicate-id trap). Two strings
change — the hub URL and the `hubUrl=` inside the one-click link — plus the same two in
`README.md`, `REPRODUCE.md`, `12_build_trackhub.sh:34`, `G3_article_guidelines.md:309` and
`project_overview.md:389`.

---

## 8. Side effects and caveats

1. **Repository size.** `.git` is already 565 MB; a 109 MB orphan branch takes it to roughly 675 MB.
   Still inside GitHub's recommended 1 GB, but each republished single-commit branch leaves the
   previous 109 MB unreachable until GitHub's own GC prunes it. Republish only when the hub data
   actually changes.
2. **Force-push by default.** `12b` without `--append` replaces `gh-pages` with a single commit.
   That is right for a build product and wrong for anything else — it is why the branch name is a
   parameter and the push is a flag.
3. **`.nojekyll` matters.** Without it Pages runs Jekyll over 109 MB of binaries; with it the tree
   is copied verbatim. It is written by `12b`, not by hand.
4. **The `trackhub/` sub-directory is load-bearing.** Publishing the hub at the site root would give
   `…/T2T_transposons_genes/hub.txt` and invalidate the URL already in the manuscript. `12b` hard-codes
   `SITE_SUBDIR=trackhub` and the comment says why.
5. **CI cannot rebuild the hub**, by design (§3.3). Anyone adding a build step to the workflow will
   need the gitignored 155 MB input and will be tempted to re-download RepeatMasker — which would
   silently change the annotation the paper reports.
6. **Pages propagation.** The site is not instant: the first build takes a minute or two, and the
   CDN caches responses for about 10 minutes, so run `12c` again before concluding something is
   wrong. UCSC additionally caches remote hub files briefly (M9).
7. **Bandwidth.** GitHub Pages has a soft 100 GB/month limit. A bigBed hub streams only the ranges a
   viewer looks at, so realistic reviewer traffic is megabytes.
8. **The 16 uncommitted supplementary deletions and the modified `.gitignore` are untouched** by this
   plan and still need their own decision (M3).
9. **No Figma work, no manuscript rebuild, no notebook execution** — this plan changes nothing that
   feeds a figure or a number, so `11_results_numbers.py` output is unaffected.

---

## 9. Verification commands

```bash
# --- local, before publishing -------------------------------------------------
bash revision_G3/12b_publish_trackhub.sh --preflight-only
#   expect: 10 bigBed, 11 description pages, magic ok, ~109 MB total, 10 relative bigDataUrl

# --- after M4 (push) but before M5 (Pages enabled) ----------------------------
git ls-remote --heads origin | grep gh-pages          # expect one line
curl -s -o /dev/null -w '%{http_code}\n' \
  https://nikit357.github.io/T2T_transposons_genes/trackhub/hub.txt
#   expect: still 404 — pushing does not enable Pages. This is the current bug, isolated.

# --- after M5 ------------------------------------------------------------------
curl -s https://api.github.com/repos/Nikit357/T2T_transposons_genes \
  | grep '"has_pages"'                                 # expect: true
bash revision_G3/12c_verify_trackhub_live.sh
#   expect: PASS on hub.txt, genomes.txt, hs1/trackDb.txt, 10 x (206 + magic + size),
#           11 description pages, site root; "0 failed"; exit 0

# the single decisive check, by hand
curl -s -r 0-3 \
  https://nikit357.github.io/T2T_transposons_genes/trackhub/hs1/TEs_LINE.bb | od -An -tx1
#   expect: eb f2 89 87

# --- CI ------------------------------------------------------------------------
# Actions -> verify-trackhub -> Run workflow  (green, artifact trackhub-live-check)
```

Expected `output/trackhub_live_check.json`: `"passed": 25, "failed": 0` (3 text + 10 bigBed +
11 description pages + site root, plus `hubCheck` if installed).

---

## 10. Rollback

Nothing in this plan is hard to undo. `git push --delete origin gh-pages` removes the branch;
*Settings → Pages → Source → None* unpublishes the site; deleting `.github/workflows/*.yml` removes
the CI. The analysis branch is never checked out or modified by the publish path, so a failed publish
cannot damage the working tree.

---

## 11. TODO

Implementation (mine, after approval):

- [x] `revision_G3/12b_publish_trackhub.sh` — new, per §4.1; `chmod +x`
- [x] `revision_G3/12c_verify_trackhub_live.sh` — new, per §4.2; `chmod +x`
- [x] `.github/workflows/verify-trackhub.yml` — new, per §4.3
- [x] `.github/workflows/deploy-trackhub.yml` — new, per §4.4 (dormant; dispatch only)
- [x] `revision_G3/12_build_trackhub.sh` — add the "=== 6. publish ===" block (§5.1)
- [x] `revision_G3/12_build_trackhub.sh` — *optional* gene-set `longLabel` wording (§5.2), plus the
      `sed` on the built `trackDb.txt` so the published hub matches without a rebuild
      (taken; 2,075 / 2,271 re-verified from the bigBed headers before editing the labels)
- [x] `REPRODUCE.md` §6 — replace the hand-typed publish/verify recipe (§5.3)
- [x] `revision_G3/CLAUDE.md` — track hub section, directory layout, run order, "what is still open" (§5.4)
      (the plan assumed a "The UCSC track hub" section already existed here; it did not — that section
      is in the parent `CLAUDE.md` — so the paragraph was added as a **new** section after the
      directory layout)
- [x] `CLAUDE.md` (parent) — track hub paragraph (§5.5)
- [x] `revision_G3/README.md` — file inventory + deliverable note (§5.6)
- [x] `revision_G3/project_overview.md` §8 — open items + the 2026-08-14 diagnosis (§5.6)
- [x] `README.md` (root) — landing page + weekly-check line (§5.7)
- [x] Run `bash revision_G3/12b_publish_trackhub.sh --preflight-only` and paste the report here —
      passed 2026-08-16: 10 bigBed, 11 description pages, magic ok on all ten, 109 MB total
      (largest 41 MB), 10 relative `bigDataUrl`, `MANIFEST.json` + `CHECKSUMS.md5` over 24 files
- [x] `bash -n` both new scripts, and `python3 -c "import yaml, sys; [yaml.safe_load(open(f)) for f in sys.argv[1:]]"`
      on both workflow files — both clean
- [x] Staging exercised end-to-end without `--push`: orphan worktree on git 2.43.0, 28 files staged
      (105 MB), commit made, worktree and temp branch cleaned up, `git worktree list` back to one
      entry. One fix on top of §4.1: the staged-size line used `git count-objects`, which from a
      worktree reports the **whole repository's** pack (633 MiB) rather than the staged tree, so it
      now reports `du` of the worktree
- [x] Diagnosis re-confirmed 2026-08-16, unchanged from 2026-08-14: `git ls-remote --heads origin`
      has no `gh-pages`, and `hub.txt` returns **404**

Daniil's manual steps, in order (§6):

- [x] **M1** approve this plan — approved 2026-08-16
- [ ] **M2** `git push --dry-run origin HEAD` — confirm write credentials
- [ ] **M3** commit and push the new scripts + workflows on `main` (selective `git add`, **not** `-A`)
- [ ] **M4** `bash revision_G3/12b_publish_trackhub.sh` then `--push`
- [ ] **M5** Settings → Pages → Deploy from a branch → `gh-pages` / `/ (root)` → Save
- [ ] **M6** Actions → `pages-build-deployment` green
- [ ] **M7** `bash revision_G3/12c_verify_trackhub_live.sh` → 0 failed; also dispatch `verify-trackhub`
- [ ] **M8** open the one-click URL, confirm all ten tracks draw, screenshot for the response letter
- [ ] **M9** if a stale error shows, reconnect the hub in *My Data → Track Hubs*
- [ ] **M10** response-letter sentence for the editor comment + confirm the repository link in the
      Data availability statement resolves
- [ ] **M11** *(optional)* propose the hub to UCSC's public hub list after publication
