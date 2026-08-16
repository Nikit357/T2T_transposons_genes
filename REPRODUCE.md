# REPRODUCE.md — running this analysis from a clean checkout

Everything in the paper is reproducible from this repository plus four downloads. This file
gives the exact order, the resolved paths, and what to do about the files the repository
deliberately does not carry.

Read it in order the first time. The relative paths inside the original notebooks are
inconsistent because they were written across several sessions with different working
directories (`CLAUDE.md` documents which ones are wrong); every path below is given relative
to the repository root and has been checked.

---

## 0. What you need before you start

| Requirement | Version used | Notes |
|---|---|---|
| Python | 3.11 | `pandas`, `numpy`, `scipy`, `statsmodels`, `goatools`, `matplotlib`, `seaborn`, `networkx`, `pyvis`, `supervenn`, `plotly`, `adjustText`, `statannotations`, `scikit-learn`, `tqdm`, `zstandard` |
| Python (manuscript only) | 3.11 + `python-docx` + `lxml` | a **separate** environment; see the warning below |
| `bedtools` | 2.31.1 | on `PATH` |
| `zstd` | any | the permutation store is zstd-compressed |
| `pandoc` | any | only for producing a plain-text view of the manuscript |
| UCSC `bedToBigBed`, `hubCheck`, `bigBedInfo` | — | track hub only; see step 6 |
| Disk | ~15 GB free | the permutation background is regenerated into a temporary directory |

**The two environments are not interchangeable.** The analysis environment does **not** have
`python-docx`, and the manuscript scripts need it. On the machine this was developed on:

```bash
source ~/venvs/Retroelements_3_11/bin/activate     # analysis: steps 1-5, 7
~/venvs/collagen_3_11/bin/python <script>          # manuscript: step 8
```

Substitute your own environments; the point is that steps 1–7 and step 8 have different
dependency sets.

---

## 1. Inputs the repository does not carry, and how to get them

Four files are excluded from git. Each is either over GitHub's 100 MB per-file hard limit or
is a large third-party download. None of them is a result — they are all inputs or
regenerable intermediates, and every table needed to *check* a number in the paper is
tracked.

### 1.1 `T2T_repeat_masker_processed_sorted.bed` (155 MB) — required by almost everything

The processed RepeatMasker annotation: 3,709,429 rows, tab-separated, columns
`chrom, start, end, divergence, subfamily, family, class`.

Two ways to get it:

- **Copy it** from the sibling project, which owns the master table:
  `../../T2T_article/T2T_repeat_masker_processed.csv` — convert to the seven BED columns above
  and `sort -k1,1 -k2,2n`.
- **Rebuild it** from UCSC: Table Browser → *assembly* `Jan 2022 (T2T CHM13v2.0/hs1)`,
  *group* `Variation and Repeats`, *track* `RepeatMasker`, *table* `T2T RepeatMasker`, output
  as all fields. Keep `genoName, genoStart, genoEnd, milliDiv, repName, repFamily, repClass`
  in that order, then sort as above.

Sanity check once you have it:

```bash
wc -l T2T_repeat_masker_processed_sorted.bed        # 3709429
cut -f7 T2T_repeat_masker_processed_sorted.bed | sort | uniq -c | sort -rn
# 1706485 SINE / 1005214 LINE / 531410 LTR / 458177 DNA / 6274 Retroposon / 1869 RC
```

`Retroposon` is SVA and `RC` is Helitron; the paper uses the second name of each pair, the
files use the first.

### 1.2 `goa_human.gaf` (190 MB) — required for any Gene Ontology step

```bash
wget https://current.geneontology.org/annotations/goa_human.gaf.gz
gunzip goa_human.gaf.gz
```

The version used was downloaded on 2025-12-31. `go-basic.obo` **is** tracked, so the ontology
itself is pinned; only the annotations need downloading.

### 1.3 `revision_G3/external/HOM_MouseHumanSequence.rpt` (15 MB) — WP4 only

MGI mouse–human homology, needed only for the comparison with Lu et al. 2020. Exact URL,
size and md5 are in `revision_G3/external/PROVENANCE.md`, which also records the four
download routes that no longer work and why the comparison needs orthology mapping at all
(their published categories are **mouse** gene sets).

### 1.4 The 5 kb and 20 kb permutation backgrounds — window sensitivity only

`revision_G3/output/permutation_counts_5kb/` (81 MB) and `_20kb/` (135 MB) are excluded.
The **10 kb** store, which every headline number in the paper rests on, **is** tracked at
`revision_G3/output/permutation_counts_10kb/` (106 MB, 500 seeds plus `MANIFEST.json`).
Regenerate the other two with step 3 below (2–4 h each). Every table derived from them is
tracked, so no claim in the paper is unverifiable without them.

---

## 2. The published analysis (frozen)

These three notebooks and one script produced the originally submitted results. **They are
frozen for the G3 revision and must not be edited** — `CLAUDE.md` records their MD5s and the
reason. Run them only if you want to regenerate the published outputs from scratch; every
output they produce is already tracked.

```bash
jupyter nbconvert --execute download_and_process_files_UCSC_genes.ipynb   # annotations, windows
jupyter nbconvert --execute TEs_mapped_on_TSS_analysis.ipynb              # proximity + enrichment
jupyter nbconvert --execute Gene_ontology_analysis.ipynb                  # GO + networks
python GO_subfamilies.py                                                  # per-subfamily GO
```

Two things to know before you run them:

- **The permutation count.** Cell 34 of the download notebook sets `NUM_PERMUTATIONS = 1000`.
  That line was **never executed at 1,000**. What actually ran was **N = 500**: there are 501
  files in the original `permutation_results/`, the empirical p-value floor is
  2/501 = 0.004, and the revised manuscript says 500 throughout. The line stays in the
  notebook because the notebook is frozen; do not treat it as authoritative. See caveat C20
  in `G3_revision_implementation_plan_260803.md`.
- **A broken path.** `TEs_mapped_on_TSS_analysis.ipynb` refers to
  `../T2T_article/T2T_repeat_masker_processed.csv`, which does not resolve from here. The
  working replacement is the local `T2T_repeat_masker_processed_sorted.bed` from step 1.1.
  The revision scripts all use the local file.

`GO_subfamilies.py` is deliberately left at `fdr_threshold=0.1` because it belongs to a
companion manuscript whose figures have not been regenerated. The G3 revision gets subfamily
GO at 0.05 from `revision_G3/06_go_rerun_fdr005.py` instead, which writes to
`revision_G3/output/` and leaves the shared `GO_tables/` alone.

---

## 3. The revision: long compute first

```bash
source ~/venvs/Retroelements_3_11/bin/activate
cd <repository root>

# Compact the legacy 11 GB permutation store to ~106 MB, verifying every seed before
# anything is deleted. Skip if you have the tracked permutation_counts_10kb/ store.
python revision_G3/01b_compact_permutation_results.py
python revision_G3/01c_expand_counts.py --seed 1 --check      # proves the compaction is lossless

# Build the 5 kb and 20 kb TSS windows and their backgrounds (step 1.4).
bash   revision_G3/05a_build_windows.sh
bash   revision_G3/01_permutations_stream.sh --window 5kb      # N = 500, 2-4 h
bash   revision_G3/01_permutations_stream.sh --window 20kb     # N = 500, 2-4 h
python revision_G3/01a_consolidate_counts.py --check-legacy    # regression against legacy 10 kb

# Re-run every Gene Ontology level at FDR 0.05 (reviewer minor comment 2).
python revision_G3/06_go_rerun_fdr005.py --level all           # 1-2 h
```

`01_permutations_stream.sh` streams per-seed counts instead of writing intermediate BEDs,
which is why the store is ~106 MB rather than 11 GB. To get a legacy-format BED back for any
single seed, use `01c_expand_counts.py --seed N`.

---

## 4. The revision: analysis

```bash
python revision_G3/02_ifna_domain_test.py          # 4 interferon-alpha tests, ~20 min
python revision_G3/04a_lu2020_geneset_overlap.py   # overlap with Lu et al. 2020 (needs 1.3)
python revision_G3/04b_newly_resolved_regions.py   # bound on the assembly contribution, ~10 min
python revision_G3/05b_window_sensitivity.py       # 6 classes + 44 families x 3 windows
python revision_G3/05c_percentile_sensitivity.py   # GO at top/bottom 10 % vs 5 %, ~10 min
# The GO grid: 3 TSS windows x 2 gene-set percentiles, all three GO levels.
# 07a's --verify-10kb is a GATE, not a formality: it refuses to write the 5 kb and 20 kb
# gene tables unless rebuilding the 10 kb one reproduces the published TEs_on_genes.csv
# exactly (rows, genes, empty windows, per-class totals, mean divergence to 1e-9, and the
# per-window multiset of family / class / subfamily names).
python revision_G3/07a_build_gene_tables.py --verify-10kb     # ~15 s
python revision_G3/07a_build_gene_tables.py --window 5kb      # ~10 s
python revision_G3/07a_build_gene_tables.py --window 20kb     # ~10 s
python revision_G3/07b_go_grid.py --window 5kb --window 20kb  # 4 new cells, ~45 min
python revision_G3/07b_go_grid.py --check-reuse               # the 10 kb cells ARE the published ones

python revision_G3/10_tables.py                    # Table1.csv + Table2.csv
python revision_G3/11_results_numbers.py           # every number the manuscript text quotes
```

`11_results_numbers.py` is the one to run if you only want to check the paper's numbers: it
re-derives all of them from the persisted outputs into
`revision_G3/output/results_numbers.{json,txt}`, so any disagreement between a figure, a
table and the text shows up there.

---

## 5. The revision: figures

**Order matters: `nb07` before `nb05`.** `nb05`'s panel S13D and its headline-claim table both
read `output/go_grid_headline_by_condition.csv`, which `nb07` writes, and `nb05` is the single
writer of `robustness_headline_claims_by_condition.csv`.

Every figure panel comes from a notebook, so each subfigure can be inspected before it is
used. **Nothing writes to Figma**; the notebooks emit SVG to `revision_G3/svg/` and
`revision_G3/svg/PLACEMENT.md` maps each file to its target frame for manual placement.

```bash
for nb in nb01_permutation_convergence nb02_ifna_domain nb03_relabelled_figures \
          nb06_go_networks_fdr005 nb07_go_grid_robustness \\
          nb05_sensitivity_robustness; do
  jupyter nbconvert --execute --to notebook --inplace \
      --ExecutePreprocessor.kernel_name=retroelements_3_11 \
      --ExecutePreprocessor.timeout=3600 "revision_G3/$nb.ipynb"
done
```

`nb06` asserts that no two labels overlap before it writes any SVG, and records the outcome
per panel in `revision_G3/output/network_qc.json`. If a panel cannot be drawn without
collisions it climbs a documented fallback ladder rather than writing an unreadable figure.
`network_qc.json` also records `canvas_area_vs_baseline`, `compaction_target_met` and the best
result any compaction rung reached. **Read Figure 6A's term count off that file**, never from
memory: the 1.2× font increase costs it four terms (5 per family, not 9), and the caption has to
quote what was achieved.

`nb03` asserts that Figures 2D–2F reproduce the three published significance sentences. If that
assertion fires, the panels have been put back at subfamily level or the D/E letters have been
swapped again — do not "fix" it by loosening the assertion.

---

## 5b. The revision: the supplementary package

```bash
python revision_G3/08_build_supplementary.py            # five workbooks, ~1 min
python revision_G3/08_build_supplementary.py --verify
```

Writes `revision_G3/supplementary/`: **Files S1–S5** as thematic Excel workbooks (8.7 MB total),
each opening with a README sheet, plus `README.md`, `INVENTORY.json`, `CHECKSUMS.sha256` and
`figures/PLACEHOLDERS.md`. `--verify` checks five workbooks, every sheet within Excel's limits,
every GO sheet at FDR < 0.05, and the checksums.

Two things to know. The script reads the **published** gene sets from
`manuscript_figures_supplementary_old/Supplementary File {3,5,7}.xlsx`, which is why those three
files alone are un-gitignored — without them a clean checkout cannot rebuild workbook S2. And the
April `Supplementary File 3.xlsx` shipped **two empty sheets** (`TE top`, `TE bottom`); 08
reconstructs them from the same construction `06_go_rerun_fdr005.py` uses, after checking that
reconstruction against the six non-empty sheets. Verify the checksums from inside the directory,
since they are bare filenames:

```bash
cd revision_G3/supplementary && sha256sum -c CHECKSUMS.sha256 && cd -
```

---

## 6. The revision: UCSC track hub

```bash
mkdir -p ~/bin && cd ~/bin
for t in bedToBigBed hubCheck bigBedInfo; do
  curl -sSfLO "https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/$t" && chmod +x "$t"
done
cd <repository root>
bash revision_G3/12_build_trackhub.sh          # ~5 min, writes revision_G3/trackhub/
```

The build converts the TE annotation, the TSS windows, the two gene sets and the
interferon-alpha domain to bigBed, writes `hub.txt` / `genomes.txt` / `hs1/trackDb.txt` and
eleven per-track description pages, then runs `hubCheck` and confirms every bigBed opens with
the expected record count.

Two expected messages:

- `hubCheck` reports `Can't get default spec from host genome.ucsc.edu` on any machine
  without a local `hg.conf`. The script tolerates that one message and fails on anything
  else, so a track-level problem cannot slip through.
- `12a_trackhub_beds.py` reports a small number of intervals clipped at a chromosome end.
  `T2T_genes.bed` extends 5 kb either side of each TSS without checking the chromosome
  boundary, so a few windows run past the end. That is harmless for `bedtools intersect` but
  `bedToBigBed` rejects it, so the hub clips those windows to the chromosome length and says
  how many.

The hub is **not** tracked on the main branch — it is published to `gh-pages`, so that a
clone does not carry 105 MB twice. Publishing and verifying are scripted:

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

GitHub Pages is required rather than Zenodo because UCSC needs HTTP range requests. The first
build takes a minute or two; watch it under the repository's Actions tab as
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

One-click load, opening at the interferon-alpha domain:

```
https://genome.ucsc.edu/cgi-bin/hgTracks?db=hs1&hubUrl=https://nikit357.github.io/T2T_transposons_genes/trackhub/hub.txt&position=chr9:21150692-21370055
```

---

## 7. Verifying the revision without re-running it

```bash
# The four frozen files are untouched.
md5sum -c <<'SUMS'
6d59a2a735b8d0f4fcf6d9dddbb8bb39  TEs_mapped_on_TSS_analysis.ipynb
a75ceaf51c0a0d221f53357bb0040b55  Gene_ontology_analysis.ipynb
3e8aec87bd9e78fce53463a2073d968b  download_and_process_files_UCSC_genes.ipynb
cfd78a7eb38b8f5bbc76dd0fba75dc01  GO_subfamilies.py
SUMS

# The permutation store reads back.
python -c "
import sys; sys.path.insert(0, 'revision_G3')
import revision_lib as rl
d = rl.load_permutation_counts(window='10kb', seeds=[1, 7, 500])
print(len(d), 'rows', sorted(d['seed'].unique()))"

# Every number the manuscript quotes, re-derived.
python revision_G3/11_results_numbers.py
```

---

## 8. The manuscript

Needs the `python-docx` environment, not the analysis one.

```bash
~/venvs/collagen_3_11/bin/python revision_G3/13_manuscript_tracked_edits.py
~/venvs/collagen_3_11/bin/python revision_G3/15_house_style.py          # must follow 13
~/venvs/collagen_3_11/bin/python revision_G3/14_build_extensive_discussion.py

# align the text with the figures as exported on 2026-08-07 — runs LAST, and unlike 13 and 15
# it writes a new file rather than editing in place
~/venvs/collagen_3_11/bin/python revision_G3/16_figure_alignment_edits.py --report   # dry run
~/venvs/collagen_3_11/bin/python revision_G3/16_figure_alignment_edits.py
```

**`16_figure_alignment_edits.py` reads `…260804_manual.docx` and writes
`…260807.docx`**, which is now the current state of the paper. Its input is never modified —
the verification block checks the 260804 md5 is unchanged — so a bad run is repaired by deleting
the output and re-running with `--force`. Because the output is a deterministic rebuild rather
than an accumulation of in-place edits, every edit must report `applied`; `already present` is a
warning and `not found` is fatal. It asserts 129 Mendeley content controls before and after.
Options: `--in` / `--out` / `--report` / `--force`. Report goes to
`revision_G3/output/figure_alignment_edit_report.json`.

**The working file is `T2T_genes_article_G3_revision_260804_manual.docx`**, which carries
Daniil's manual edits and his acceptance of 594 of the 1,124 tracked deletions. `13_…` edits it
**in place** (`15_house_style.py` targets the same file); it no longer rebuilds from the baseline,
because doing so would discard that work. `--from-baseline` still exists for a clean rebuild onto
a different filename and refuses to overwrite the working file.

Idempotency is now per-edit rather than by rebuilding: each edit reports `applied`,
`already present` (skipped) or **`not found`**, only the last is fatal, and every structural
insertion carries a marker sentence so a re-run cannot add a second copy of it. Outcome counts go
to `revision_G3/output/manuscript_edit_report.json`. Run `--dry-run` first to see what would
change without writing.

`14_build_extensive_discussion.py` still builds from the pristine baseline and is unaffected.

**Accept All gives the clean revision, but Reject All no longer restores the baseline** —
accepted deletions are permanent by design. For a full tracked diff against the submitted
version, use Word's **Compare** on the baseline and the final clean file, not the revision marks.

To read the result as text:

```bash
pandoc -t plain --wrap=none \
  revision_G3/Revised_manuscript/T2T_genes_article_G3_revision_260804_manual.docx -o /tmp/ms.txt
```

---

## 9. Files near GitHub's limits

| File | Size | Status |
|---|---|---|
| `T2T_repeat_masker_processed_sorted.bed` | 155 MB | gitignored, over the 100 MB hard limit — rebuild per step 1.1 |
| `goa_human.gaf` | 190 MB | gitignored — download per step 1.2 |
| `plots/length_vs_score_repeats_by_class.svg` | 659 MB | gitignored — regenerable by `draw_length_divergence_corr.py` |
| `revision_G3/trackhub/` | 105 MB | gitignored on this branch, published to `gh-pages` |
| `revision_G3/output/permutation_counts_{5,20}kb/` | 81 + 135 MB | gitignored — regenerable per step 3 |
| `revision_G3/output/TEs_on_genes_{5,10,20}kb.csv` | 9 + 15 + 27 MB | gitignored — rebuilt in ~10 s each by `07a`, whose gate proves the 10 kb one reproduces the published table |
| `revision_G3/output/GO_grid/*_fdr01_reference.csv` | 61 MB | gitignored — the FDR 0.05 files (55 MB) and all five derived tables are tracked; the 0.1 twins are the retrieval cut only, rebuilt by `07b` |
| `revision_G3/supplementary/File_S1_TE_TSS_map_and_enrichment.xlsx` | 7.0 MB | **tracked** — the largest workbook, comfortably under the limit; rebuilt by `08` |
| `TEs_on_genes_counts_subfamilies.csv` | 85 MB | **tracked, at 85 % of the hard limit** — see below |

`TEs_on_genes_counts_subfamilies.csv` compresses **75.5×**, from 85.1 MB to 1.1 MB, and
`TEs_on_genes_counts_subfamilies.csv.gz` is provided. Swapping the tracked raw file for the
gzip would take the repository from 85 % of GitHub's per-file limit to 1 %, but it changes a
file that the Data availability statement names, so it is left as an explicit decision:

```bash
git rm --cached TEs_on_genes_counts_subfamilies.csv
echo TEs_on_genes_counts_subfamilies.csv >> .gitignore
git add TEs_on_genes_counts_subfamilies.csv.gz
# then change the Data availability statement to name the .gz
```

---

## 10. Known discrepancies, stated deliberately

1. **N = 1000 vs N = 500.** The frozen download notebook contains `NUM_PERMUTATIONS = 1000`
   and was never run at that value. N = 500 is authoritative: 501 files in the original
   store, empirical p floor 2/501 = 0.004, and the corrected Methods say 500. The notebook is
   frozen, so the discrepancy is documented here, in `CLAUDE.md` and in `revision_G3/README.md`
   rather than edited away.
2. **GO FDR 0.1 vs 0.05.** `GO_subfamilies.py` stays at 0.1 for the companion subfamilies
   manuscript. The G3 revision uses 0.05 everywhere via `revision_G3/06_go_rerun_fdr005.py`.
3. **`GO_tables/` is the FDR 0.1 output.** The 0.05 equivalents are in
   `revision_G3/output/GO_tables_fdr005/`, and the pre-revision snapshot is preserved in
   `revision_G3/output/legacy_fdr01_n500/`.
4. **Element counts per gene are per-TSS.** A gene with several annotated TSS contributes
   several windows, so an element within 10 kb of two TSS of the same gene is counted twice.
   This is a property of the published design, is flagged in the manuscript's Limitations,
   and is reproduced deliberately by the revision scripts and the track hub.
