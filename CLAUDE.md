# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Genomics research analyzing transposable elements (TEs) and their proximity to genes in the
T2T-CHM13 (Telomere-to-Telomere CHM13) human genome assembly. The analysis maps 3,709,429 TEs to
28,738 genes using a 10 kb TSS window, computes enrichment statistics for 6 TE classes, 44
families and 1,143 subfamilies, and performs Gene Ontology (GO) enrichment to interpret TE–gene
co-localization and evolutionary arms races (notably young L1 clusters in the IFNA domain).

Two manuscripts come out of this directory:
- The class/family-level paper (preprint DOI [10.32942/X2FM2M](https://doi.org/10.32942/X2FM2M)),
  **conditionally accepted at *G3: Genes|Genomes|Genetics* as G3-2026-406828** — supplementary
  deliverables are the root-level `Supplementary Figure *.pdf` / `Supplementary File *`.
- The subfamily-level follow-up, currently being prepared for PLOS ONE. Its figures and drafted
  Results text live in `T2T_genes_subfamilies_article_figures/` (which has its own `CLAUDE.md`
  describing every figure panel).

## The G3 revision — read this before touching anything

The revision lives in **`revision_G3/`**, which is authoritative for everything it changes, and
is planned in `G3_revision_implementation_plan_260803.md`. `REPRODUCE.md` at the repository root
is the run order from a clean checkout and is the first thing to read if you want to re-run rather
than read.

**Four files are frozen for the revision and must not be edited:**

```
6d59a2a735b8d0f4fcf6d9dddbb8bb39  TEs_mapped_on_TSS_analysis.ipynb
a75ceaf51c0a0d221f53357bb0040b55  Gene_ontology_analysis.ipynb
3e8aec87bd9e78fce53463a2073d968b  download_and_process_files_UCSC_genes.ipynb
cfd78a7eb38b8f5bbc76dd0fba75dc01  GO_subfamilies.py
```

Their helper functions were **copied** into `revision_G3/revision_lib.py` and modified there
(caveat C19), so the same function exists twice on purpose. What came from where:

| Frozen source | Copied into `revision_lib.py` as |
|---|---|
| `Gene_ontology_analysis.ipynb` cell 6 | `run_goatools_enrichment`, `run_goatools_ordered_enrichment` (now `fdr_threshold=0.05`), `save_go_network_svg` |
| GO nb cells 36, 175 | `visualize_go_class_network`, `save_go_network_svg_families_by_classes` (plus `min_shared_genes`, `max_term_genes`, the collision check) |
| `TEs_mapped_on_TSS_analysis.ipynb` cell 40 | `load_permutation_counts(window=...)`, the compacted-store reader |
| GO nb cell 3 | the TE class palette and font settings |

Four discrepancies are documented rather than edited away, because editing them would break the
freeze or another paper's figures:

1. **N = 500 is authoritative, not 1,000.** Cell 34 of the download notebook sets
   `NUM_PERMUTATIONS = 1000` and **was never executed at that value**: 501 files exist in the
   original `permutation_results/`, the empirical p floor is 2/501 = 0.004, and the manuscript
   Methods were corrected *down* to 500 — not the other way round (caveat C20). The line stays in
   the frozen notebook.
2. **`GO_subfamilies.py` stays at `fdr_threshold=0.1`** because it belongs to the companion
   subfamilies manuscript, whose figures have not been regenerated (caveat C3). The G3 revision
   gets subfamily GO at 0.05 from `revision_G3/06_go_rerun_fdr005.py`, which writes to
   `revision_G3/output/` and leaves the shared `GO_tables/` untouched.
3. **The frozen `TEs_mapped_on_TSS_analysis.ipynb` still contains the broken
   `../T2T_article/…` path.** The row in "Working Directory Gotchas" below is therefore
   *superseded*, not *fixed*: the working replacement is the local
   `T2T_repeat_masker_processed_sorted.bed`.
4. **Element counts are per-TSS, not per-gene.** A gene with several annotated TSS contributes
   several windows, so an element within 10 kb of two TSS of the same gene is counted twice. This
   is a property of the published design; it is flagged in the manuscript's Limitations and
   reproduced deliberately by the revision scripts and the track hub.

### The compacted permutation store

`permutation_results/` (11 GB, 501 legacy BEDs) was replaced by a per-seed count store,
**~100× smaller**, and removed only after every seed was verified:

| | |
|---|---|
| Location | `revision_G3/output/permutation_counts_10kb/` (106 MB, tracked), `_5kb/` (81 MB) and `_20kb/` (135 MB) (both gitignored, regenerable) |
| Schema | `seed, score, subfamily_name, family_name, class_name, n`, one zstd-19 file per seed |
| Metadata | `MANIFEST.json` per store: script, version, creation time, source directory, window, seed list |
| Read it | `revision_lib.load_permutation_counts(window="10kb", seeds=[...])` |
| Reconstruct a legacy BED | `revision_G3/01c_expand_counts.py --seed N` |

### The manuscript files

All three live in `revision_G3/Revised_manuscript/`, moved out of
`T2T_genes_subfamilies_article_figures/` where they were originally filed with the other paper:
the read-only submitted baseline, the tracked revision, and the extended-discussion supplementary
file. `revision_G3/13_manuscript_tracked_edits.py` and `14_build_extensive_discussion.py` both
read the baseline and write beside it, so they are idempotent.

Citations are **Mendeley Cite** content controls (128 `<w:sdt>` tagged `MENDELEY_CITATION` plus one
`MENDELEY_BIBLIOGRAPHY`), with the payload in `word/webextensions/webextension1.xml` and **no**
legacy field codes. Consequences: any tracked edit must operate on whole `<w:sdt>` elements —
never rewrite a paragraph from its concatenated text, which is what the `word_rewrite` skill's
`tracked_replace` does and why `13_…` implements its own citation-safe replacement. The
tracked-changes helpers need `python-docx` + `lxml`, which are in **`~/venvs/collagen_3_11`**, not
in `~/venvs/Retroelements_3_11` (caveat C15).

### Figures are made in Figma, not by Python

The plotting scripts and notebooks produce **panels only, as SVG**. Final figures are composed,
lettered and exported by hand in Figma: file key `WRNeTzKZObdmAQ8QG1EZlq`, single page `0:1`, one
frame per figure. `revision_G3/svg/PLACEMENT.md` maps every SVG to its target frame ID, and
`G3_figure_pvalue_labels_260803.md` gives the correct raw-vs-FDR label for every statistic in every
panel. **Nothing in `revision_G3/` writes to Figma.** Note the `Figure 9` name collision: the same
canvas carries the subfamilies paper's labels, so scope any Figma work **by node ID, not by name**
(caveat C9).

### The UCSC track hub

`revision_G3/12_build_trackhub.sh` builds a hs1 track hub from the canonical inputs: one bigBed per
TE class coloured with the project palette, the 38,704 TSS windows, the TE-top and TE-bottom gene
sets, and the interferon-alpha domain. It is gitignored on this branch and published to `gh-pages` by
`revision_G3/12b_publish_trackhub.sh`, then checked live by `12c_verify_trackhub_live.sh`. GitHub
Pages must be enabled once by hand on the `gh-pages` branch; until that is done the URL the
manuscript prints returns 404. Publishing on a separate branch is what keeps UCSC's HTTP range
requests satisfied without adding 105 MB to every clone.

## Working Directory Gotchas (read this before running anything)

Relative paths inside the notebooks and scripts are **inconsistent** — they were written across
several sessions with different working directories. Verify a path exists before trusting it.

| Referenced path | Where the data actually is |
|---|---|
| `./epigenomic_files/...` | `../epigenomic_files/` — one level up, in `T2T_genes_article/` (54 mapped BED/bedGraph files) |
| `./epigenomic_files/permutation_results/` | `../epigenomic_files/permutation_results/` (501 files) |
| `../T2T_article/T2T_repeat_masker_processed.csv` | **RESOLVED** — superseded by the local `T2T_repeat_masker_processed_sorted.bed` (below). The broken path still stands in `draw_length_divergence_corr.py:30` and in one cell of the frozen `TEs_mapped_on_TSS_analysis.ipynb`; every revision script uses the local file instead |
| `chm13.genome`, `T2T_genes_sorted.bed` | `chm13.genome` and `T2T_genes.bed` **are** present locally now |

### Canonical inputs — use these and only these for TE ↔ TSS intersections

| File | Rows | Columns | Why it is safe |
|---|---|---|---|
| `T2T_repeat_masker_processed_sorted.bed` | 3,709,429 | `chrom, start, end, score(divergence), subfamily, family, class` | Byte-size-identical to `../epigenomic_files/repeats_all.bed`, the file the legacy permutations shuffled. **155 MB, gitignored** — over GitHub's per-file limit; rebuild per `REPRODUCE.md` §1.1 |
| `T2T_genes.bed` | 38,704 | `chrom, start, end, gene` | Interval set exactly identical to the `mapped_on_TSS` bedGraph the legacy pipeline used as its `-b` file |
| `chm13.genome` | 24 | `chrom, size` | chr1–22, X, Y; used by `bedtools shuffle` and `bedToBigBed` |

Because the first two are provably interchangeable with what the legacy pipeline consumed, the
retained 10 kb N = 500 background and the new 5 kb / 20 kb runs are directly comparable — there is
no hidden geometry change to explain in the Methods. Note that `T2T_genes.bed` extends 5 kb either
side of each TSS **without clipping at chromosome ends**, so a few windows run past the end of their
chromosome; harmless for `bedtools intersect`, but `bedToBigBed` rejects it and the track hub build
clips and reports them.

The master RepeatMasker table (`T2T_repeat_masker_processed.csv`) is **shared with the sibling
project** `../../T2T_article/`. Do not duplicate it — that project owns it.

## Running the Analysis

### Notebooks (run in order)
```bash
jupyter notebook download_and_process_files_UCSC_genes.ipynb   # Data acquisition + generates the permutation script
jupyter notebook TEs_mapped_on_TSS_analysis.ipynb              # Core proximity mapping, enrichment, Figures 1-2 + Suppl. 1,2,3,5,6
jupyter notebook Gene_ontology_analysis.ipynb                  # GO enrichment, Figures 3,4,5,8
```

### Python scripts
```bash
python GO_subfamilies.py                    # Per-subfamily GO enrichment -> GO_tables/ (230 files)
python genes_subfamilies_network.py         # "Ring of Power" network, Figure 6 (JACCARD_THRESHOLD = 0.025)
python genes_subfamilies_network_clusters.py # MCL-clustered L1 subnetwork, Figure 7 (JACCARD_THRESHOLD = 0.01, MCL_INFLATION = 1.14)
python draw_length_divergence_corr.py       # TE length vs divergence correlation (fix the ../T2T_article path first)
```

### The G3 revision (`revision_G3/`) — compute, figures, manuscript

Compute lives in `.py`/`.sh` so it is cacheable and backgroundable; every figure lives in a
notebook so each subfigure can be inspected before use. Full run order with timings is in
`REPRODUCE.md`; `revision_G3/README.md` has the directory layout.

```bash
# compute
python revision_G3/01b_compact_permutation_results.py   # 11 GB -> ~106 MB, verified per seed
bash   revision_G3/01_permutations_stream.sh --window 5kb   # and --window 20kb, N=500, 2-4 h each
python revision_G3/06_go_rerun_fdr005.py --level all    # every GO level at FDR 0.05
python revision_G3/02_ifna_domain_test.py               # 4 interferon-alpha domain tests
python revision_G3/04a_lu2020_geneset_overlap.py        # overlap with Lu et al. 2020
python revision_G3/04b_newly_resolved_regions.py        # bound on the assembly contribution
python revision_G3/05b_window_sensitivity.py            # enrichment at 5/10/20 kb
python revision_G3/05c_percentile_sensitivity.py        # GO at top/bottom 10 % vs 5 %
python revision_G3/10_tables.py                         # Table1.csv + Table2.csv
python revision_G3/11_results_numbers.py                # every number the manuscript quotes

# figures (SVG panels only; placement is manual, see svg/PLACEMENT.md)
jupyter nbconvert --execute --inplace revision_G3/nb0{1,2,3,5,6}_*.ipynb

# manuscript (needs python-docx: the collagen venv, NOT Retroelements_3_11)
~/venvs/collagen_3_11/bin/python revision_G3/13_manuscript_tracked_edits.py
~/venvs/collagen_3_11/bin/python revision_G3/14_build_extensive_discussion.py

# deliverable
bash revision_G3/12_build_trackhub.sh                   # bigBeds + hub + hubCheck, ~5 min
```

### Permutation test (statistical background)
`run_permutation_test.sh` is **not checked in** — it is generated by a cell in
`download_and_process_files_UCSC_genes.ipynb` (search for `bash_script_path`). That cell writes a
script which, per seed, runs `bedtools shuffle -i repeats_all.bed -g chm13.genome -seed $SEED`
followed by `bedtools intersect -b <TSS bedGraph> -wa`, fanned out over `os.cpu_count()` cores via
`xargs -P`. Inputs are resolved relative to `epigenomic_files/`.

```bash
# after generating it from the notebook, from the directory holding epigenomic_files/
bash run_permutation_test.sh
```

**Permutation count is inconsistent across the codebase.** The generator sets
`NUM_PERMUTATIONS = 1000`; the consuming cell in `TEs_mapped_on_TSS_analysis.ipynb` sets
`NUM_PERMUTATIONS = 500`; `permutation_results/` holds 501 files; the manuscript and `results.md`
report 1,000 permutations. Resolve this with Daniil before quoting a number in any figure or text.

## Architecture & Data Pipeline

```
download_and_process_files_UCSC_genes.ipynb
  -> T2T RefSeq annotations, TSS coordinates, chm13.genome, repeats_all.bed
  -> generates run_permutation_test.sh (1,000x bedtools shuffle + intersect)

TEs_mapped_on_TSS_analysis.ipynb
  -> reads the 10 kb TSS neighborhood BEDs from ../epigenomic_files/*.mapped_on_TSS_with_TEs.bed
  -> consolidates permutation_results/ into consolidated_random_data.csv
  -> Fisher exact test (observed) vs permutation background (random OR, empirical p)
  -> TEs_on_genes.csv (23 MB), TEs_on_genes_counts_subfamilies.csv (86 MB)
  -> enrichment_families_with_random.csv, enrichment_subfamilies_with_random*.csv

Gene_ontology_analysis.ipynb + GO_subfamilies.py
  -> goatools GOEnrichmentStudy against go-basic.obo + goa_human.gaf (190 MB, not in git)
  -> foreground = top 5% genes by TE count per subfamily/family; FDR threshold 0.1
  -> GO_tables/*.csv (230 files, one per subfamily) + top_5_perc_genes_by_*.csv

genes_subfamilies_network*.py
  -> Jaccard similarity between subfamily gene sets -> NetworkX graph
  -> custom force-directed layout + collision relaxation (all parameters at the top of each script)
  -> plots/*.svg
```

The permutation background is what makes the enrichment numbers meaningful: raw Fisher odds ratios
scale with element length, so SINEs (~300 bp) would be systematically under-called and
LINEs/LTRs (~6 kb) over-called without it. Always report the observed-to-random OR ratio, not the
raw OR.

**The correlation was re-measured on 2026-08-04 and this file previously recorded it wrongly
as R = 0.661.** Mean element length against the mean random OR gives **Pearson R = 0.985**
across the 44 families (n = 44, p = 8.4 × 10⁻³⁴, mean lengths 122–6,357 bp) and R = 0.983 across
the 1,130 subfamilies. Concretely: Alu averages 316 bp and a random OR of 1.535, L1 averages
6,357 bp and a random OR of 2.664. The figures are persisted to
`revision_G3/output/length_bias_correlation.json` and quoted in the revised Methods; use them
rather than 0.661.

## Key Files

| File | Purpose |
|------|---------|
| `TEs_on_genes.csv` | Master TE–gene intersection table (10 kb TSS window) |
| `TEs_on_genes_counts_subfamilies.csv` | Per-gene TE counts by subfamily; the input to every subfamily script |
| `individuals_by_classes_TE.csv` / `families_by_classes_TE.csv` | Subfamily -> class and family -> class maps; required by the network and GO scripts |
| `enrichment_families_with_random.csv` | Enrichment stats for 44 families |
| `enrichment_subfamilies_with_random.csv` | Per-subfamily enrichment stats (`_length`, `_and_gene_sets_stats` variants add element length and gene-set size columns) |
| `top_5_perc_genes_by_subfamilies.csv` / `_families.csv` | Foreground gene sets used for GO enrichment |
| `GO_Classified_Table_Ordered_Gemini.tsv` | Manually curated GO -> metagroup classification |
| `Subfamilies_classified - Classification.csv` | Manual subfamily -> functional-group annotation used in figures |
| `sequence_report.tsv` | T2T chromosome lengths; converted to `chm13.genome` for `bedtools shuffle` |
| `go-basic.obo` | Gene Ontology DAG (31 MB, tracked) |
| `goa_human.gaf` | Human GO annotations (190 MB, gitignored — download separately) |
| `T2T_genes_subfamilies_article_figures/results.md` | Drafted Results section for the subfamilies paper; the authoritative narrative for what each figure shows |

## Plotting Conventions

Every plotting script sets these — match them in new code:

```python
plt.rcParams["svg.fonttype"] = "none"   # keeps text editable in Illustrator/Inkscape
GLOBAL_FONT_SIZE = 10
```

Shared TE class palette (keep identical across all figures in both manuscripts):

```python
{"LINE": "#cc660b", "LTR": "#70453c", "SINE": "#ab1f20",
 "DNA": "#195f90", "Retroposon": "#765297", "RC": "#238023"}
```

## Environment

Python 3.11. Key packages: `pandas`, `numpy`, `scipy`, `statsmodels`, `goatools`, `matplotlib`,
`seaborn`, `networkx`, `pyvis`, `supervenn`, `plotly`, `adjustText`, `statannotations`,
`scikit-learn` (silhouette scoring in the MCL script), `tqdm`. External tool: `bedtools`.

`goa_human.gaf` (190 MB) must be downloaded separately — it is the only real gitignore entry.
Note that `plots/` **is** tracked (270 files), including several ~30 MB SVGs; the repo is already
heavy, so prefer writing new large intermediates outside it or adding them to `.gitignore`.

## Key Parameters to Tune

- `JACCARD_THRESHOLD` — `0.025` in `genes_subfamilies_network.py`, `0.01` in
  `genes_subfamilies_network_clusters.py`; controls edge density
- `MCL_INFLATION` (`1.14`) in `genes_subfamilies_network_clusters.py` — cluster granularity;
  `MCL_BENCHMARK_RANGE` sweeps it for the silhouette benchmark plot
- FDR threshold in `run_goatools_enrichment()` in `GO_subfamilies.py` — stays at `0.1` for the
  companion paper. **The G3 revision uses 0.05 everywhere, with no "suggestive" band** (decision
  D2), via `revision_lib.run_goatools_enrichment(fdr_threshold=0.05)`
- TSS window size — 10 kb for every published number. The revision additionally reports
  **5 kb and 20 kb** (`revision_G3/05a_build_windows.sh`, `05b_window_sensitivity.py`)
- Top-N percentile cutoff for foreground gene sets — top/bottom **5 %** (`CLASS_TOP_N = 1436`) for
  every published number, with **10 %** (2,872 genes) as the sensitivity arm
  (`revision_G3/05c_percentile_sensitivity.py`)
- Network parameters in `revision_lib.save_go_network_svg*`: main-text panels use
  `top_n=10, jaccard_threshold=0.2, min_shared_genes=5, max_term_genes=500`; supplementary full
  versions use `top_n=30, jaccard_threshold=0.1, min_shared_genes=0, max_term_genes=1000`.
  `assert_no_label_collisions()` refuses to write an SVG with overlapping labels, so a colliding
  figure cannot be produced silently — outcomes per panel are in `revision_G3/output/network_qc.json`
- `SUBFAMILY_LIMIT` — `1143`, i.e. all annotated subfamilies (1,129 of which appear in at least one
  TSS window)
- `N_PERMUTATIONS` — **500** (not 1,000; see the freeze notes above). Empirical p floor 2/501 = 0.004

## Numbers worth not re-deriving

| Quantity | Value | Where it comes from |
|---|---|---|
| Length bias the permutation control removes | mean element length vs random OR, **Pearson R = 0.985** across 44 families (n = 44, p = 8.4 × 10⁻³⁴, lengths 122–6,357 bp); Alu 316 bp → OR 1.535, L1 6,357 bp → OR 2.664 | `revision_G3/output/length_bias_correlation.json`. **An earlier version of this file said 0.661 — that was wrong** |
| GO terms at FDR 0.1 → 0.05 | classes 504 → 425; divergence 516 → 414; families 195 → 140; subfamilies 1,231 → 1,003 | `revision_G3/output/results_numbers.json` |
| Families with ≥ 1 GO term at 0.05 | **13 of 44** (Dong-R4 is the one lost) | same |
| Interferon-alpha domain | 175 TEs, 77 L1 from 36 subfamilies, mean L1 divergence 135.7 vs 188.2 genome-wide, 1.94× L1 density | `revision_G3/output/ifna_qc.json` |
| Newly resolved sequence (bound on the assembly contribution) | 208.8 Mb (6.70 %), but only 0.41 % of TEs, 0.49 % of windows, 0.55 % of genes — and **0 bp in the IFNA domain** | `revision_G3/output/assembly_bound_summary.json` |

## Housekeeping

`TEs_mapped_on_TSS_analysis-Copy1.ipynb` is a 29 MB stale duplicate of the main notebook. Treat
`TEs_mapped_on_TSS_analysis.ipynb` as authoritative and confirm with Daniil before deleting the copy.
