# CLAUDE.md — `revision_G3/`

Guidance for Claude Code working in this directory — read it before running or editing anything here.
Background, per-work-package detail, every measured result and all 20 caveats are in
**`project_overview.md`** beside this file; this file is the operating manual.

## What this directory is

The complete revision package for manuscript **G3-2026-406828** — *"Telomere-to-telomere co-mapping
of transposable elements and human genes identifies a cluster of young L1 elements in the
interferon-alpha domain"* — conditionally accepted at *G3: Genes|Genomes|Genetics* (letter
2026-07-29, 30-day window, target submission on or before **2026-08-28**).

`revision_G3/` is **authoritative** for everything the revision changes. Where a helper exists both
here and in a parent-directory notebook, the copy here is what the revised manuscript reports.

| Input document | What it is |
|---|---|
| `../G3_revision_implementation_plan_260803.md` | the plan: WP1–WP16, decisions D1–D9, caveats C1–C20, and a per-item completion log |
| `../G3_reviewer_report_260802.md` | the decision letter (6 major + 6 minor + editor items) |
| `../G3_article_guidelines.md` | G3 house style + gap list G1–G18 + 10 items to verify in a browser |
| `../G3_figure_pvalue_labels_260803.md` | correct raw-vs-FDR label for every statistic in every panel |
| `../REPRODUCE.md` | run order from a clean checkout; read this before re-running anything |
| `README.md` · `svg/PLACEMENT.md` | the human-facing version of this file · every SVG → target Figma frame ID and panel position |
| `external/PROVENANCE.md` | third-party inputs: URLs, sizes, md5s, and the download routes that fail |

## Hard rules

1. **Four files in the parent directory are frozen — zero edits, not a cell, not a label, not a
   threshold.** Verify with `md5sum -c` (below) before submission.

   ```
   6d59a2a735b8d0f4fcf6d9dddbb8bb39  ../TEs_mapped_on_TSS_analysis.ipynb
   a75ceaf51c0a0d221f53357bb0040b55  ../Gene_ontology_analysis.ipynb
   3e8aec87bd9e78fce53463a2073d968b  ../download_and_process_files_UCSC_genes.ipynb
   cfd78a7eb38b8f5bbc76dd0fba75dc01  ../GO_subfamilies.py
   ```

2. **Never import from a frozen notebook** — no `nbimporter`, no `import_ipynb`. That reintroduces
   the coupling the freeze removes and silently picks up FDR 0.1. Helpers were **copied** into
   `revision_lib.py` and modified there; the only legitimate cross-script import is
   `importlib.import_module("06_go_rerun_fdr005")` in `04a`/`05c`, which keeps their gene-set
   constructions identical to the published ones.
3. **Match every panel to its published caption, and assert it.** Three panel-identity
   defects have now been found by the same method — Figure 7, Supplementary 8B, and Figures
   2D–2F — each one a panel matched to a letter by *content* rather than against the caption
   the journal already has. Figure 2D–2F was at the wrong analysis level (subfamily, not
   family) *and* had D/E swapped. When a panel's identity is settled, encode the check as an
   assertion in the notebook: `nb03` cell 17 asserts the three published significance
   statements, so putting the panels back at the wrong level cannot pass silently.
4. **Nothing here writes to Figma.** Notebooks emit SVG panels only; placement, lettering and export
   are Daniil's manual step. Do not add `use_figma` calls.
5. **`13_…` now edits the working file IN PLACE**, and `15_house_style.py` targets the same
   file. The old rule ("15 must follow 13, every time, because 13 rebuilds from the pristine
   baseline") is **superseded**: 13 no longer copies the baseline over the target, so it no
   longer destroys the Phase 7 pass — and `15_…`'s refusal to run twice is now correct
   behaviour rather than an obstacle. See "The manuscript files" below.
6. **Font sizes go through `rl.fs()`, never as bare literals.** Matplotlib sizes are points and
   SVG consumers render at 96 dpi, so 1 pt = 96/72 px: `FONT_SCALE = 1.2` makes the base 12 pt
   land at exactly **16 px** in Figma. A bare `fontsize=8` would sit at 10.67 px against 16 px
   titles, which is the defect the scaling exists to remove. 56 literals across the five
   notebooks were converted; zero remain.
7. **Never delete a source file whose verification did not pass** (the compaction rule; `MANIFEST.json`
   plus `01c_expand_counts.py` are what make the 11 GB deletion defensible).
8. **GO FDR is 0.05 everywhere, with no "suggestive" band** (decision D2) — main text, supplementary
   figures and supplementary tables alike.
9. **`goa_human.gaf` and `../go-basic.obo` must not be re-downloaded.** A silent refresh changes term
   memberships and invalidates every measured FDR shift.

## Environment

| Item | Value |
|---|---|
| Analysis / plotting venv | `~/venvs/Retroelements_3_11` — `goatools`, `adjustText` **1.3.0**, `networkx`, `pyvis`, `supervenn`, `statannotations` |
| **docx / tracked-changes venv** | **`~/venvs/collagen_3_11`** — the only venv with `python-docx` 1.2.0 + `lxml`, which `~/.claude/skills/word_rewrite_trackchanges.py` needs. The two venvs stay independent on purpose |
| `bedtools` | 2.31.1 at `/usr/local/bin/bedtools` |
| `zstd` · cores | 1.5.7 at `/opt/conda/bin/zstd` · 8, so `xargs -P 8` |

`adjustText` 1.3.0 renamed the layout kwargs: `force_points` → `force_static`, `expand_text` →
`expand` — use the 1.3.0 names. pandas 3 also broke two patterns in the copied frozen code,
`float(Series)` and `np.fill_diagonal(df.values, …)` (a DataFrame's `.values` is read-only under
copy-on-write); both are already fixed in the notebooks.

## Canonical inputs — use these and only these for TE ↔ TSS intersections

| File | Rows | Columns | Why it is safe |
|---|---|---|---|
| `../T2T_repeat_masker_processed_sorted.bed` | 3,709,429 | `chrom, start, end, score(divergence), subfamily, family, class` | **md5-identical** to `../../epigenomic_files/repeats_all.bed`, which the legacy permutations shuffled. 155 MB, **gitignored** (over GitHub's limit) — rebuild per `../REPRODUCE.md` §1.1 |
| `../T2T_genes.bed` | 38,704 | `chrom, start, end, gene` | interval set exactly identical to the `mapped_on_TSS` bedGraph used as the legacy `-b` file |
| `../chm13.genome` | 24 | `chrom, size` | `bedtools shuffle` / `bedToBigBed` |

They supersede the broken `../../T2T_article/T2T_repeat_masker_processed.csv` path. Because both are
provably interchangeable with what the legacy pipeline consumed, the retained 10 kb N = 500 background
and the new 5 kb / 20 kb runs are directly comparable. `T2T_genes.bed` extends 5 kb either side of
each TSS **without clipping at chromosome ends** — harmless for `bedtools intersect`, but
`bedToBigBed` rejects it, so the track hub build clips and reports the count.

`external/` holds the only three files this repository does not produce (Lu et al.'s supplementary
tables, MGI mouse–human homology, UCSC hs1→hg38 chains). Their categories are **mouse** gene sets,
which is why WP4 has to map through orthology.

## Run order

Compute lives in `.py`/`.sh` (cacheable, backgroundable); every figure lives in a notebook so each
subfigure can be inspected before use. Full version with timings: `../REPRODUCE.md`.

```bash
source ~/venvs/Retroelements_3_11/bin/activate
cd /home/jovyan/Projects/Retroelements/T2T_genes_article/T2T_transposons_genes

# --- Phase 1: free the disk FIRST, then launch the long jobs ------------------
python revision_G3/01b_compact_permutation_results.py      # 11 GB -> ~106 MB, ~12 min
python revision_G3/01c_expand_counts.py --seed 1            # lossless proof
bash   revision_G3/05a_build_windows.sh                     # 5 kb / 20 kb windows
bash   revision_G3/01_permutations_stream.sh --window 5kb    # N=500, ~16 min
bash   revision_G3/01_permutations_stream.sh --window 20kb   # N=500, ~18 min
python revision_G3/01a_consolidate_counts.py --window 5kb    # and 20kb
python revision_G3/06_go_rerun_fdr005.py --level all         # FDR 0.05, ~20 min w/ the DAG cache
python revision_G3/07a_build_gene_tables.py --verify-10kb    # gate: must reproduce the published table
python revision_G3/07a_build_gene_tables.py --window 5kb     # and --window 20kb, ~10 s each
python revision_G3/07b_go_grid.py --window 5kb --window 20kb  # 4 new GO cells, ~45 min, backgroundable
python revision_G3/07b_go_grid.py --check-reuse              # the 10 kb cells ARE the published ones

# --- Phase 2: analysis -------------------------------------------------------
python revision_G3/02_ifna_domain_test.py                   # 4 IFNA tests
python revision_G3/04a_lu2020_geneset_overlap.py            # Lu et al. overlap
python revision_G3/04b_newly_resolved_regions.py            # assembly bound (C6)
python revision_G3/05b_window_sensitivity.py                # 6 classes + 44 families x 3 windows
python revision_G3/05c_percentile_sensitivity.py            # GO at 10 % vs 5 %
python revision_G3/10_tables.py                             # Table1.csv + Table2.csv

# --- Phase 4: figures (notebooks; SVG only, no Figma) ------------------------
# nb07 must run BEFORE nb05: nb05's S13D and its claim table both read nb07's
# go_grid_headline_by_condition.csv, and nb05 is the single writer of the claim table.
jupyter nbconvert --execute --inplace revision_G3/nb0{1,2,3,6}_*.ipynb
jupyter nbconvert --execute --inplace revision_G3/nb07_go_grid_robustness.ipynb
jupyter nbconvert --execute --inplace revision_G3/nb05_sensitivity_robustness.ipynb

# --- Phase E: the supplementary package --------------------------------------
python revision_G3/08_build_supplementary.py                # five workbooks, ~1 min
python revision_G3/08_build_supplementary.py --verify

# --- Phase 5/7: manuscript (collagen venv — NOT Retroelements) ---------------
python revision_G3/11_results_numbers.py                    # every number the text quotes
~/venvs/collagen_3_11/bin/python revision_G3/13_manuscript_tracked_edits.py
~/venvs/collagen_3_11/bin/python revision_G3/15_house_style.py    # MUST follow 13
~/venvs/collagen_3_11/bin/python revision_G3/14_build_extensive_discussion.py

# --- Phase 6: deliverable ----------------------------------------------------
bash revision_G3/12_build_trackhub.sh                       # bigBeds + hub + hubCheck
```

Most scripts take `--summary` or `--report` to print their result without recomputing.
`13_…` **edits the working file in place** and is idempotent by per-edit detection rather than
by rebuilding: every edit reports `applied` / `already present` / **`not found`**, and only the
last is fatal, so a silent skip is still impossible while a re-run is no longer destructive. It
also tolerates the one rename `15_…` performs (`Supplementary File 4` → `File S4`), so its own
search strings still match a document 15 has already touched. `15_…` edits the same file,
refuses to start if the Phase 5 edits are absent, and refuses to run twice.

## Directory layout

```
revision_lib.py                     shared helpers; AUTHORITATIVE for the revision
01_permutations_stream.sh           N=500 streaming permutations, --window 5kb|10kb|20kb
01a_consolidate_counts.py           per-permutation summaries + pipeline-equivalence check
01b_compact_permutation_results.py  11 GB -> 106 MB, verified per seed before any delete
01c_expand_counts.py                reconstruct a legacy BED from counts (--seed N)
02_ifna_domain_test.py              four interferon-alpha domain tests
04a_lu2020_geneset_overlap.py       gene-set overlap with Lu et al. 2020 (via mouse orthology)
04b_newly_resolved_regions.py       ceiling on the assembly-attributable component (C6)
05a_build_windows.sh                5 kb / 20 kb TSS neighbourhoods + TE mapping
05b_window_sensitivity.py           enrichment at 5/10/20 kb
05c_percentile_sensitivity.py       GO at top/bottom 10 % vs 5 %
06_go_rerun_fdr005.py               every GO level at FDR 0.05  (imported by 04a, 05c, 07b, 12a)
07a_build_gene_tables.py            per-window TE-TSS tables; --verify-10kb is the regression gate
07b_go_grid.py                      GO across 3 windows x 2 percentiles; reuses the two 10 kb cells
08_build_supplementary.py           the five thematic Excel workbooks, Files S1-S5
10_tables.py                        Table1.csv + Table2.csv (reformat; refuses if a value moves)
11_results_numbers.py               re-derives every number the manuscript quotes
12_build_trackhub.sh + 12a_…_beds.py  hs1 track hub: 10 bigBeds + trackDb + 11 description pages
13_manuscript_tracked_edits.py      Phase 5 edits D–K, tracked, citation-safe, idempotent
15_house_style.py                   Phase 7 G3 house style G2–G16; MUST run after 13
14_build_extensive_discussion.py    Extensive_discussion docx by copy-and-delete
nb01…nb07 (6 notebooks)             the only figure surface; 44 SVGs -> svg/
supplementary/                      the deliverable: 5 workbooks + README + CHECKSUMS, 8.7 MB
Revised_manuscript/                 baseline (read-only) + revision + extensive discussion
external/                           third-party inputs; see PROVENANCE.md
output/                             every derived table, log and QC JSON
  permutation_counts_10kb/          compacted legacy store, 106 MB, TRACKED
  permutation_counts_{5,20}kb/      new N=500 backgrounds, gitignored, regenerable
  legacy_fdr01_n500/                pre-revision snapshot (caveat C3), 57 MB
  GO_grid/                          18 cells + INDEX.csv + MANIFEST.json; the 0.05 files
                                    are TRACKED (55 MB), the 0.1 retrieval twins gitignored
  TEs_on_genes_{5,10,20}kb.csv      per-window gene tables, gitignored, regenerable by 07a
svg/                                one SVG per panel + PLACEMENT.md
trackhub/                           105 MB, gitignored; published to gh-pages
```

## `revision_lib.py`

Constants: `CLASS_PALETTE` (`LINE #cc660b`, `LTR #70453c`, `SINE #ab1f20`, `DNA #195f90`,
`Retroposon #765297`, `RC #238023`), `GLOBAL_FONT_SIZE = 10`, `svg.fonttype="none"`,
`FDR_THRESHOLD = 0.05`, `N_PERMUTATIONS = 500`, `EMPIRICAL_P_FLOOR = 2/501`, `IFNA_DOMAIN`, and
resolved paths (`REPEATS_BED`, `GENES_BED`, `OBO_PATH`, `GAF_PATH`, `OUTPUT_DIR`, `SVG_DIR`).

| Function | Copied from | Changed how |
|---|---|---|
| `run_goatools_enrichment`, `run_goatools_ordered_enrichment` | GO nb cell 6 | `fdr_threshold` 0.1 → **0.05**; DAG/GAF loaded once per process and cached |
| `save_go_network_svg` | GO nb cell 6 | `+min_shared_genes`, `+max_term_genes`, `+title` (suppressible, G11), enforced collision check |
| `visualize_go_class_network` | GO nb cell 36 | same filters; palette switched from Tableau to the shared TE class palette |
| `save_go_network_svg_families_by_classes` | GO nb cell 175 | same filters; `family_to_class` is an explicit argument, not a notebook global |
| `CLASS_PALETTE` | GO nb cell 3 | verbatim |

New here, with no frozen original: `find_label_collisions` / `assert_no_label_collisions` /
`save_svg_collision_checked` (a colliding figure cannot be written silently — C16);
`find_offpage_labels` / `assert_labels_on_page` / `record_offpage_labels` and
`assert_svg_geometry` (the pinned-canvas guards: a panel cannot come out the wrong shape, and a
label cannot be cropped off a page that no longer grows to fit it);
`load_permutation_counts(window=…, seeds=…)` / `permutation_totals` / `read_counts_file` /
`load_manifest` (the compacted-store reader — **`n` is a weight column, not a row count**).

Both network plotters also gained presentation-only arguments, every one of them defaulting to
the published behaviour: `label_wrap`, `label_max_move`, `leader_lines`, `show_colorbar`,
`show_size_legend`, `edge_alpha`, `edge_width_cap`, `tight_bbox`, plus `show_jaccard_legend`,
`same_class_weight` and `legend_rect` on the families version. None of them changes a statistic;
`edge_width_cap` clips drawn stroke width only and the Jaccard legend is capped identically, so
the legend cannot become a lie (the Figure 7H lesson, avoided in advance).

The DAG/GAF cache is the one behavioural change beyond FDR and the new arguments: the frozen helper
reparses 31 MB + 190 MB *per call*, and the GO re-run makes ~1,200 calls (20 min instead of ~9 h).
GO results are unchanged — the same DAG and association dict reach `GOEnrichmentStudy` either way,
and the 0.1 retrieval cut reproduces the published tables exactly.

## Numbers worth not re-deriving

| Quantity | Value | Source |
|---|---|---|
| Permutations | **N = 500**; empirical p floor 2/501 = **0.004** | `permutation_counts_10kb/MANIFEST.json` |
| Length bias the permutation control removes | mean element length vs random OR, **Pearson R = 0.985** (n = 44, p = 8.4 × 10⁻³⁴, lengths 122–6,357 bp); Alu 316 bp → OR 1.535, L1 6,357 bp → 2.664 | `output/length_bias_correlation.json`. **Not 0.661** — an earlier `../CLAUDE.md` said that and was wrong |
| Convergence | at N = 250 the running mean is within 0.06 SD (worst class) / 0.10 SD (worst family) of its N = 500 value, SD within 4 %; at N = 100 drift is still 0.18 SD | `output/permutation_convergence_checkpoints.csv` |
| GO 0.1 → 0.05 | classes-by-count 504 → **425**; classes-by-divergence 516 → **414**; families **195 → 140**; subfamilies 1,231 → **1,003** (31 subfamilies lose every term) | `output/results_numbers.json`. Plan §3.2's family row (196 → 160) was wrong |
| Families with ≥ 1 term at 0.05 | **13 of 44** (Dong-R4 is the one lost) | same |
| IFNA domain (chr9:21,150,692–21,370,055) | 175 TEs, **77 L1** (44 %) from 36 subfamilies, 12 genes; mean L1 divergence **135.7** vs 188.2 genome-wide; 351 vs 181 L1/Mb = 1.9× | `output/ifna_qc.json` |
| IFNA tests | T1 p = 0.022 · T2 (≥ 40 L1) p = 0.0061 · T3 (≥ 10 genes) p = 0.0017 · T4 OR = 3.01, Fisher p = 3.2 × 10⁻⁶ | `output/ifna_test_results.csv` |
| Assembly ceiling | 208.75 Mb newly resolved (6.70 %), but only **0.41 % of TEs, 0.49 % of windows, 0.55 % of genes — and 0 bp in the IFNA domain** | `output/assembly_bound_summary.json` |
| Window concordance | classes ρ = 0.943 / 0.943 / **1.000**, 1 of 6 flips; families 0.891 / **0.828** / 0.941, 10 of 44 flips | `output/window_sensitivity_concordance.csv` |
| Percentile robustness | **8 of 9** headline claims survive at 10 %; 0.85–0.93 of published terms preserved | `output/percentile_sensitivity_headline.csv` |
| GO grid, 3 windows × 2 percentiles | 18 cells. Widening the **percentile** always finds more terms (9/9 level × window combinations); widening the **window** does not, and not even in the same direction: classes-by-count **falls** 510 → 425 → 299, classes-by-divergence **peaks at 10 kb** (209 / 414 / 303), families **rises** 137 → 140 → 201 | `output/go_grid_summary.json` |
| GO grid preservation | of the terms the paper reports, the fraction still significant in another cell is at worst **0.440**, median **0.677** — markedly weaker than the percentile-only result (0.85–0.93), so **the TSS window matters more to GO than the gene-set percentile does** | `output/go_grid_preservation.csv` |
| GO grid concordance | per-group term counts vs the published cell: lowest Spearman ρ **0.614**, every label-shuffling permutation p ≤ **0.022** | `output/go_grid_concordance.csv` |
| Headline claims across all six conditions | **3 of 9** survive all six; all 9 hold in the published cell. `SVA / termination of RNA polymerase II transcription` survives 2 of 6 and `hAT-Charlie / MHC class I protein complex` 1 of 6 | `output/go_grid_headline_by_condition.csv` |
| Group survival across all six | classes-by-count 5 of 7, classes-by-divergence 6 of 9, families 9 of 18 keep ≥ 1 term | `output/go_grid_group_survival.csv` |
| The grid is GO only | **no permutations were re-run for it**, so a difference between cells is a gene-set effect, not a background effect. The enrichment odds ratios of Table 1 and Figure 2 are unchanged (caveat S8) | — |
| Supplementary package | five workbooks, **8.7 MB** total: S1 TE-TSS map + enrichment (7.0 MB, 4 sheets), S2 gene sets (0.9 MB, 3), S3 GO at FDR 0.05 (0.2 MB, 3), S4 IFNA + prior work (0.2 MB, 7), S5 sensitivity (0.4 MB, 13) | `supplementary/INVENTORY.json` |
| Two April gene-set sheets were **empty** | `Supplementary File 3.xlsx`'s `TE top` and `TE bottom` sheets shipped with no rows. Reconstructed by 08 from the same construction 06 uses, after checking it against the six non-empty sheets: every disagreement sits exactly **on** the tie boundary (LINEs: 1,257 genes above the boundary count of 9, 1,033 tied at it for 179 places), and no gene above the boundary differs in any set | `supplementary/README.md` |
| Element counts per window | 293,652 (5 kb) / **582,540** (10 kb) / 1,157,235 (20 kb); merged span 144,952,895 / **272,233,268** / 494,969,139 bp | `output/window_sensitivity_ntss.json` |

Counts are **per-TSS, not per-gene**: a gene with several annotated TSS contributes several windows,
so an element within 10 kb of two TSS of one gene is counted twice. This is a property of the
published design, flagged in Limitations and reproduced deliberately everywhere here.

## Discrepancies documented rather than fixed

1. **`NUM_PERMUTATIONS = 1000`** still stands in the frozen download notebook (cell 34) and **was
   never executed**: 501 files existed, p floor 2/501 = 0.004, and the Methods were corrected *down*
   to 500 (C20). Wording to use verbatim: *"500 permutations were run; the generator's 1000 was never
   executed."*
2. **`GO_subfamilies.py` stays at `fdr_threshold=0.1`** because it belongs to the companion
   subfamilies manuscript, whose figures are not regenerated (C3). The revision's subfamily GO at
   0.05 comes from `06_go_rerun_fdr005.py --level subfamilies` → `output/GO_tables_fdr005/`; the
   shared `../GO_tables/` is untouched.
3. **The frozen `TEs_mapped_on_TSS_analysis.ipynb` still contains the broken `../T2T_article/…`
   path** and reads the now-deleted 6.37 GB `consolidated_random_data.csv`. Both are superseded, not
   fixed: use the local BED and `load_permutation_counts()`; `01c_expand_counts.py` rebuilds the
   legacy files if that notebook must be re-run.
4. **The legacy 10 kb background is not byte-reproducible.** `bedtools shuffle -seed N` is
   deterministic for a given binary but the December 2025 run used a different build, so
   `01a --check-legacy` tests **distributional** equivalence instead: every class Mann-Whitney and
   KS p > 0.15, largest standardised mean difference 0.16 legacy SD (RC), pooled KS D = 0.00031
   (`output/pipeline_equivalence_check.json`). That is the right standard because D1 means the
   published background is never regenerated.
5. **`flavone metabolic process` (GO:0051552) does not simply leave the paper.** Same raw p at two
   levels, different BH correction: **L1 family FDR 0.031 — kept**; LINE class FDR 0.088 — removed.
   Deleting it everywhere, as the plan said, would have dropped a still-significant result. The five
   overlapping genes are UGT1A6/7/8/9/10.

## Figures

44 SVGs in `svg/`, from six notebooks (`nb07` adds S12A–C); `svg/PLACEMENT.md` maps each to its Figma frame. Figma file
key `WRNeTzKZObdmAQ8QG1EZlq`, single page `0:1`, one frame per figure. **Scope any Figma work by node
ID, never by name** — the same canvas carries the subfamilies paper's labels and `Figure 9` exists
twice (C9). Frames: Fig 1 `856:7`, 2 `856:8`, 3 `856:9`, 4 `859:25`, 5 `861:28`, 6 `861:33`,
7 `861:34`, schematics `766:1208` → Fig 9 and `861:35` → Fig 10; S1–S3 `856:10/11/12`,
S4/S5 `860:26/27`, S6 `861:29` (promoted to new main Fig 8A), S7 `861:30`, S8 `861:32`.

**The figures were placed and exported on 2026-08-07 to `current_figures_260807/`, and layout
changed two things the SVG names still reflect.** The exports are authoritative:

- **The supplementary set is S1–S13, not S1–S14.** Old S6 became Figure 8A, so S7–S14 shifted
  down to S6–S13. `S13E_geneset_and_rank_stability.svg` is panel **S13A**; `S14A/B/C_*.svg` are
  panels **S13B/C/D**. The SVG filenames are deliberately **not** renamed — the notebooks write
  them and renaming breaks re-execution.
- **Figure 8 has four panels**: A browser view, B null distributions, C subfamily composition,
  **D** leave-one-out (promoted from an optional inset).
- **The colourbar situation is the reverse of `PLACEMENT.md` §4.4**: Supplementary S10 carries a
  colourbar and legends, Figure 6A carries neither.

`svg/PLACEMENT.md` §0 holds the full mapping and the open Figma items;
`figures_text_alignment_plan_260807.md` is the reconciliation with the manuscript.

Facts that bite when regenerating panels:
- **Figure 5 does not group by `class_name`.** Frozen cells 108/114 rewrite it to
  `"{class}_{highest|lowest}"` with a 21-key palette in cell 111. Grouping by `class_name` alone
  silently merges the highest- and lowest-divergence sets — the whole contrast the figure exists to
  show. Correct is **8 groups**; under FDR 0.05 the DNA class drops out of that level entirely.
- **Figure 6A carries 5 GO terms per group** after the 1.2× font increase — it was 9. Read the
  number off `output/network_qc.json`, never from memory: at 12 pt the collision checker cannot
  place 9 labels for 44 family groups at any canvas size the ladder tries, and the overlaps it
  rejects are real (93–348 px², not float noise). 4A and 5A keep 10, **S9 and S10 both keep 30**
  (S10 gained one over its earlier 29 on the pinned canvas). The legend must quote the achieved
  value (edit M10). `S11` is the one panel where the collision check is **waived**, now at **7
  overlapping pairs instead of 11** and still with all 30 terms per family, recorded in
  `output/network_qc.json`, never silent.
- **S9, S10 and S11 are drawn on a pinned canvas, not by the canvas ladder** (2026-08-07):
  815.04 × 1608.48 pt at aspect 0.507 for S9/S10 (half the published width, same height) and
  804.24 × 1267.92 pt at aspect 0.634 for S11, which is the aspect of reference frame `861:33`.
  They are saved with `bbox_inches` off so `figsize × 72` is exactly the page;
  `rl.assert_svg_geometry` refuses a panel more than 3 % off target and
  `rl.assert_labels_on_page` refuses a label that crosses the page edge, because a pinned page
  crops what a tight bbox would have grown to include. Two rungs replace canvas growth and both
  cost no information, so both are exhausted before a term is dropped: `WRAP_RUNGS` (label text
  wrapped at 26/22/18 characters — nothing here wrapped labels before) and `K_RUNGS`
  (`spring_layout` k at 0.6/0.9/1.2). S9 needed only k = 0.9, S10 k = 1.2, S11 wrap 18 + k 1.2.
  **S11 also drops `same_class_weight` 10.0 → 2.0**, which is why its GO circles used to clump,
  and carries no colourbar and no legends — reclaiming the legend strip is what buys the full
  width, and the legends then land on node labels; Figure 6A carries the same scales.
- **The review's "30–50 % denser packing" is not achievable for 4A/5A/6A and this is measured, not
  conceded** — but the first ladder never tried the two levers that do work, so do not read it as
  exhausted. A fixed number of labels needs a minimum area; 1.2× fonts inflate every label box by
  1.44× in area, so 30 % *less* area means the same labels in 0.49× the space they already needed.
  The ladder was inverted to try compaction first anyway, and it records
  `canvas_area_vs_baseline`, `compaction_target_met` and the best any compaction rung reached
  per panel. Raising the `adjust_text` forces was measured to make packing **worse** — the solver
  overshoots — and more `spring_layout` iterations changes almost nothing; neither is worth
  retrying. **Label wrapping and a raised k at fixed canvas are worth retrying** — they are what
  got S9/S10 clean at half the width with every term kept.
- **Nothing about these panels is byte-reproducible, and no md5 gate on them can pass.**
  `adjust_text` is stochastic (two renders with the same `PYTHONHASHSEED` differ) and matplotlib
  stamps a `<dc:date>` into every SVG. Compare with `svg/compare_panels.py` — page size, mark
  count and the exact set of label strings — not with `md5sum`. A consequence: re-running `nb06`
  re-lays out 4A/5A/6A as well, and on 2026-08-07 Figure 6A twice came out at 7 terms per family
  instead of the approved 5. The three main-text SVGs **and their `network_qc.json` records** were
  restored from the approved run so M10 still reads 5; taking the extra two terms means re-placing
  6A in Figma and updating M10.
- **Figure 7H filters its ribbons at ≥ 5 GO terms, and Supplementary S8C is the unfiltered
  twin.** The published Figure 7 caption already promised the filter; the first version applied
  it only to the count *labels*, so every ribbon was still drawn. Filtering hides 36 class →
  group and 50 group → family ribbons, 146 GO terms (`output/sankey_ribbon_filter.json`).
  Because the filter is visual only, bar heights stay unfiltered, so retained ribbons do not
  fill their bars — correct, and the caption says so. **`Figure S8C` is not `Figure 8C`** (the
  IFNA subfamily composition panel); every reference must carry the `S`.
- **Figure 7B/E/F/G are half-size and their statistics live in the caption**, not the artwork
  (decision D-a, edit M9). `nb06` prints the four sentences to paste. If M9 is skipped those
  panels report a test with no p-value anywhere.
- **Figure 7 is entirely family-level** and Supplementary 8B is one combined class+family clustermap
  (18 × 24 at 0.05). Both were misidentified once; swap all eight Figure 7 panels together.
- **Figure 2D–2F are family-level, and D = obs/random OR, E = observed OR, F = random OR.**
  The first version plotted the 1,129 subfamilies and had D/E swapped; the published caption is
  unambiguous on both counts. At family level the panels reproduce the three published
  significance sentences exactly (2D: SINE–DNA only; 2E: LINE–DNA and SINE–DNA; 2F: LINE–SINE,
  LINE–DNA, LTR–SINE) and `nb03` cell 17 asserts it. Retroposon and RC have one family each, so
  they are points without a box and their 9 class pairs read `n/a`, not `ns`.
  Figure 7's A–H are still inferred from content — verify against the frame before swapping.
- No `fig.suptitle(...)` in any exported panel (G11); colorbars ship as vector
  (`Fig456A_colorbar_vector.svg`, `Fig4B_colorbar_vector.svg`) to replace pasted bitmaps (G12).

## The manuscript files

All three live in `Revised_manuscript/`: the read-only submitted baseline
(md5 `1dbcbd4419987fd28ddf803129487cfd`), the tracked revision, and the extended-discussion
supplementary file. Citations are **Mendeley Cite content controls** — 128 `<w:sdt>` tagged
`MENDELEY_CITATION` plus one `MENDELEY_BIBLIOGRAPHY`, payload in
`word/webextensions/webextension1.xml`, **no** legacy field codes.

Consequences: any tracked edit must operate on whole `<w:sdt>` elements — never rewrite a paragraph
from its concatenated text, which is exactly what the `word_rewrite` skill's `tracked_replace` does
and why `13_…` implements `tracked_replace_safe` / `delete_paragraph_safe`. Matching must tolerate
non-breaking and zero-width spaces, curly quotes, en dashes and superscript exponents stored as
plain digits (`10-40`, not `10⁻⁴⁰`), and must strip `<w:proofErr>` markers. A second document must be
built by **copying** the manuscript and deleting, never by pasting into a fresh file — the payload
lives in a part a fresh document does not have.

**The working file changed on 2026-08-04 and with it this script's contract.** Daniil edited the
revision by hand and **accepted 594 of the 1,124 tracked deletions**, so
`T2T_genes_article_G3_revision_260804_manual.docx` — not the baseline, and not the 260803 output —
is the current state of the paper. `13_…` and `15_…` both target it, `13_…` edits **in place**, and
`--from-baseline` refuses to overwrite it. `T2T_genes_article_G3_revision_260803.docx` is kept as
the record of what the scripts produced before that pass: do not edit or delete it.
`14_build_extensive_discussion.py` still builds from the pristine baseline and is unaffected.

**The current state of the paper moved again on 2026-08-07: it is now
`T2T_genes_article_G3_revision_260807.docx`.** `16_figure_alignment_edits.py` aligned the text
with the exported figures, and unlike `13_…` and `15_…` it does **not** edit in place — it reads
the 260804 file and writes a new one, so 260804 joins 260803 and the baseline as read-only input
of record (its md5 is asserted unchanged by the verification block). Consequences: any further
manuscript script targets **260807**, and `16_…` may be re-run at will because its output is a
deterministic rebuild from its input rather than an accumulation of in-place edits. All 129
content controls survive; the script asserts the count before and after.

**The current state of the paper moved again on 2026-08-09: it is now
`T2T_genes_article_G3_revision_260809.docx`.** `18_review_edits_260809.py` applied the
`scientific-review` skill to `…_260807_manual.docx` (Daniil's hand-edited 260807 file, which joins
the read-only inputs of record). Like `16_…` it is a deterministic rebuild from its input, so it
may be re-run at will; 48 operations, all reported, none silent. What it changed:

- the three inline `'Daniil to Claude:'` comments, resolved by doing what they ask and removed;
- **three subsections rewritten** — the interferon-alpha domain as a self-contained five-passage
  analysis, the prior-work comparison citing File S4 sheet by sheet plus a new **Table 3**, and the
  sensitivity subsection citing File S5 sheets and individual S11–S13 panels;
- **four numbers that did not reproduce** from the source tables (see the table in
  `../Daniil_review_of_Claude_examples_principles_260809.md` §IV.6 and below);
- **supplementary figures renumbered by first citation** — old S8→S6, S9→S7, S6→S8, S7→S9,
  S12→S11, S13A→S12, S11→S13, and old S13B–D merged into S1 as panels A–C with the length ridge
  plots becoming S1D/E. All 13 legends rewritten and physically reordered. The Figma-side moves are
  in `../Figures_renaming_260809.md`; nothing here writes to Figma.
- the seven remaining `[EDITOR NOTE]` blocks and the `(ref)` placeholder.

After it, every series is cited in ascending order: Figures 1–10, Figures S1–S13, Tables 1–3,
Files S1–S6. The supplementary package it cites is **`supplementary_260809/`**, built by
`17_build_supplementary_260809.py`.

**Two docx mechanics this pass discovered, both now guarded in the script:**

1. **`delete_paragraph` is a silent no-op on a paragraph made of `<w:ins>`.** It only wraps
   direct-child runs, so deleting an editor note or a legend that an earlier pass had *inserted*
   left it fully visible in the Accept All view. `delete_paragraph_smart` drops our own unaccepted
   insertions outright — which is what Word does and is invisible to Reject All — and refuses to
   touch another author's.
2. **An anchor that straddles an existing `<w:del>` relocates it.** `safe_tracked_replace` rebuilds
   the span from the first matched run, so a pre-existing deletion sitting between two matched runs
   survives but moves, and Reject All then returns that already-deleted fragment in the wrong
   place. `spans_revision()` refuses such an anchor; shorten it until it lies inside one run.

**The current state of the paper moved again on 2026-08-10: it is now
`T2T_genes_article_G3_revision_260810.docx`.** `19_review_edits_260810.py` is the final review pass,
run against the re-exported `current_figures_260810/` and `supplementary_260809/`. Like `16_…` and
`18_…` it is a deterministic rebuild from its input, so it may be re-run at will; 68 operations, all
reported. Its findings are in `review_report_260810.md`. Measured: `<w:sdt>` 138 → 138,
`MENDELEY_CITATION` 137 → 137, Reject All byte-identical to the 260809 input, 0 duplicate revision
ids, 0 notes in the docx. The P3 numbering audit went from FAIL to PASS.

Six facts from that pass worth not rediscovering:

1. **The class-level summary sentence quoted the ERVK family row for "all TEs".** 1.9509 / 1.7782 /
   1.0971 is ERVK in `enrichment_families_with_random.csv`. For all TEs the numbers are **1.95 /
   1.79 / 1.086**: 582,540 of 3,709,429 elements in a window against a 500-seed mean of 542,973.1,
   and the enrichment score is a ratio of count ratios so it is independent of the window-geometry
   constant. 1.94 / 1.78 / 1.097 was also arithmetically impossible.
2. **Five GO terms in the divergence subsection sat in the rejected 0.05 < FDR ≤ 0.1 band** and were
   reported as results (spindle localization 0.071, fatty acid catabolic process 0.080, cell-cell
   adhesion 0.066, cell population proliferation 0.052, VLDL particle 0.056), plus one term that
   does not exist at that level ("spermatogenesis") and one that is not a GO term at all ("synaptic
   chemical transition" → *chemical synaptic transmission*). Term lists inherited from the submitted
   baseline were never re-checked against FDR 0.05 — check any that remain.
3. **The document carried a duplicate Table 1 and Table 2** after Funding — the 10_tables.py pair
   that Daniil superseded, left as unaccepted insertions — and a 7 × 11 table whose text was struck
   but whose rows were not, so Accept All kept an empty grid. **A table row is only removed by Accept
   All when its `w:trPr` carries a `w:del`;** striking the `w:t` is not enough.
4. **`word_rewrite_trackchanges` numbers revisions from a fixed base of 1000**, so a fourth pass over
   the same file produces duplicate `w:id` values. Set `wrt._rid[0]` above the highest id already in
   the input before editing.
5. **28 orphan citation-only paragraphs carry 83 references, 68 of them cited nowhere else.**
   Deleting them takes the reference list from 126 to ~58 entries and renumbers the whole paper.
   That is correct — `Extensive_discussion_260803.docx` (File S6) carries those references in its own
   bibliography — but it must be one pass: delete all 28, refresh Mendeley, re-read for moved numbers.
6. **G3 allows an abstract under 250 words**, so the 232-word abstract is compliant;
   `nar_review_tools validate` flags it only because it hard-codes NAR's 200.

**`Reject All` no longer restores the baseline.** Accepted deletions are permanent by design, so
the clean-plus-tracked submission pair is produced from the working file forward, and the diff a
reviewer sees covers the edits made *after* Daniil's acceptance pass rather than the whole
revision. If the journal wants a full tracked diff against the submitted version, produce it as a
Word **Compare** of the baseline against the final clean file, not from the revision marks.

Measured state of all four documents (`<w:sdt>` / `<w:ins>` / `<w:del>` / editor notes): baseline
129 / 0 / 0 / 0; 260803 output 129 / 395 / 1,124 / 10; **working file 129 / 380 / 530 / 9**. All
129 Mendeley content controls survive in the working file, so nothing is broken.

Earlier state, for the 260803 file: 56/56 Phase 5 edits + 21/21 Phase 7 edits applied,
no `<w:t>` inside a `<w:del>`, all 129 controls intact.
Two known limits: the two section relocations are structural moves rather than Word tracked moves,
so Reject All restores the original text but not the original order (validation compares paragraphs
as a multiset); and `run_blocks` cannot see inside a `<w:ins>`, so edits to our own inserted text
need the insertion-aware pass.

## Verification

```bash
md5sum ../*.ipynb ../GO_subfamilies.py                 # match the four sums above
grep -rn "nbimporter\|import_ipynb" revision_G3/       # expect no hits
python revision_G3/01c_expand_counts.py --seed 1       # lossless round-trip
python revision_G3/11_results_numbers.py               # re-derive every quoted number
```

Plan §10 has the full sweep — the Mendeley integrity check, the manuscript grep block (no
`T2T_genes_evolution`, no `1000 random`, no `FDR < 0.1`, no class-level `flavone`) and the track-hub
range-request check. `project_overview.md` §7 lists it with expected values.

## What is still open

The figures were placed, renamed per `../Figures_renaming_260809.md` and re-exported to
`current_figures_260810/` on 2026-08-10. All of §1 and §2 of that plan verified against the PDFs, and
two of the five §3 corrections landed (Figure 4's panel A, the colourbars on 4A/5A/6A/S7). **Three
artwork items remain: Figure 8 has no panel letter C** although the text cites Figure 8C twice,
Figure 3's A overlaps its y-axis tick `1.0`, and Figure 8C writes `3.2*10-6` where the text writes
`3.2 × 10-6`. `review_report_260810.md` §4 has the measured detail.

Everything computational is done; the remainder is Daniil's, one action each — the three figure fixes
above, publishing `trackhub/` to `gh-pages`,
the Zenodo DOI, the Mendeley style switch and preprint metadata, refreshing Mendeley in both docx
files, deleting the 28 orphan citation paragraphs, adding the five missing prior-work citations
(`review_report_260810.md` §3.2), the 10 browser-verify items, and the response letter plus
submission package. Two decisions
and two cheap improvements are itemised in `project_overview.md` §8.
