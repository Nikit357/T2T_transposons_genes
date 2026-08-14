# Supplementary material — G3-2026-406828 (2026-08-09)

*Telomere-to-telomere co-mapping of transposable elements and human genes identifies a cluster of young L1 elements in the interferon-alpha domain*

Built by `revision_G3/17_build_supplementary_260809.py` from the 2026-08-05 package
(`revision_G3/supplementary/`), which is left untouched. Five thematic workbooks, one sheet per
table, each workbook opening with a README sheet that describes its sheets.

## What changed relative to the 2026-08-05 package

**File numbers S1–S6 are unchanged.** After the 2026-08-09 rewrite added File S4 citations to the
interferon-alpha and prior-work subsections, the files are cited in ascending order (S1, S2, S3,
S4, S5, S6), so no renumbering was required. Three things did change.

1. **File S3 gains a `GO_borderline` sheet — 236 rows.** The Results name several GO terms as
   narrowly *not* significant and quote their FDR values (flavone metabolic process at 0.088,
   glutamatergic synapse and positive regulation of lipopolysaccharide-mediated signaling at
   0.086, the MIR metal-ion terms at 0.078). Every shipped GO sheet is cut at FDR < 0.05, so
   none of those values could be checked anywhere in the package. The band
   0.05 < FDR ≤ 0.1 is now shipped for all three levels of analysis.
   *These terms are not results.* They are supplied only so that every FDR value printed in the
   manuscript is traceable.

2. **File S4 sheets are reordered to the order the text cites them, and the file is renamed** from
   `File_S4_IFNA_domain_and_prior_work.xlsx` to
   `File_S4_IFNA_domain_prior_work_and_assembly.xlsx`. `assembly_bound` previously sat between
   the two Results subjects; it is cited in the Discussion and now comes last.

3. **File S5: `headline_by_condition` moved** next to the percentile sheets, which is where the
   rewritten sensitivity subsection cites it.

## The workbooks

### `File_S1_TE_TSS_map_and_enrichment.xlsx` — 6.96 MB, 5 sheets

| Sheet | Rows | Cols | Contents |
|---|---|---|---|
| `TSS_TE_intersections` | 38,704 | 22 | every TSS window with the elements within 10 kb: counts and mean divergence per class, and the comma-joined element lists |
| `enrichment_classes` | 6 | 17 | the six TE classes: observed and random odds ratios, raw **and** adjusted Fisher p, empirical permutation p (N = 500) |
| `enrichment_families` | 44 | 20 | the 44 TE families, same columns |
| `enrichment_subfamilies` | 1,143 | 17 | all 1,143 subfamilies, same columns |

### `File_S2_gene_sets.xlsx` — 0.93 MB, 4 sheets

| Sheet | Rows | Cols | Contents |
|---|---|---|---|
| `by_TE_group` | 9,708 | 3 | top 5 % of genes by element count per TE group, plus the TE-top and TE-bottom sets |
| `by_divergence` | 12,226 | 3 | highest- and lowest-divergence 5 % of genes per class |
| `by_family` | 26,551 | 3 | top 5 % of genes by element count for each of the 44 families |

### `File_S3_gene_ontology.xlsx` — 0.26 MB, 5 sheets

| Sheet | Rows | Cols | Contents |
|---|---|---|---|
| `GO_TE_groups` | 425 | 12 | GO terms enriched in each TE group's gene set, with the manual functional-group classification |
| `GO_by_divergence` | 414 | 13 | GO terms for the highest- and lowest-divergence gene sets |
| `GO_by_family` | 140 | 13 | GO terms per TE family |
| **`GO_borderline`** | **236** | **11** | **new.** Every GO term with 0.05 < FDR ≤ 0.1, at all three levels: 79 classes-by-count, 102 classes-by-divergence, 55 families. Not results — supplied so the FDR values the Results name as narrowly non-significant can be checked |

### `File_S4_IFNA_domain_prior_work_and_assembly.xlsx` — 0.21 MB, 8 sheets

Sheets are in the order the text cites them: the interferon-alpha domain (Results), the
prior-work comparison (Results), the assembly bound (Discussion).

| Sheet | Rows | Cols | Contents |
|---|---|---|---|
| `IFNA_elements` | 175 | 8 | all 175 elements in the 219,363 bp domain (chr9:21,150,692–21,370,055) |
| `IFNA_tests` | 4 | 21 | the four domain tests: observed value, null mean and SD, empirical p |
| `IFNA_subfamily_composition` | 36 | 6 | the 36 L1 subfamilies in the domain and their divergence |
| `prior_work_overlap_matrix` | 33 | 14 | Lu et al. 2020 category × TE group: overlap size, fold over expected, Fisher p (raw and adjusted), Jaccard, overlap coefficient |
| `prior_work_categories` | 5,960 | 5 | Lu et al.'s mouse categories mapped to human genes through MGI orthology |
| `prior_work_shared_genes` | 2,195 | 3 | the individual genes shared between each category and each group |
| `assembly_bound` | 75 | 9 | newly resolved T2T sequence and the share of elements, windows and genes it contributes, by class, family and chromosome |

### `File_S5_sensitivity_and_robustness.xlsx` — 0.37 MB, 14 sheets

Sheets are in the order the rewritten sensitivity subsection cites them: window size, then gene
sets, then the percentile cut-off, then the GO grid.

| Sheet | Rows | Cols | Contents |
|---|---|---|---|
| `window_classes` | 18 | 15 | class-level enrichment at 5, 10 and 20 kb |
| `window_families` | 131 | 15 | family-level enrichment at 5, 10 and 20 kb |
| `window_concordance` | 6 | 13 | Spearman rho, Pearson r and Bland-Altman per window pair, with a label-shuffling permutation test |
| `window_flips` | 11 | 5 | every group whose significance call changes between windows |
| `geneset_stability` | 21 | 10 | overlap of the top-5 % gene sets between windows, with a hypergeometric p-value |
| `rank_stability` | 21 | 7 | Kendall tau of the per-gene ranking between windows, with a bootstrap CI |
| `percentile_summary` | 28 | 13 | GO term counts and stability per group, top/bottom 5 % vs 10 % |
| `percentile_terms` | 1,046 | 7 | every GO term gained or lost at 10 %, named |
| `headline_by_condition` | 9 | 18 | each headline claim under all six window × cut-off conditions |
| `GO_grid_index` | 18 | 11 | one row per cell of the 3 windows × 2 percentiles grid |
| `GO_grid_preservation` | 18 | 11 | Jaccard and fraction of published terms preserved, per cell and level |
| `GO_grid_concordance` | 15 | 9 | Spearman of per-group term counts against the published cell, with a label-shuffling permutation test |
| `GO_grid_terms` | 4,996 | 9 | every GO term gained or lost in each cell relative to 10 kb / 5 % |

**Two term-count definitions coexist in File S5, on purpose.** `GO_grid_index.n_terms_005`
counts group–term rows (the 10 kb / 5 % cell has 425); `GO_grid_preservation.n_terms` counts
*distinct* terms (375 for the same cell). Both are correct for their question, and the direction
of every trend the Results report is the same under either.

## Conventions that apply throughout

* **Element counts are per TSS window, not per gene.** A gene with several annotated TSS
  contributes several windows, so an element within 10 kb of two TSS of the same gene is counted
  twice. This is a property of the published design and is stated in the Limitations.
* **GO FDR is 0.05** everywhere, Benjamini-Hochberg corrected. The one exception is the new
  `GO_borderline` sheet, whose entire purpose is the band above that threshold.
* **The permutation background is N = 500**, so the empirical p floor is 2/501 = 0.0040.
* **Raw and adjusted p-values are both reported** wherever a correction was applied.
* **Class naming.** The manuscript calls the two smallest classes **SVA** and **Helitrons**.
  `File S5`'s `geneset_stability` and `rank_stability` sheets, and `File S4`'s `assembly_bound`,
  carry the RepeatMasker class names for the same two classes: **Retroposon** = SVA and
  **RC** = Helitrons. The manuscript states this correspondence where it first cites those sheets.

## One column is not in the workbooks

`Full Term Gene List` — up to 125,086 characters per cell (the full human annotation of *protein
binding*), 3.8× Excel's 32,767-character cell limit. It is a property of the GO annotation
(go-basic.obo + goa_human.gaf, both cited in the Methods), not a result of this study.
`Overlapping Genes`, which is the result, is kept in full.

## Integrity

`CHECKSUMS.sha256` and `INVENTORY.json` are regenerated by the build script and record the
sha256, byte size and per-sheet row/column counts of every workbook. Re-check with

```bash
~/venvs/Retroelements_3_11/bin/python revision_G3/17_build_supplementary_260809.py --verify
```
