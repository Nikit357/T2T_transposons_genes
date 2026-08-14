# Supplementary material — G3-2026-406828

*Telomere-to-telomere co-mapping of transposable elements and human genes identifies a cluster of young L1 elements in the interferon-alpha domain*

Assembled by `revision_G3/08_build_supplementary.py`. Five thematic workbooks, one sheet per table, each workbook opening with a README sheet that describes its sheets.

## What changed relative to the April 2026 submission

1. **Every Gene Ontology result is now at FDR < 0.05** (Benjamini-Hochberg), with no "suggestive" band. The April supplementary files were at FDR 0.1. This is the change reviewer minor comment 2 asked for.
2. **Fourteen candidate files became five workbooks**, so every `File Sn` citation in the manuscript is renumbered and now names a sheet as well as a file. The mapping from the old numbering is below.
3. **The gene-set tables are long format** rather than one sheet per group. This also restores the `hAT?` and `hAT-Tip100?` family names, which Excel had mangled into sheet names in the April File S7.
4. **Workbook S1 carries both the family and the subfamily enrichment tables.** The April File S2 caption described subfamilies but the file contained the 44 families; both are now present with the full column set.

## Old numbering -> new

| April 2026 | Now |
|---|---|
| File S1 (TSS/TE coordinates) | File S1, sheet `TSS_TE_intersections` |
| File S2 (enrichment) | File S1, sheets `enrichment_families` / `enrichment_subfamilies` |
| File S3 (gene sets by TE group) | File S2, sheet `by_TE_group` |
| File S4 (GO by TE group) | File S3, sheet `GO_TE_groups` |
| File S5 (gene sets by divergence) | File S2, sheet `by_divergence` |
| File S6 (GO by divergence) | File S3, sheet `GO_by_divergence` |
| File S7 (gene sets by family) | File S2, sheet `by_family` |
| File S8 (GO by family) | File S3, sheet `GO_by_family` |
| — (Lu et al. overlap) | File S4, sheet `prior_work_overlap_matrix` |
| — ("the accompanying tables") | File S5 |

## The workbooks

### `File_S1_TE_TSS_map_and_enrichment.xlsx`

The TE-TSS map and the enrichment statistics.

| Sheet | Rows | Columns | Contents |
|---|---|---|---|
| `TSS_TE_intersections` | 38,704 | 22 | every TSS window with the transposable elements within 10 kb of it: counts and mean divergence per class, and the comma-joined element lists. Counts are per TSS window, not per gene, so a gene with several annotated TSS contributes several rows |
| `enrichment_classes` | 6 | 17 | the six TE classes: observed and random odds ratios, raw AND adjusted Fisher p, and the empirical permutation p (N = 500) |
| `enrichment_families` | 44 | 20 | the 44 TE families, same columns |
| `enrichment_subfamilies` | 1,143 | 17 | all 1,143 TE subfamilies, same columns — this is the table the File S2 caption described |

### `File_S2_gene_sets.xlsx`

The foreground gene sets used for GO enrichment.

| Sheet | Rows | Columns | Contents |
|---|---|---|---|
| `by_TE_group` | 9,708 | 3 | top 5 % of genes by element count for each TE group, plus the TE-top and TE-bottom sets (was File S3, 8 sheets) |
| `by_divergence` | 12,226 | 3 | highest- and lowest-divergence 5 % of genes per class (was File S5, 10 sheets) |
| `by_family` | 26,551 | 3 | top 5 % of genes by element count for each of the 44 families (was File S7, 44 sheets; the hAT? and hAT-Tip100? names are restored here — Excel had mangled them into sheet names) |

### `File_S3_gene_ontology.xlsx`

Gene Ontology enrichment at FDR < 0.05.

| Sheet | Rows | Columns | Contents |
|---|---|---|---|
| `GO_TE_groups` | 425 | 12 | GO terms enriched in each TE group's gene set, with the manual functional-group classification (was File S4, at FDR 0.1) |
| `GO_by_divergence` | 414 | 13 | GO terms for the highest- and lowest-divergence gene sets (was File S6, at FDR 0.1) |
| `GO_by_family` | 140 | 13 | GO terms per TE family (was File S8, at FDR 0.1) |

### `File_S4_IFNA_domain_and_prior_work.xlsx`

The interferon-alpha domain, the assembly bound, and the overlap with prior work.

| Sheet | Rows | Columns | Contents |
|---|---|---|---|
| `IFNA_elements` | 175 | 8 | all 175 transposable elements in the 220 kb interferon-alpha domain (chr9:21,150,692-21,370,055) |
| `IFNA_tests` | 4 | 21 | the four domain tests: observed value, null mean and SD, and the empirical p-value |
| `IFNA_subfamily_composition` | 36 | 6 | the L1 subfamilies represented in the domain and their divergence |
| `assembly_bound` | 75 | 9 | newly resolved T2T sequence and the share of elements, windows and genes it contributes, by TE class, family and chromosome |
| `prior_work_overlap_matrix` | 33 | 14 | Lu et al. 2020 category x TE group: overlap size, Fisher p and Jaccard |
| `prior_work_categories` | 5,960 | 5 | Lu et al.'s mouse categories mapped to human genes through MGI orthology |
| `prior_work_shared_genes` | 2,195 | 3 | the individual genes shared between each category and each group |

### `File_S5_sensitivity_and_robustness.xlsx`

Window-size, percentile and GO-grid sensitivity analyses.

| Sheet | Rows | Columns | Contents |
|---|---|---|---|
| `window_classes` | 18 | 15 | class-level enrichment at 5, 10 and 20 kb |
| `window_families` | 131 | 15 | family-level enrichment at 5, 10 and 20 kb |
| `window_concordance` | 6 | 13 | Spearman rho per window pair with a label-shuffling permutation test |
| `window_flips` | 11 | 5 | every group whose significance call changes between windows |
| `percentile_summary` | 28 | 13 | GO term counts and stability per group, top/bottom 5 % vs 10 % |
| `percentile_terms` | 1,046 | 7 | every GO term gained or lost at 10 %, named |
| `geneset_stability` | 21 | 10 | overlap of the top-5 % gene sets between windows, with a hypergeometric p-value |
| `rank_stability` | 21 | 7 | Kendall tau of the per-gene ranking between windows, bootstrap CI |
| `GO_grid_index` | 18 | 11 | one row per cell of the 3 windows x 2 percentiles GO grid |
| `GO_grid_preservation` | 18 | 11 | Jaccard and fraction of published terms preserved, per cell and level |
| `GO_grid_terms` | 4,996 | 9 | every GO term gained or lost in each cell relative to 10 kb / 5 % |
| `GO_grid_concordance` | 15 | 9 | Spearman of per-group term counts against the published cell, with a label-shuffling permutation test |
| `headline_by_condition` | 9 | 18 | each headline claim under all six conditions |

## Conventions that apply throughout

* **Element counts are per TSS window, not per gene.** A gene with several annotated TSS contributes several windows, so an element within 10 kb of two TSS of the same gene is counted twice. This is a property of the published design and is stated in the manuscript's Limitations.
* **GO FDR is 0.05** everywhere, Benjamini-Hochberg corrected.
* **The permutation background is N = 500**, so the empirical p-value floor is 2/501 = 0.0040.
* **Raw and adjusted p-values are both reported** wherever a correction was applied, per G3's statistics policy.

## One column is not in the workbooks

`Full Term Gene List` — up to 125,086 characters per cell (the full human annotation of 'protein binding'), 3.8x Excel's 32,767-character cell limit. It is a property of the GO annotation (go-basic.obo + goa_human.gaf, both cited in the Methods), not a result of this study. 'Overlapping Genes', which is the result, is kept in full.
## Two gene sets were empty in the April submission

The April `Supplementary File 3.xlsx` shipped **2 empty sheets**: `TE top`, `TE bottom`. Those gene sets drive GO results the paper reports and gene tracks the UCSC hub ships, so the sets themselves were never in doubt — only their inclusion in the supplementary file. They are reconstructed here by the same construction `revision_G3/06_go_rerun_fdr005.py` uses (top and bottom 1,436 genes by total element count), which was checked against the six non-empty April sheets first.

That check is not an exact match, and the reason matters: a top-1,436 cut lands inside a large block of genes carrying the *same* element count — for the LINEs set, 1,257 genes are strictly above the boundary count of 9 while 1,033 tie at it, competing for the remaining 179 places. Which tied genes make the cut is decided by row order and is arbitrary. The verification therefore requires every disagreeing gene to sit exactly **on** the boundary count, and it does; no gene above the boundary differs, in any set. The same arbitrariness applies to the two reconstructed sets, and to the published ones.


## Figures

The figures are exported from Figma by hand and are listed in `figures/PLACEHOLDERS.md`.
