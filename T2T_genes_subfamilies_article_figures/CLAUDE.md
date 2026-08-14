# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Directory Purpose

This directory stores the finalized publication figures (PDF + PNG) and the article manuscript
(`.docx`) for the **T2T genes subfamilies** paper. There is no runnable code here — all analysis
code lives in the parent directory `../` (see `../CLAUDE.md`).

Each figure exists in two formats: a vector PDF (print-quality) and a rasterized PNG (web/preview).

## Two papers, and which one this directory describes

**The figure inventory below is the subfamilies paper's, not the G3 paper's.** This was ambiguous
before and is worth stating plainly, because the two papers share this directory's history and a
Figma canvas:

| | Subfamilies paper | G3 paper |
|---|---|---|
| Level | 1,143 TE **subfamilies** | 6 classes, 44 **families** |
| Status | in preparation for PLOS ONE | **conditionally accepted at *G3*** (G3-2026-406828) |
| Figures | the inventory in this file | described in `../revision_G3/svg/PLACEMENT.md` |
| Manuscript | `T2T_genes_subfamilies_article.docx`, here | `../revision_G3/Revised_manuscript/T2T_genes_article_G3_revision_260803.docx` |

**The G3 manuscript files used to live here and have moved.** `T2T_genes_article_for_plos_one.docx`
was renamed `T2T_genes_article_G3_revision_260803.docx` and now sits in
`../revision_G3/Revised_manuscript/` together with the read-only submitted baseline
(`T2T_genes_article_G3_submitted_baseline_260418.docx`) and the extended-discussion supplementary
file (`Extensive_discussion_260803.docx`). Only `T2T_genes_subfamilies_article.docx` — this paper's
manuscript — remains here.

## Figures are composed in Figma, and the `Figure 9` name collides

Python and the notebooks produce **panels only, as SVG**; final figures are composed, lettered and
exported by hand in Figma (file key `WRNeTzKZObdmAQ8QG1EZlq`, single page `0:1`). Both papers'
frames share that one canvas, and **`Figure 9` exists twice on it** — once as the G3 paper's frame
`861:35` and once as a loose label belonging to this paper. Any Figma work must therefore be scoped
**by node ID, not by frame name** (caveat C9 of `../G3_revision_implementation_plan_260803.md`).
Retitling the G3 paper also implies retitling this one, since the two titles were parallel.

## Figure Contents

### Main Figures

| File | Content |
|------|---------|
| **Figure 1** | Multi-panel scatter plots: TE subfamily characteristics (divergence, enrichment, length) with Pearson correlations by TE class; includes box plots with pairwise statistical comparisons between classes |
| **Figure 2** | Three-panel scatter: (A) family-level divergence vs. odds ratio; (B) family-level with GO term grouping; (C) subfamily-level divergence vs. average TE count per TSS, with individual subfamily labels and enrichment significance markers |
| **Figure 3** | GO functional enrichment network (force-directed) for SVA/Alu subfamilies — diamond nodes = GO terms, circle nodes = gene sets; colored by TE class affiliation |
| **Figure 4** | GO functional enrichment network for LINE/MIR subfamilies — two sub-networks showing distinct functional clusters |
| **Figure 5** | (A) Jaccard index legend and example GO term network; (B) Clustermap heatmap of GO functional categories × TE subfamilies; (C) Clustermap of GO terms grouped by functional metagroups |
| **Figure 6** | "Ring of Power" full network (black background): all TE subfamilies as nodes connected by Jaccard gene-set similarity ≥ 0.025; colored by class; two connected components |
| **Figure 7** | (Top) Correlation heatmap across L1 subfamilies; (Middle) Scatter plot of L1 subgroup enrichment (LIM, LIP, LIF, LIH subtypes) with box plots; (Bottom) Network of L1/LINE subfamily relationships with GO term annotations |
| **Figure 8** | Scatter: average divergence vs. GO terms count for L1 subfamilies (Pearson R = −0.285, p = 0.050); shows young L1 elements near TSS yield more GO enrichment hits |

### Supplementary Figures

| File | Content |
|------|---------|
| **Supplementary Figure 1** | Full heatmap of enrichment statistics for all 1,143 TE subfamilies: observed odds ratio, random odds ratio, obs/random ratio, SD of random, total TE count, TSS-linked TE count, average divergence, Fisher adjusted p-value, random control p-value |
| **Supplementary Figure 2** | Chord (Circos) diagram: flow from TE classes to their constituent subfamilies, arc width proportional to subfamily size |
| **Supplementary Figure 3** | Two large heatmaps (first and second halves of subfamilies): TE count distribution per TSS neighborhood per subfamily, showing sparsity pattern |
| **Supplementary Figure 4** | Gene-subfamily network (black background, circular layout): all genes as nodes colored by TE class, showing co-occurrence structure across the full gene set |
| **Supplementary Figure 5** | 9-panel scatter grid + box plots: Spearman correlations between multiple enrichment metrics, stratified by TE class |
| **Supplementary Figure 6** | (Top) Scatter plots with Spearman correlations for additional metric pairs by TE class; (Bottom) Sankey diagram: TE class → subfamily → chromosome flow |

## Regenerating Figures

Figures are generated by notebooks and scripts in `../`:

- `../TEs_mapped_on_TSS_analysis.ipynb` → Figures 1, 2, Supplementary Figures 1, 2, 3, 5, 6
- `../Gene_ontology_analysis.ipynb` → Figures 3, 4, 5
- `../genes_subfamilies_network.py` → Figure 6 (`JACCARD_THRESHOLD = 0.025`)
- `../genes_subfamilies_network_clusters.py` → Figure 7 (L1 subnetwork, MCL clustering)
- Figure 8 is the L1-specific divergence vs. GO count scatter, produced in `../Gene_ontology_analysis.ipynb`

After regenerating, copy output SVG/PNG files here and export PDFs.

## Article Manuscript

`T2T_genes_subfamilies_article.docx` — the working manuscript. Contains figure captions that authoritatively describe each figure.
