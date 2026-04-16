# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a genomics research project analyzing transposable elements (TEs) and their proximity to genes in the T2T-CHM13 (Telomere-to-Telomere CHM13) human genome assembly. The analysis maps 3.7 million TEs across ~28,000 genes, calculates enrichment statistics for 44 TE families/subfamilies, and performs Gene Ontology (GO) functional enrichment to understand TE-gene co-localization and evolutionary arms races (with a focus on immune gene clusters like IFNA).

## Running the Analysis

### Notebooks (run in order)
```bash
jupyter notebook download_and_process_files_UCSC_genes.ipynb   # Data acquisition
jupyter notebook TEs_mapped_on_TSS_analysis.ipynb               # Core proximity mapping
jupyter notebook Gene_ontology_analysis.ipynb                   # GO enrichment
```

### Python scripts
```bash
python GO_subfamilies.py                    # Per-subfamily GO enrichment → GO_tables/
python genes_subfamilies_network.py         # "Ring of Power" network visualization
python draw_length_divergence_corr.py       # TE length vs divergence correlation plots
```

### Permutation test (statistical background)
```bash
bash run_permutation_test.sh   # 1,000x random TE shuffling
```

## Architecture & Data Pipeline

```
download_and_process_files_UCSC_genes.ipynb
  → T2T RefSeq annotations + TSS coordinates

TEs_mapped_on_TSS_analysis.ipynb
  → bedtools intersect: 3.7M TEs × 10 kb TSS windows
  → Fisher exact test + 1,000x permutation for random background
  → TEs_on_genes.csv (23 MB), TEs_on_genes_counts_subfamilies.csv (86 MB)
  → enrichment_families_with_random.csv, enrichment_subfamilies_with_random.csv

Gene_ontology_analysis.ipynb + GO_subfamilies.py
  → goatools enrichment against goa_human.gaf (190 MB, not in git)
  → Top 5% genes per subfamily as foreground set
  → GO_tables/*.csv (~100 files, one per subfamily)

genes_subfamilies_network.py
  → Jaccard similarity (threshold: 0.025) between subfamily gene sets
  → NetworkX + pyvis interactive HTML
  → plots/TE_subfamilies_clustered_connectivity_network_*.svg
```

## Key Files

| File | Purpose |
|------|---------|
| `TEs_on_genes.csv` | Master TE-gene intersection table (10 kb TSS window) |
| `TEs_on_genes_counts_subfamilies.csv` | Per-gene TE counts by subfamily |
| `enrichment_families_with_random.csv` | Enrichment stats for 44 families |
| `enrichment_subfamilies_with_random.csv` | Per-subfamily enrichment stats |
| `GO_Classified_Table_Ordered_Gemini.tsv` | Manually curated GO→metagroup classification |
| `go-basic.obo` | Gene Ontology DAG (31 MB) |
| `goa_human.gaf` | Human GO annotations (190 MB, excluded from git) |

## Environment

Python 3.11. Key packages: `pandas`, `numpy`, `scipy`, `statsmodels`, `goatools`, `matplotlib`, `seaborn`, `networkx`, `pyvis`, `supervenn`, `plotly`, `adjustText`, `statannotations`. External tool: `bedtools`.

The `goa_human.gaf` file (190 MB) must be downloaded separately — it is excluded from git via `.gitignore`, along with `plots/`.

## Key Parameters to Tune

- `JACCARD_THRESHOLD` in `genes_subfamilies_network.py` — controls edge density in the network graph
- FDR threshold in `run_goatools_enrichment()` in `GO_subfamilies.py` — currently `0.1`
- TSS window size — currently 10 kb, set during the bedtools intersect step in `TEs_mapped_on_TSS_analysis.ipynb`
- Top-N percentile cutoff for foreground gene sets — currently top 5% by TE count per subfamily
