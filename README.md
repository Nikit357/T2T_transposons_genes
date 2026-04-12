# T2T_transposons_genes
Data and code for the T2T-CHM13 analysis of TE-gene proximity. Features proximity mapping for 3.7M TEs across 28k genes, enrichment statistics for 44 families, and GO functional networks. Includes Jupyter notebooks for the IFNA cluster arms race and regulatory innovation analysis.

## Repository Structure

The repository is organized into several key components: the core processing notebooks, specialized scripts for visualization and network analysis, and a comprehensive collection of epigenomic and Gene Ontology (GO) data.
### 1. Core Analysis & Notebooks

- download_and_process_files_UCSC_genes.ipynb: Initial pipeline for acquiring T2T RefSeq gene annotations and preparing genomic intervals.
- TEs_mapped_on_TSS_analysis.ipynb: The primary analysis notebook performing proximity mapping of 3.7 million TEs to transcription start sites (TSS).
- Gene_ontology_analysis.ipynb: Performs functional enrichment analysis using goatools to link specific TE groups to biological processes.

### 2. Scripts & Frameworks
- run_permutation_test.sh: Bash script for executing the 1,000x random shuffling of TEs to establish a statistical background.
- draw_length_divergence_corr.py: Python script for generating correlations between TE evolutionary age (divergence) and sequence length.
- genes_subfamilies_network.py: Implementation of network analysis using networkx to visualize the "Ring of Power" functional co-associations.
- GO_subfamilies.py: Specialized script for processing subfamily-specific functional enrichments.

### 3. Data & Results Tables
- TEs_on_genes.csv: The master mapping file containing every TE intersection with a 10 kb TSS neighborhood.
- enrichment_families_with_random.csv: Statistical table containing Fisher exact test results and permutation-based Odds Ratios for 44 TE families.
- GO_Classified_Table_Ordered_Gemini.tsv: Curated list of GO terms manually classified into biological metagroups (e.g., Embryogenesis, Nervous System).

### 4. Specialized Directories

- GO_tables/: A library of CSV files containing specific GO term results for individual TE subfamilies (e.g., GO_terms_by_subfamilies_AluY.csv, GO_terms_by_subfamilies_L1PA2.csv).
- plots/: Comprehensive collection of SVG and PNG visualizations, including:
-   TE_GO_Sankey.svg: Visualizes the flow from TE classes to functional gene groups.
-   Divergence_score_ridge_plot.svg: Shows the bimodal age distribution of TEs.
-   TE_Process_Hierarchical_Network_Refined.svg: The "Ring of Power" diagram illustrating the co-evolutionary functional network.
  ## Installation & Requirements
  The project utilizes a Python 3.11 environment.
  
  ### Key dependencies include:
  
  Bioinformatics: bedtools , goatools
  
  Data Processing: pandas, numpy, scipy
  
  Visualization: matplotlib, seaborn, networkx, pyvis, supervenn
  
Author: Daniil Nikitin
Affiliation: Institute of Molecular Biology, National Academy of Science of the Republic of Armenia 
