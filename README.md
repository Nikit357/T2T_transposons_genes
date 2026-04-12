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
    - TE_GO_Sankey.svg: Visualizes the flow from TE classes to functional gene groups.
    - Divergence_score_ridge_plot.svg: Shows the bimodal age distribution of TEs.
    - TE_Process_Hierarchical_Network_Refined.svg: The "Ring of Power" diagram illustrating the co-evolutionary functional network.
  ## Installation & Requirements
  The project utilizes a Python 3.11 environment.
  
  ### Key dependencies include:
  
  *Bioinformatics*: bedtools , goatools
  
  *Data Processing*: pandas, numpy, scipy
  
  *Visualization*: matplotlib, seaborn, networkx, pyvis, supervenn
  
**Author**: Daniil Nikitin

**Affiliation**: Institute of Molecular Biology, National Academy of Science of the Republic of Armenia 

## Related publications

- Nikitin D. Joint Analysis of Human Retroelements-Linked Histone Modification Profiles Reveals Quickly Evolving Molecular Processes Connected with Cancer. published online 27 Sep. 2025. https://doi.org/10.1101/2025.09.24.677146.
- Nikitin D. Transposable element–host genome evolutionary arms race revealed by multi-modal epigenomic profiling in a telomere-to-telomere human genome reference. BioRxiv 23 Mar. 2026:2026.03.19.712972. https://doi.org/10.64898/2026.03.19.712972.
- Nikitin D. Retroelements-Driven Regulatory Evolution of Human Genes and Molecular Processes: Analysis of Genome Binding Profiles of Transcription Factors and Histone Modifications. n.d. https://doi.org/10.5281/ZENODO.19052416.
- Nikitin D, Garazha A, Sorokin M et al. Retroelement-Linked Transcription Factor Binding Patterns Point to Quickly Developing Molecular Pathways in Human Evolution. Cells 2019;8(2). https://doi.org/10.3390/CELLS8020130.
- Nikitin D, Kolosov N, Murzina A et al. Retroelement-Linked H3K4me1 Histone Tags Uncover Regulatory Evolution Trends of Gene Enhancers and Feature Quickly Evolving Molecular Processes in Human Physiology. Cells 2019, Vol 8, 2019;8(10). https://doi.org/10.3390/CELLS8101219.
- Nikitin D, Penzar D, Garazha A et al. Profiling of Human Molecular Pathways Affected by Retrotransposons at the Level of Regulation by Transcription Factor Proteins. Front Immunol 2018;9(JAN). https://doi.org/10.3389/FIMMU.2018.00030.
- Nikitin D, Sorokin M, Tkachev V et al. RetroSpect, a New Method of Measuring Gene Regulatory Evolution Rates Using Co-mapping of Genomic Functional Features with Transposable Elements. Evolution, Origin of Life, Concepts and Methods 1 Oct. 2019:85–111. https://doi.org/10.1007/978-3-030-30363-1_5.

