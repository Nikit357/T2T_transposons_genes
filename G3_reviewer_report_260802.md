July 29, 2026
RE: G3-2026-406828


Dear Dr. Nikitin:

We are pleased to conditionally accept your manuscript titled "Evolutionary arms race between transposable elements and human genes: telomere-to-telomere genome comprehensive analysis identifies young L1 clusters in the interferon-alpha domain" for publication in G3: Genes|Genomes|Genetics. We expect all remaining minor revisions can be completed within 30 days. If you require additional time, please let us know. In addition to the reviewer comments, we also note that the GitHub page link provided in the data availability statement is broken; please ensure that the data availability statement is updated in your revision.

Follow this link to submit the revised manuscript: https://g3.msubmit.net/cgi-bin/main.plex?el=A7NQ1IFP2A1gOY1I4A9ftdolTB6k66m5lMLcBNIDNmgZ

Prior to uploading your revised manuscript please format it according to G3 style and ensure you have all required elements. Author Guidelines are at https://academic.oup.com/g3journal/pages/author-guidelines#section-11. These guidelines have been updated with new requirements and your careful attention is required to avoid delays.

In your final submission, please include:
1. A clean version of your manuscript, formatted for G3
2. A highlighted or tracked version of your manuscript that links your response to reviewers via the current text
3. A separate document with a response to each of the editor's/reviewers' comments

Thank you for submitting your research to G3. As a fully open access journal of the Genetics Society of America (GSA), our mission is to publish peer-reviewed and peer-edited reproducible science with high-quality data. Thank you for your contribution.

Sincerely,

Arun Sethuraman
Associate Editor
G3: Genes
Genomes
Genetics



Andrew Kern
Senior Editor
G3: Genes
Genomes
Genetics



----------------------------------------------------------------------------

Reviewer comments:
Reviewer #1 :

In this manuscript, Daniil Nikitin leveraged the complete T2T human genome assembly to systematically map transposable elements (TEs) within 10 kb of gene transcription start sites, quantified their enrichment or depletion relative to random expectations, and performed GO term analyses to identify functional associations between specific TE classes/families and nearby genes. The work is data-rich and identifies several intriguing associations, most notably a cluster of young L1 elements in the interferon-alpha gene region. Overall, while this study provides a potentially valuable resource for the field, there are several major issues that need to be addressed.

Major concerns:
1. The title and Abstract frame the study as evidence of an ongoing "evolutionary arms race" and suggest that young L1 elements "influence" innate immune responses. However, the data are purely correlative and static-they show co-localization, not causation. The "arms race" narrative implies a dynamic process that this study does not directly measure. Please revise throughout the manuscript to reflect the descriptive, hypothesis-generating nature of the findings.

2. The authors provide the average L1 divergence (95-161.7) for the interferon-alpha region, but do not compare this to genome-wide background or to matched random regions. Please perform a statistical test (e.g., permutation test) to demonstrate that the low divergence in this region is significantly different from random expectation. Also, please report the number of L1 elements in this region and their specific subfamilies to confirm that the observation is not driven by a few outliers.

3. The study relies solely on a single genome assembly (T2T-CHM13) and does not integrate orthogonal data (e.g., epigenetic marks, expression data, eQTL) to support functional relevance. It would be helpful to add such analyses using public resources (e.g., ENCODE, GTEx, TCGA).

4. The author notes methodological differences between this study and the prior landmark study by Lu et al. (14), but does not directly compare results on the same dataset. Is this difference due to methodology or the updated genome assembly?

5. Although the manuscript notes the choice of a 10 kb window and the top/bottom 5% cutoffs, it would be better to perform sensitivity analyses for key findings (e.g., the IFNA-L1 association, the SVA-termination association) using alternative window sizes (5 kb and 20 kb) and alternative percentiles (e.g., 10%).

6. The Discussion is disproportionately long and reads more like a review article on TE evolution than a focused discussion of the present results. It is suggested that the author substantially shorten and refocus the Discussion on: (1) what the author found; (2) how it compares to prior work; (3) limitations; and (4) specific hypotheses for future testing.

Minor concerns:
1. The number of random permutations is stated as 500 in the Results and Discussion but as 1000 in the Methods. Please unify this and justify the chosen number.

2. The GO analysis uses an FDR threshold of 0.1. Please justify this relatively lenient threshold or, ideally, tighten it to 0.05.

3. For all figures reporting statistical significance, please clearly indicate in axis labels or legends whether p-values are raw or FDR-adjusted, and explicitly state "adjusted p-value" if applicable. In Figure 1D, for example, the third vertical bar plot labeled "-log10(p-value)" should read "-log10(FDR-corrected p-value)".

4. Table 1 appears too wide, with the left margin extending beyond the row names. Please reformat it for better readability-consider splitting into two tables, reducing column widths, or moving less critical columns to supplementary materials.

5. Line 256: "flavone metabolism" - please verify the term (did you mean flavonoid?).

6. Figures 4A, 5A, and 6A are dense network visualizations that are difficult to interpret in print. Please provide simplified versions in the main text (e.g., filtering by edge weight or minimum shared genes) and move the full versions to supplementary materials.


Associate Editor comments:

This is an interesting and important contribution towards the human genomics literature; a complete annotation of TEs is imperative towards understanding the complex landscape of functional variation in transposable elements. We appreciate the diligence in this work performed by the author.

Please ensure that code/data links are made available in the revision.

Additionally, it would be very helpful to have a browser/track instance for the TE annotations made available with the GitHub repository. We hope that this would be a straightforward addition to the manuscript.