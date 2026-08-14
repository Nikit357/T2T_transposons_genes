# Final scientific review — 2026-08-10

Object reviewed: `Revised_manuscript/T2T_genes_article_G3_revision_260809.docx`
against the exported figures in `current_figures_260810/` and the supplementary package in
`supplementary_260809/`.

Output: **`Revised_manuscript/T2T_genes_article_G3_revision_260810.docx`**, written by
`19_review_edits_260810.py`. 68 tracked operations, all reported, none silent.

| | |
|---|---|
| Accept All | the corrected text |
| Reject All | byte-identical to the 260809 input — verified, not assumed |
| `<w:sdt>` / `MENDELEY_CITATION` | 138 / 137 before, 138 / 137 after, contents identical |
| `<w:t>` inside `<w:del>` | 0 |
| duplicate revision ids | 0 |
| notes inserted into the docx | 0 — everything for you is in this file instead |
| P3 numbering audit | was **FAIL**, now **PASS** for all four series |

Word and LibreOffice are not installed here, so the revision display was verified structurally
rather than visually. Open the file in Word and eyeball the marks before accepting anything.

---

## 1. Numbers that did not reproduce

Every quantitative claim was recomputed from the source table. Fifteen did not reproduce; all
fifteen are corrected in the output.

| # | Claim | As written | As recomputed | Source |
|---|---|---|---|---|
| 1 | class-level summary, all TEs | observed OR **1.94**, random **1.78**, enrichment **1.097** | observed **1.95**, random **1.79**, enrichment **1.086** | see §1.1 |
| 2 | LINE count threshold of the top-5 % set | at least **10** elements | **9** | `File S1/TSS_TE_intersections` |
| 3 | TE-top set range | **27** to 42 per TSS | **26** to 42 | same |
| 4 | L2 maximum copies per TSS | **17** | **12** | same |
| 5 | hAT-Charlie maximum copies per TSS | **19** | **10** | same |
| 6 | family order in the same sentence | Alu, MIR, L1, L2, hAT-Charlie | L2 (30,376 TSS) precedes L1 (30,184), and the sentence says the families are ordered by TSS count | same |
| 7 | TSS with SVA elements | "1,140 **genes**, 2.9 %" | 1,140 **TSS**; they carry 962 distinct genes | same |
| 8 | hAT-hAT19 outlier | divergence "**1.39 higher**" | **1.39-fold** (314.0 against a family mean of 226.2) | `File S1/enrichment_families` |
| 9 | least-present functional group, Figure 6B | protein ubiquitination alone, **1** group at 2 terms | protein ubiquitination **and embryogenesis**, both at 2 | recomputed Figure 6B |
| 10 | all TEs, lowest divergence | rRNA binding, **spermatogenesis**, **mitotic spindle localization**, subcortical maternal complex, **fatty acids catabolism** | ribonucleoprotein complex, rRNA binding, subcortical maternal complex, establishment of the Sertoli cell barrier | see §1.2 |
| 11 | LINEs, highest divergence | voltage-gated K channel complex, **cell adhesion**, differentiation, **proliferation**, **VLDL components** | potassium ion transport, voltage-gated K channel complex, protein kinase binding, ventricular cardiac muscle cell differentiation, ECM structural constituent | see §1.2 |
| 12 | olfactory sharing | "**all TEs** and LTRs of high divergence sharing olfactory receptors" | all TEs of highest divergence carry **no** olfactory or smell term at FDR 0.05 **or** 0.1; only LTRs do | `GO_classes_divergence_fdr005.csv` |
| 13 | three most significant divergence terms | nervous system development, postsynaptic membrane, **"synaptic chemical transition"**, FDR 10-16 – 10-11 | postsynaptic membrane (3.0 × 10-16), chemical synaptic transmission (2.8 × 10-12), axon guidance (1.4 × 10-10) | see §1.3 |
| 14 | SINEs of lowest divergence | "**no** significant GO enrichments" | two significant terms: protein binding and extracellular exosome, FDR 8.2 × 10-4 | same |
| 15 | raw-odds-ratio reversal | "several groups" | **four of the six classes** (LINE, LTR, DNA, Helitrons are above 1 on the raw OR and below 1 on the ratio) | Table 1 |

### 1.1 The all-TE enrichment sentence quoted the wrong row

`In general, the observed OR for all TEs was 1.94 (1.78 the random one), which was showing an
enrichment by a factor of 1.097`.

The triple **1.9509 / 1.7782 / 1.0971** is the **ERVK family** row of
`enrichment_families_with_random.csv`. It is not the all-TE figure.

For all TEs: 582,540 of 3,709,429 elements fall in a 10 kb window, against a mean of 542,973.1
over the 500-seed background (`output/permutation_counts_10kb/`). The enrichment score is a ratio
of two count ratios and therefore does not depend on the window-geometry constant at all:

```
(582540 / 3126889) / (542973.1 / 3166455.9) = 1.0865
```

so **1.086**, and the two odds ratios round to **1.95** and **1.79**. The three numbers as written
were also mutually impossible: 1.94 / 1.78 = 1.090, which cannot be 1.097 under any rounding.

*This is the one correction that rests on a reconstruction of the odds-ratio denominator rather
than on a stored table — the reconstruction reproduces all six published class odds ratios to
three decimals, and the 1.086 is independent of it, but the two odds ratios are worth confirming
against the frozen notebook's own all-TE row if you can find it.*

### 1.2 Five terms in the divergence subsection are not significant at FDR 0.05

The paper states one threshold, FDR 0.05, with no suggestive band. These five sat in the rejected
0.05 < FDR ≤ 0.1 band and were reported as results:

| Term named | Group | Actual FDR |
|---|---|---|
| establishment of spindle localization ("mitotic spindle localization") | all TEs / lowest | 0.071 |
| fatty acid catabolic process ("fatty acids catabolism") | all TEs / lowest | 0.080 |
| cell-cell adhesion ("cell adhesion") | LINE / highest | 0.066 |
| cell population proliferation ("proliferation") | LINE / highest | 0.052 |
| very-low-density lipoprotein particle ("VLDL components") | LINE / highest | 0.056 |

"spermatogenesis" does not exist as a term at this level at either threshold; the term the sentence
was probably reaching for is **establishment of the Sertoli cell barrier** (FDR = 0.042), which is
what the corrected sentence names.

### 1.3 "synaptic chemical transition" is not a GO term

The real term is **chemical synaptic transmission**. Two further problems in the same sentence:
*nervous system development* has 642 annotated genes, so the 500-gene cap that the paragraph itself
declares keeps it out of Figure 5A; and under that cap the third most significant shared term is
*axon guidance*, not the one named. The FDR range therefore closes at 10-10, not 10-11.

### 1.4 Everything else reproduced

For the record, these all check out exactly and were left alone: the entire interferon-alpha
subsection (175 elements, 77 L1, 36 subfamilies, 351 vs 181.4 per Mb, 1.94-fold, 135.7 vs 188.2,
all four tests with their z-scores and p-values, the leave-one-out means, the 565,459 total against
545,659 clade-assignable L1 copies — which are two different quantities and both right); the whole
assembly bound (208.8 Mb, 6.70 %, 0.41 / 0.49 / 0.55 %, 4.7 % vs 82.9 %, 0 bp in the domain); all
eight rows of Table 3 and the three MGI mapping rates; Tables 1 and 2 in full; the 7 enriched and 9
depleted families with all 16 fold changes; every one of the ~40 GO FDR values quoted at family and
class level; the 10-term MIR band and the 3-term LTR band; 13 families with terms and 31 without
with their exact class breakdown; all seven Figure 7 statistics; the window concordance (0.828–1.000,
largest p 0.0090, Bland-Altman 0.027–0.084, 11 of 50 flips with the exact list); the gene-set
stability (0.285–0.931, 3.5 × 10-205, Kendall 0.482–0.793, and both "weakest cell" attributions);
the percentile results (90.3 % median, 144 lost, 902 gained, 8 of 9); the whole GO grid (0.440,
0.677, 0.614, 0.022, 3 of 9); the functional-group counts 22 of 25, 21 and 24 and every starred
cell in Figures 4B, 5B, 6B and S9B; the Helitron K-S p of 0.222; the 21 and 13 family
divergence/length comparisons with all nine percentages; and the SSU72L cluster's 3 SVA_B copies.

---

## 2. Structural defects only visible in the Accept-All view

| # | Defect | Action taken |
|---|---|---|
| S1 | **A second Table 1 and a second Table 2**, after the Funding section. They drop the unadjusted Fisher p and merge mean ± SD into one uncopyable cell — the two things you fixed in the Results pair. | Removed. Both were entirely our own tracked insertions, so Reject All is unaffected. The Results pair, placed at first citation, is now the only one. |
| S2 | **An empty 7 × 11 table** between Table 1 and the Table 2 caption. Its text was already struck through, but a table row is only removed by Accept All when its `w:trPr` carries a `w:del`, so the accepted document kept a blank grid. | The seven rows are now marked deleted, so Accept All removes them. |
| S3 | **The Ethical statement stopped mid-clause** — "…obtained from the T2T genome assembly and the RepeatMasker track" with no period. The sentences it needed were sitting in the supplementary block instead, four sections earlier. | Section completed; the misplaced copy removed. |
| S4 | **The Discussion roadmap listed the subsections in an order the Discussion does not follow** (cancer and null model before the mechanistic framework; the document has them after). | Reordered to match, and the "makes the positional pattern interpretable" clause dropped. |
| S5 | **`the same edge and term-size filters as Figure S8`** in the legends of Figures S7 and S10. The renumbering moved the first full network to S6; S8 is now the log-scaled TSS distributions, which have no filters. | Both repointed to Figure S6. |
| S6 | **The SVA class-level association was pointed at Figure 5A**, the divergence-stratified panel. It is a count-based result and belongs to Figure 4A. | Repointed. |
| S7 | **Figures S6 and S7 were cited only from the legends of Figures 4 and 5**, never from prose. | One prose citation added to each, at the position they already occupied. P3 now passes. |
| S8 | **The five workbooks were described twice**, in Supplementary material and again in Data availability. | The second copy now refers to the first. |

---

## 3. What you have to do in Word — I did not touch any of it

### 3.1 The 28 orphan citation paragraphs, and what deleting them does to the reference list

Twenty-eight paragraphs consist of nothing but citation markers, left behind when the extended
discussion was moved out. They are at these places in the Accept-All text:

- one block of 2 after "Hypotheses for future testing" — `(9)(10)(11)(8)(5,6,18,42,45,46)` and `(47)(48)`
- eleven between "Mechanistic framework" and the Figure 9 legend
- fourteen between the Figure 10 legend and "Connection of TE enrichment with cancer"
- one after the File S6 sentence — `(115,116)(117–119)(120,121)(122–124)(125)(126)`

**They carry 83 distinct references, and 68 of those are cited nowhere else in the manuscript.**
Deleting the paragraphs and refreshing Mendeley therefore takes the reference list from 126 entries
to about 58, and **renumbers every citation in the paper**. That is the correct outcome — those 68
references belong to the extended discussion, and `Extensive_discussion_260803.docx` carries them in
its own bibliography (78 Mendeley controls, 117 entries), so nothing is lost. But do it in one pass:
delete all 28, refresh Mendeley, then re-read the Results and Discussion for citation numbers that
moved.

### 3.2 Citations that are missing and that only you can add

Five passages describe the literature with no citation attached. Each needs one inserted in Mendeley:

| Location | Uncited claim |
|---|---|
| Principal findings | "earlier estimates that roughly a quarter of human promoters are TE-derived" — also a rounded quantity; a citation would let you state the published figure exactly |
| Comparison with prior work, first sentence | "The previous repeat-based categorization of genes used 20 kb and 2 kb windows…" — this is Lu et al., cited as (43) in Results but nowhere in this paragraph |
| same subsection | "the SVA association with premature transcription termination that was **predicted from the internal polyadenylation signal** of these elements" |
| same subsection | "An early survey that found Alu elements preferentially near metabolism, transport and signaling genes **examined only chromosomes 21 and 22**" |
| same subsection | "A proximity study that placed SVA elements near zinc finger clusters **used hg19 and 1 Mb bins**" |

### 3.3 Mendeley formatting

Five in-text citations render as author-date instead of numbers, so Mendeley has not been refreshed
since they were inserted:

| Citation | Section |
|---|---|
| `(Nikitin, 2026)` | Introduction, second paragraph |
| `(Klopfenstein et al. 2018)` | Methods, Gene Ontology |
| `(GitHub - gecko984/supervenn, n.d.)` | Methods, Visualization |
| `(Hagberg et al. 2008)` | Methods, Visualization |
| `(Pandas team, n.d.)` | Methods, Data analysis |

Refresh in both this document and `Extensive_discussion_260803.docx`.

### 3.4 `[ZENODO DOI]`

Still a placeholder in Data availability. It is the only piece of scaffolding left, and it has to
stay until you deposit.

---

## 4. The exported figures — three items are Figma-side

Everything in `Figures_renaming_260809.md` §1 and §2 checks out in `current_figures_260810/`. I read
the panel letters and text out of all 23 PDFs:

- **New S1 has five panels A–E**, convergence first, length distributions as D and E ✓
- **New S12 is a single panel with no letter** ✓
- The S6 ⇄ S8 and S7 ⇄ S9 swaps landed correctly: S6 is the class-count network, S7 the
  divergence-stratified one, S8 the log-scaled TSS-by-family histograms, S9 the three family panels ✓
- S10 full family network, S11 four concordance panels, S13 three GO-grid panels ✓
- §3 item 1 fixed: **Figure 4 now carries its panel letter A** ✓
- §3 item 3 fixed: **Figures 4, 5, 6 and S7 all carry the `-log10(FDR)` colourbar** their captions
  promise ✓

Three things are still open, and all three are artwork:

| | Item |
|---|---|
| **F1** | **Figure 8 has no panel letter C.** A, B and D are all present at x ≈ 14 in the left margin; the subfamily-composition panel in the right column has none. The manuscript cites Figure 8C twice — in the Results and in the Figure 8 legend — so this one is load-bearing. |
| **F2** | **Figure 3's panel letter A overlaps the y-axis tick.** It extracts as the single token `A1.0`, i.e. the letter is flush against the `1.0` label. B, C and D are clear. |
| **F3** | Figure 8C prints `Fisher p = 3.2*10-6` with an asterisk where the manuscript uses `3.2 × 10-6`. Cosmetic, but it is the same statistic in two places. |

---

## 5. Things I deliberately left alone

- **The abstract is 232 words.** `nar_review_tools` flags it against NAR's 200-word cap; G3 allows
  **under 250** (`G3_article_guidelines.md` §4), so it is compliant. No action.
- **`P-value` with a capital P** in Data availability. It names the literal column header of the
  `GO_tables/` CSVs, so the capital is correct there.
- **`(FDR corrected p-value 10-16 – 10-10)`** is still an order-of-magnitude range rather than three
  printed values. The three exact FDRs are in `File S3, sheet GO_by_divergence`, which the same
  paragraph cites, so the claim is traceable. Say the word if you want the values inline.
- **"Ring of Power"** in the Discussion. It is your figure name from the companion paper and reads as
  deliberate; I am flagging it rather than removing it.
- **Divergence quoted as a percentage** in the gene-set paragraph (9.8 %, 27.9 % …) while the rest of
  the paper uses the score (135.7, 188.2). Both are right — 9.8 % is a score of 98 — and all eight
  values reproduce exactly, so this is a terminology choice, not an error. Figure S5B's axis decides
  which is better.
- **Reference numbering.** Mendeley owns it. §3.1 will move most of it.

---

## 6. Files

| File | |
|---|---|
| `Revised_manuscript/T2T_genes_article_G3_revision_260809.docx` | untouched; now a read-only input of record |
| `Revised_manuscript/T2T_genes_article_G3_revision_260810.docx` | **the file to work with** |
| `19_review_edits_260810.py` | the edit script; a deterministic rebuild from its input, so re-runnable |
| `review_report_260810.md` | this file |
