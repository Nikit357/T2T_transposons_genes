# Daniil's review of Claude's manuscript output — worked examples and editing principles

**Date:** 2026-08-09
**Manuscript:** G3-2026-406828, *"Telomere-to-telomere co-mapping of transposable elements and human
genes identifies a cluster of young L1 elements in the interferon-alpha domain"*
**Compared documents:**

| Role | File |
|---|---|
| Claude's output | `Revised_manuscript/T2T_genes_article_G3_revision_260807.docx` |
| Daniil's manual edit of it | `Revised_manuscript/T2T_genes_article_G3_revision_260807_manual.docx` |

**Purpose of this document.** It is a training reference for future Claude sessions (and for a
dedicated `scientific-review` skill) working on Daniil Nikitin's manuscripts. Part I catalogues
every manual edit Daniil made to Claude's output. Part II derives the editing principles from that
catalogue, in a form that can be applied prospectively rather than only recognised after the fact.
Part III gives the verification protocol. Part IV lists the defects that this comparison exposed and
that are still open. Part V is the pre-hand-off checklist.

The governing observation: **almost none of Daniil's edits were about taste.** Of 134 word-level
edit sites, roughly 100 were substantive, and each one removed a statement that could not be
checked, redirected a number to the object that actually contains it, renumbered something so a
reader could find it, or deleted a sentence about the revision process that a reader of the
published paper has no way to interpret. Treat every item in Part I as a defect class, not a
preference.

---

## 0. How the comparison was made (reproduce this before any future review pass)

Both `.docx` files contain tracked changes, so a naive text dump compares the wrong things. Extract
the **final** text — all insertions accepted, all deletions removed — from each, then diff. The four
scripts used are kept in `revision_G3/review_tools/`.

```bash
# venv with python-docx + lxml
cd revision_G3/review_tools
V=~/venvs/collagen_3_11/bin/python
M=../Revised_manuscript
$V extract_docx.py "$M/T2T_genes_article_G3_revision_260807.docx"        a.txt final  # accept-all
$V extract_docx.py "$M/T2T_genes_article_G3_revision_260807_manual.docx" b.txt final
$V extract_docx.py "$M/T2T_genes_article_G3_revision_260807.docx"        a0.txt orig  # reject-all
$V diff_para.py  a.txt b.txt > diff_sentences.txt   # what argument changed
$V diff_words.py a.txt b.txt > diff_words.txt       # spelling, citations, single clauses
$V refs_order.py a.txt b.txt                        # sequential-numbering audit
```

`extract_docx.py` walks `word/document.xml` with lxml, emits one line per `<w:p>`, and for the
*final* view skips every `<w:t>` that has a `<w:del>` ancestor. Then:

* **sentence-level diff** — split each paragraph on `(?<=[.:;])\s+(?=[A-Z0-9(])`, run
  `difflib.SequenceMatcher` over the sentence list. Good for seeing what argument changed.
* **word-level diff** — tokenise on whitespace, same matcher, print 8 tokens of context either
  side. Good for seeing spelling, citation numbers and single-clause deletions.
* **citation-order scan** — regex every `Figure N`, `Figure SN`, `File SN`, `Table N` token and
  print the order of first appearance. This is the only reliable way to audit sequential numbering.

Measured state of the two documents:

| | 260807 (Claude) | 260807_manual (Daniil) |
|---|---|---|
| Non-empty paragraphs | 384 | 500 |
| Words | 18,323 | 18,318 |
| `<w:ins>` | 398 | 302 |
| `<w:del>` | 570 | 487 |
| `<w:sdt>` (content controls) | 129 | 137 |
| Mendeley citation controls | 128 | 136 |

Word count is nearly identical, paragraph count is up by 116: Daniil deleted prose and inserted two
full Word tables plus a section. **A flat word count is not evidence that little changed.**

Edit sites by type:

| Type | Sites |
|---|---|
| Substantive text edits | 98 |
| Citation renumbering | 22 |
| Superscript exponent → plain text | 8 |
| British → American spelling | 6 |
| **Total** | **134** |

---

# Part I — Catalogue of the manual edits

## I.1 Abstract

**Edit A1 — replaced a methods sentence with a scope sentence.**

> **Claude:** The complete telomere-to-telomere (T2T) assembly makes it possible to measure where TEs
> sit relative to every annotated gene.
> **Daniil:** The complete telomere-to-telomere (T2T) assembly enables comprehensive investigation of
> TE contributions to gene regulation.

Claude described *what the analysis does*. Daniil stated *what the assembly makes possible*, which is
the claim an abstract's second sentence is for. "measure where TEs sit relative to" is also
conversational — see P4.

**Edit A2 — moved the statistical battery out of the Abstract, added the age qualifier.**

> **Claude:** 77 of its 175 TEs are L1 copies from 36 subfamilies, at mean divergence 135.7 against
> 188.2 genome-wide, a deficit that survives 10,000 random windows matched for L1 count (empirical
> p = 0.006) or gene density (p = 0.0017).
> **Daniil:** 77 of its 175 TEs are **evolutionary young** L1 copies from 36 subfamilies, at mean
> divergence 135.7 against 188.2 genome-wide.

The two p-values are correct (`output/ifna_test_results.csv`: T2 = 0.006099, T3 = 0.001674) and are
retained in Results. The Abstract keeps the two numbers a reader can hold — 77/175 and 135.7 vs
188.2 — and drops the null-model description, which cannot be evaluated in an abstract anyway.

**Edit A3 — replaced a negative hedge with a forward-looking hypothesis.**

> **Claude:** This is consistent with recent L1 activity; whether these elements affect innate immune
> gene regulation remains untested.
> **Daniil:** This is consistent with recent L1 activity and could indicate a recent evolutionary
> arms-race affecting innate immune gene regulation.

Note carefully, because this looks like a violation of P1 and is not. "remains untested" is an
assertion about the *literature*, and Claude had no citation survey to support it. "could indicate"
is explicitly modal and is followed, at the end of the Introduction, by a sentence that states the
design cannot detect an arms race (Edit I3). Daniil replaces an unsupported negative claim with a
labelled hypothesis plus a stated limitation. See **P1c**.

**Edit A4.** "candidate loci for follow-up" → "candidate loci for **further investigation of
evolutionary arms race**". Names the follow-up rather than gesturing at one.

## I.2 Introduction

**Edit I1 — deleted a mid-paragraph defensive disclaimer.**

> **Deleted:** Throughout this paper the term "arms race" is used only when reporting that cited
> literature: the design applied here is correlative and static, and cannot establish an arms race
> for our own observations.

Claude inserted this into the *opening paragraph*, between the definition of the arms race and the
mechanism sentence that follows it. The same statement already exists in Limitations ("the causal
vocabulary of an evolutionary arms race is used in this paper only when reporting other people's
conclusions"), where it survives untouched. **The caveat was not weakened — it was de-duplicated
and moved to the section that exists for caveats.** See **P8**.

**Edit I2 — precision of a technical noun.**

> **Claude:** how close each TE group **sits** to each **class of gene**
> **Daniil:** how close each TE group **locates** to each **functional gene group**

"class of gene" collides with "TE class", the load-bearing term of the paper. "functional gene
group" is the term used in Results and in Figure 4B. Terminology is single-valued throughout.

**Edit I3 — replaced Claude's editor note with the scope sentence.**

The end of the Introduction carried:

> `[EDITOR NOTE] G7: Materials and methods is moved ahead of Results, which is G3's printed order.
> This is a structural relocation, not a Word tracked move, so Reject All restores the original text
> but leaves the sections in this order.`

Replaced with:

> Albeit the current analysis is purely correlative and cannot detect the ongoing evolutionary arms
> race, it points at genome loci and molecular processes that could be impacted by it for future
> investigation.

Two things at once: process scaffolding out (**P7**), and the arms-race caveat lands where the
Introduction states scope, one sentence, once (**P8**).

## I.3 Materials and methods

**Edit M1 — every citation renumbered by order of first appearance.** 22 sites. The Methods
citations went from the 98–123 range down to 12–41; Discussion citations went from 20–31 up to
49–60 and from 84–85 to 113–114. Claude had appended new references at the end of the bibliography;
Daniil renumbered the whole list so that citation *n* is the *n*-th cited work. See **P3**.

**Edit M2 — added citations to a literature-survey passage that had none.** The window-size
justification listed five prior studies in prose with no citations attached:

> **Claude:** A study of interferon-inducible enhancers spread by LTR elements used a 10 kb window;
> an analysis of TE content around duplicated and singleton genes used 4 kb and 20 kb windows; …
> **Daniil:** … used a 10 kb window **(9)**; … used 4 kb and 20 kb windows **(10)**; … 5 kb either
> side of coding exons **(11)**; … 1 kb downstream) **(8)**; and our own … applied here **(3,18)**.

Five claims about five specific papers, zero citations. This is the most direct possible violation
of **P2**: a reader cannot check any of it.

**Edit M3 — deleted a subsection heading and folded the content into the paragraph it belongs to.**

> **Deleted:** all key results are additionally reported at 5 kb and 20 kb (below). **Sensitivity to
> the window and to the gene-set cut-off.** The full enrichment analysis was repeated…
> **Daniil:** Because the **10kb window** choice is a convention rather than a measurement, the full
> enrichment analysis was repeated for TSS neighborhoods of 5 kb…

One forward reference removed ("below"), one redundant pseudo-heading removed, one sentence doing
the work of three. Also: `the choice` → `the 10kb window choice` — the antecedent was two sentences
back.

**Edit M4 — deleted a redundant clause about the 10 kb background.**

> **Deleted:** , with 500 fresh permutations per window; the 10 kb background is the one described
> above.

Stated in the next subsection. See **P6**.

**Edit M5 — bounded a gene-set size.** "(2,872 genes)" → "(2,872 genes **at maximum**)", same for
"(1,436 genes at maximum)". The top-5 % set is a *ceiling*: ties at the boundary mean the realised
set is smaller for some TE groups. Claude stated the nominal number as if it were the realised one.

**Edit M6 — named what the permutation control is for, in the sentence that introduces it.**

> **Claude:** To control TE enrichment near genes against a random background, we performed 500
> random permutations…
> **Daniil:** To control TE enrichment near genes against a random background **and remove a
> potential length bias**, we performed…

**Edit M7 — the largest single deletion: the length-bias justification paragraph.**

> **Deleted:** The permutation background is what makes these enrichment values interpretable,
> because it removes a length bias that the Fisher exact test cannot. The probability that an element
> intersects a fixed window grows with the length of the element, so the random odds ratio scales
> almost linearly with mean element length across the 44 families (Pearson R = 0.985, n = 44, mean
> lengths 122-6,357 bp): Alu elements average 316 bp and a random OR of 1.54, whereas L1 elements
> average 6,357 bp and a random OR of 2.66. Reporting the observed odds ratio alone would **therefore**
> systematically under-call short elements and over-call long ones, **and** every enrichment statement
> in this work is consequently made on **the ratio of the observed to the random odds ratio rather
> than on the observed odds ratio itself**.
> **Daniil:** **Since** reporting the observed odds ratio alone would systematically under-call short
> elements and over-call long ones, every enrichment statement in this work is consequently made on
> **this enrichment score**.

Six sentences to one. What was cut and why:

* *"is what makes these enrichment values interpretable"* — an editorial self-assessment. **P1.**
* the R = 0.985 correlation — it is a real measured result (`output/length_bias_correlation.json`),
  but it is a *result*, and Methods is not where results are reported. It has no figure or table in
  this manuscript, so under **P2** it cannot stand in the text as written.
* *"the ratio of the observed to the random odds ratio rather than on the observed odds ratio
  itself"* — the definition was given three sentences earlier as "enrichment score". Saying it twice
  in two vocabularies invites the reader to wonder whether they are two things. **P4.**
* The clause was also reversed into a `Since …, …` construction, which puts the reason before the
  consequence and removes both "therefore" and "and".

**Edit M8 — rewrote the convergence justification into figure-readable numbers.**

> **Claude:** Five hundred permutations **are sufficient for that purpose, and this was verified
> rather than assumed.** Recomputing the random odds ratio from truncated **prefixes** of the
> permutation set shows that by N = **250** the running mean **is already within 0.06 standard
> deviations** of its N = 500 value for the worst-behaved class and within 0.10 standard deviations
> for the worst-behaved family, with the standard deviation itself estimated to within 4 %; at
> N = 100 the drift is still 0.18 standard deviations. The convergence trajectories are shown in
> **Figure S13B-D** and the checkpoint values are provided in the repository.
> **Daniil:** Five hundred permutations **were considered sufficient as** recomputing the random odds
> ratio from truncated **parts** of the permutation set shows that by N = **140** the running mean
> **differs from the N = 500 value by less than 1 %** (**Figure S1A**), and standard deviation at
> N = 250 differs by less than 6 % of its N = 500 value (**Figure S1B**). At the level of families,
> after N = 315 standard deviation of the least stable family was below 10 % of the final N = 500 one
> (**Figure S1C**). The checkpoint values for the random permutation curves are provided in the
> repository.

Four separate principles fire here:

1. *"and this was verified rather than assumed"* — deleted. A claim about the authors' diligence,
   not about the data. **P1.**
2. Units changed from **standard deviations of drift** to **per cent of the N = 500 value**. The
   figure's y-axis is a fraction of the N = 500 value with a ±1 % shaded band. Claude's numbers were
   in units the reader cannot read off the panel. **P2.**
3. N values changed from the coarse checkpoint grid (50/100/200/250/300/400) to the actual crossing
   points (140, 250, 315). Claude quoted the nearest grid row; Daniil quoted where the curve
   actually crosses the threshold. **P5.**
4. Every clause now carries its own panel pointer — S1A, S1B, S1C — instead of one range
   "Figure S13B-D" for three different claims. **P2/P3.**

**Edit M9 — renumbered this supplementary figure from S13 to S1** because after the Methods/Results
reorder it is the first supplementary figure cited. The length-distribution panels that used to be
S1A/S1B became S1D/S1E. **P3.** (See Part IV — the legends have not been updated to match.)

**Edit M10 — moved the "Statistical tests" subsection after the interferon-alpha permutation test.**
Claude had it between the permutation background and the IFNA test, splitting one subject across
two non-adjacent passages. Daniil placed the two permutation-based methods together and the generic
statistics block after them. Also added "The 220 kb domain … was **separately** tested" so the
reader knows this is a different null model from the genome-wide shuffle. **P4.**

**Edit M11 — moved the concordance-methods sentence into "Statistical tests", where it belongs, and
cited every method.**

> **Claude (in the sensitivity subsection):** Agreement between conditions was measured by Spearman
> and Pearson correlation of the observed to random odds ratio, by Bland-Altman comparison, by the
> overlap coefficient and a hypergeometric test for the gene sets, by Kendall tau with a bootstrap
> confidence interval for the gene rankings, and by a label-shuffling permutation test of the
> observed correlation, **so that "concordant" is a measured rather than a visual statement.**
> **Daniil (in Statistical tests):** Agreement **between window sizes (5, 10 and 20 kb) and gene
> percentiles (top 5 and 10 %)** was measured by Spearman and Pearson correlation **(24)** of the
> observed to random odds ratio, by Bland-Altman comparison **(26)**, by the overlap coefficient and
> a hypergeometric test for the gene sets **(27)**, by Kendall tau **(28)** with a bootstrap
> confidence interval for the gene rankings, and by a label-shuffling permutation test of the
> observed correlation.

Five statistical methods, zero citations, closing on a self-congratulatory clause. All five now
cited; the clause deleted; "between conditions" replaced by the actual conditions. **P1 + P2.**

**Edit M12 — p-value labelling scope narrowed to what is true.**

> **Claude:** p-values that are raw are labelled as such in the **figures, tables and captions**.
> **Daniil:** … labelled as such in the **figure captions and tables**.

The panels themselves do not carry raw/adjusted labels; the captions do. Claude asserted a property
of the artwork that the artwork does not have. **P1/P2.**

**Edit M13 — pointed the supplementary CSV names at a resolvable location.** Twice, Daniil added
the repository URL next to the filenames (`enrichment_families_with_random.csv` etc.), because
those filenames are not in the supplementary workbooks — they are in the repository. Also added a
standalone sentence after the Human TEs subsection: *"These files are deposited in the GitHub
repository accessible via the link https://github.com/Nikit357/T2T_transposons_genes."* **P2.**

**Edit M14 — AI usage subsection expanded and made specific.**

> **Claude:** Gemini PRO (122) was used for code refining… Chat GPT was used for grammar corrections
> of the manuscript (123).
> **Daniil:** Gemini **2.5** PRO (40) **and Claude Code Opus 5.0 (ref)** was used for code refining…
> Chat GPT was used for grammar corrections of the manuscript (41). **All code writing was done in a
> plan-then-act mode: the AI agent first wrote a research and implementation plan, the author then
> reviewed and commented on it, and after two to five review cycles the agent implemented the
> requested feature. AI specifications can be found in the project repository
> (https://github.com/Nikit357/T2T_transposons_genes) as CLAUDE.md files in all the subfolders. The
> intermediate research outputs and implementation plans, as well as additional documentation are
> stored in other markdown files.**

Version numbers pinned, the human-in-the-loop protocol stated with a number ("two to five review
cycles"), and the prompt/specification artefacts pointed at a public location. The same
traceability standard is applied to the AI method as to a statistical method. **P2.**

**Edit M15 — heading and terminology fixes.** "Mapping of TEs on gene TSS **10 kb** neighborhoods"
→ "…gene TSS neighborhoods" (the subsection now covers 5/10/20 kb). "Data and code availability**.**"
promoted from a run-in bold lead to a real subsection heading.

## I.4 Results

**Edit R1 — the editor note was replaced by the actual tables, placed at first citation.**

Claude left:

> `[EDITOR NOTE] Tables 1 and 2 replace the original 11-column Table 1 (reviewer minor comment 4).
> The values are unchanged; the unadjusted Fisher p-values moved to File S1, sheet
> enrichment_classes. Apply the journal's table style to both when formatting.`

Daniil replaced it with two full Word tables, immediately after the sentence that first cites them,
and **restored the raw Fisher p-value column that Claude had exiled to the supplementary workbook**:

| Table 1 — Claude's columns | Table 1 — Daniil's columns |
|---|---|
| Class | Class name |
| TEs in TSS windows | TE count in TSS |
| TEs total | TE count total |
| Observed OR | Odds ratio |
| Observed / random OR | Mean of random odds ratios |
| — | SD of random odds ratios |
| — | Observed to random odds ratio fold change |

| Table 2 — Claude's columns | Table 2 — Daniil's columns |
|---|---|
| Class | Class name |
| Adjusted Fisher p | **Fisher p-value** |
| Random OR (mean ± SD) | Adjusted Fisher p-value |
| Empirical p | Empirical p-value |
| Adjusted empirical p | Adjusted empirical p-value |

Every value is identical to `output/Table1.csv` / `output/Table2.csv`; only the presentation changed.
The reasoning is **P2**: a reader holding the printed table should not have to download an Excel
workbook to see the unadjusted p-value that the adjusted one was derived from, and "mean ± SD"
collapsed into one cell cannot be sorted, copied or recomputed. **A supplementary file is a place to
put more detail, never a place to put a number the main text depends on.**

**Edit R2 — precision of an enrichment statement.** "the significance of genes TSS neighborhoods by
TEs" → "the significance of genes TSS neighborhoods **enrichment** by TEs". The noun that carries
the test was missing.

**Edit R3 — supplementary panel pointers corrected for the renumbering.** "(Figure S1A, 1B)" →
"(Figure S1D, 1E)", consequence of Edit M9.

**Edit R4 — every superscript exponent flattened to plain text.** 8 sites: `2.3 × 10⁻⁹¹` →
`2.3 × 10-91`, `6.4 × 10⁻⁴¹` → `6.4 × 10-41`, `1.2 × 10⁻⁶` → `1.2 × 10-6`, etc. Unicode superscripts
do not survive Word/journal typesetting reliably and cannot be found by a text search for `10-91`.
The manuscript's own convention (documented in `revision_G3/CLAUDE.md`) is plain digits.

**Edit R5 — the single largest class of edit: removal of every "what changed since the last
version" statement.** Six sites, all deleted or rewritten. Counts across the whole document:

| Phrase | Claude | Daniil |
|---|---|---|
| "no longer" | 6 | 0 |
| "not retained" | 2 | 0 |
| "previous threshold" | 3 | 0 |
| "at the tightened threshold" | 3 | 1 |
| "rather than assumed" | 1 | 0 |
| "measured rather than [visual]" | 1 | 0 |
| "is what makes" | 1 | 0 |
| "cannot establish" | 1 | 0 |

Worked examples:

> **Claude:** At the **tightened** threshold the LINE class **loses** flavone metabolic process
> (FDR = 0.088) and the LTR class **loses** glutamatergic synapse and positive regulation of
> lipopolysaccharide-mediated signalling (both FDR = 0.086).
> **Daniil:** At the **0.05 FDR** threshold the LINE class **had insignificant term of** flavone
> metabolic process (FDR = 0.088) and the LTR **had** glutamatergic synapse and positive regulation
> of lipopolysaccharide-mediated signaling (both FDR = 0.086), **that could be deemed significant
> under a more relaxed FDR threshold of 0.1 but are rejected here.**

The information content is *higher*, not lower: the reader now learns what the alternative threshold
would have done, without being told a story about a draft they never saw.

> **Claude:** Two associations that were significant at the previous threshold do not survive at 0.05
> and are no longer claimed: lipid metabolism in LINE-adjacent genes (FDR = 0.021 at 0.1) and DNA
> replication and recombination in the all-TE top group (FDR = 0.036 at 0.1).
> **Daniil:** *[deleted]*

> **Claude:** Dong-R4, whose single term the previous threshold retained, no longer reaches
> significance, which is consistent with the random nature we ascribed to it.
> **Daniil:** *[deleted]*

> **Claude:** the metals metabolism association of MIR elements, significant at the previous
> threshold (FDR = 0.045), **is not retained at 0.05**.
> **Daniil:** *[deleted]*

> **Claude:** ERV1 elements had a statistically significantly overrepresented other metabolism group
> (FDR = 0.011) **but no longer a significant lipids metabolism group (FDR = 0.045 at 0.1)**
> **Daniil:** ERV1 elements had a statistically significantly overrepresented other metabolism group
> (FDR = 0.011)

> **Claude:** cell adhesion **is no longer represented at all**.
> **Daniil:** cell adhesion **was not represented**.

> **Claude:** and **at the tightened threshold** the difference in GO term count … **is no longer
> significant** either
> **Daniil:** and the difference in GO term count … **is not significant** either

**Edit R6 — an editorial verdict replaced by a count.**

> **Claude:** **The metal claim is now zinc-only:** cellular response to cadmium ion and detoxification
> of copper ion **both have** FDR = 0.078 **and are not retained**, as are macrophage activation
> (FDR = 0.053), phosphatidylinositol-4,5-bisphosphate binding, negative regulation of response to
> cytokine stimulus (both FDR = 0.062) and the complement component C1q complex (FDR = 0.078).
> **Daniil:** **6 GO terms with FDR between 0.05 and 0.1 were rejected as insignificant:** cellular
> response to cadmium ion and detoxification of copper ion (FDR = 0.078 for both), as are macrophage
> activation (FDR = 0.053), …

"The metal claim is now zinc-only" is a statement about the paper. "6 GO terms with FDR between 0.05
and 0.1 were rejected" is a statement about the data, and it is countable — which is how the error
in Part IV.6 became visible.

**Edit R7 — a figure citation replaced with the object that actually contains the number.**

> **Claude:** LTR families demonstrated 33 unique GO terms compared to 22 in LTRs as a class
> **(Figure 5A)**.
> **Daniil:** … **(File S3)**.

**This is a corrected error, not a style change.** Figure 5A is the connection map for the
*divergence-stratified class-level* gene sets. It cannot show a family-level term count. The number
is verifiable in File S3:

```python
fam = pd.read_csv('output/GO_families_fdr005.csv')
cls = pd.read_csv('output/GO_classes_count_fdr005.csv')
fam[fam.class_name=='LTR']['Term ID'].nunique()   # 33  ✓
cls[cls.class_name=='LTR']['Term ID'].nunique()   # 22  ✓
```

**Edit R8 — a claim about a figure that the figure cannot carry.**

> **Claude:** beta-2-microglobulin binding and antigen processing … both have FDR = 0.051 and are no
> longer retained, **so the MHC association rests on one term rather than three (Figure 6A)**.
> **Daniil:** … both have FDR = 0.051 and **were rejected by FDR**.

Figure 6A shows up to five terms per family; "one rather than three" is a statement about the GO
table, and the parenthetical pointed at the wrong object.

**Edit R9 — added the missing pointer to the full network.**

> **Daniil:** Visualisation of up to five GO terms per family in Figure 6A showed a high degree of
> functional distinction between processes by families **(the full network is shown in Figure S10)**.

Figure 6A is a five-terms-per-family subset. The reader is told, at the point of the claim, where
the unabridged version is. **P2.**

**Edit R10 — pronoun and antecedent repairs.**

> **Claude:** **The latter term** genes that had SVA elements in their vicinity were POLR2A…
> **Daniil:** **The termination of RNA polymerase II transcription** genes that had SVA elements…

> **Claude:** These immune system GO terms all **share** a single core interferon gene set, which is
> characterised in detail at the end of **this section**.
> **Daniil:** … all **shared** a single core interferon gene set, which is characterized in detail at
> the end of **the Results section**.

"The latter" required the reader to hold two GO term names in working memory across a sentence
boundary. "this section" is ambiguous once Methods precedes Results. Tense unified to past.

**Edit R11 — set-membership language made exact.**

> **Claude:** Two of **the remaining** families (hAT-Tip100 and L2) have only 1 significant GO term
> each, **again** reflecting…
> **Daniil:** Two of the families **with significant GO terms** (hAT-Tip100 and L2) have only 1 term
> each, reflecting…

"the remaining" only had a referent because of the deleted Dong-R4 sentence. Verified: 13 families
have ≥ 1 term at FDR 0.05; hAT-Tip100 and L2 have exactly 1 each.

**Edit R12 — statistics-reporting hygiene.** "(16)" → "(16 terms)" so a bare integer in parentheses
is not read as a citation. "Mann-Whitney U, raw p" → "Mann-Whitney U raw p" (the comma made "U" look
like a variable). "medians 12 and 2 for 11 and 2 families" → "… for 11 and 2 families,
**respectively**". "falling **just** outside the 0.05 threshold" → "falling outside" ("just" is an
editorial nudge toward significance). "according to **the FDR-adjusted** Fisher exact test" →
"according to Fisher exact test" (the FDR adjustment is already stated globally in Methods and the
value quoted is already an FDR).

## I.5 Discussion

**Edit D1 — citation renumbering**, 12 further sites (20–31 → 49–60; 84–85 → 113–114; 8,10,11,14,15
→ 8,10,11,43,44).

**Edit D2 — 28 orphan citation-only paragraphs appeared.** Paragraphs consisting of nothing but
citation markers, e.g. `(9)(10)(11)(8)(5,6,18,42,45,46)`, `(50,51)(60)(59)(61)(62)(63)`,
`(104)(105,106)(107)(108)(109)(110)(111)(112)`. These are Mendeley content controls left behind when
the passages moved to File S6 were deleted and the deletions accepted. **Not an editorial choice — a
mechanical artefact.** See Part IV.2.

## I.6 Back matter

**Edit B1 — an "Ethical statement" section was added** before Literature cited:

> Ethical statement. This study represents a purely computational analysis based on publicly
> available genomic datasets. All data analyzed were obtained from the T2T genome assembly and the
> RepeatMasker track

The sentence is truncated mid-clause. See Part IV.5.

**Edit B2 — seven editor notes remain in the document** (Claude wrote nine; two were consumed by
Edits I3 and R1). All seven must be removed before submission. See Part IV.3.

## I.7 Global mechanical edits

**Spelling normalised to US English**, 6 sites: `neighbourhoods` → `neighborhoods`, `characterised`
→ `characterized`, `signalling` → `signaling`, `randomised` → `randomized`, `co-localised` →
`co-localized`. (`Visualisation` survives in one place — see Part IV.7.)

**Capitalisation:** `P-value` → `p-value`.

---

# Part II — The editing principles

Each principle is stated, justified from the evidence above, and given a **detection rule** (how to
find violations mechanically) and a **rewrite recipe** (what to do instead).

---

## P1 — No proofless, vague, or self-promoting statements

**Statement.** Every sentence in the manuscript must assert something about the data, the method, or
the literature that a reader can check. Sentences that assert something about the *quality of the
work* — its rigour, its robustness, its sufficiency, its having been verified — are deleted, or
replaced by the measurement that would justify them.

**Why.** A reader cannot falsify "this was verified rather than assumed". They can falsify "the
running mean differs from the N = 500 value by less than 1 % by N = 140". The second sentence
carries all the persuasive force of the first *and* is checkable. The first adds nothing except a
target for a referee.

### P1a — The banned constructions

Delete on sight, or replace with the measurement:

| Construction | Example from this pass | Replacement |
|---|---|---|
| "X rather than assumed / rather than guessed" | "this was verified rather than assumed" | the verification number |
| "X rather than a visual statement" | "so that 'concordant' is a measured rather than a visual statement" | *delete* — the cited tests already say it |
| "X is what makes Y interpretable / meaningful" | "The permutation background is what makes these enrichment values interpretable" | state the bias it removes, once |
| "high / strong robustness", "reliable baseline", "solid evidence" | — | the robustness statistic |
| "notably", "strikingly", "importantly", "remarkably" | — | *delete* |
| "clearly shows", "demonstrates convincingly" | — | "shows" |
| Editorial verdicts about the paper | "The metal claim is now zinc-only" | the count and the terms |
| Diligence claims | "so that a silent skip is impossible" | *delete from the manuscript* — belongs in code docs |

### P1b — Unsupported claims about the literature are also proofless

"whether these elements affect innate immune gene regulation **remains untested**" asserts a
negative about the entire literature. Unless a systematic search was performed and can be cited,
this is a proofless statement in the same sense. The five prior window-size studies described with
no citations (Edit M2) are the same failure in the positive direction.

### P1c — What is *not* banned: labelled hypotheses

"could indicate a recent evolutionary arms-race" survived. The distinction:

| Banned | Allowed |
|---|---|
| Unfalsifiable self-assessment of the work | Explicitly modal claim about the biology |
| "high robustness", "reliable" | "could indicate", "is consistent with", "we hypothesise" |
| Asserted without a stated limitation | Paired with an explicit statement of what the design cannot show |

The test: **would a referee know what evidence would refute it?** "could indicate an arms race,
though the design is correlative and cannot detect one" is refutable in principle and is honest
about its status. "high robustness" is not.

### P1d — Detection rule

```bash
grep -n -i -E "rather than (assumed|guessed|a visual)|is what (makes|justifies)|robustness|reliab|\
notably|strikingly|remarkably|importantly|clearly (shows|demonstrates)|convincing|it is worth noting|\
we emphasi[sz]e|it should be noted|serves as a (solid|strong|reliable)" manuscript.txt
```

**Rewrite recipe.** For each hit, ask: *what number would make this sentence unnecessary?* If the
number exists, put the number in and delete the sentence. If the number does not exist, delete the
sentence.

---

## P2 — Every number, GO term and named result is traceable to a specific source object

**Statement.** Every quantitative claim carries a pointer to exactly one object that contains it:
a named supplementary file **and sheet**, a specific figure **panel**, or a numbered table. The
pointer must be to the object that *actually holds* the number, not to a related one.

**Why.** This is what makes a paper auditable. It is also, empirically, the way errors get caught:
Edit R7 was a wrong pointer, and the wrong pointer was only visible because the correct one was
demanded.

### P2a — Granularity

| Too coarse | Correct |
|---|---|
| "(File S3)" for a specific term | "(File S3, sheet GO_by_family)" |
| "(Figure S13B-D)" for three claims | "(Figure S1A)", "(Figure S1B)", "(Figure S1C)" — one per claim |
| "in the supplementary tables" | "`enrichment_families_with_random.csv` in the GitHub repository (URL)" |
| "the checkpoint values are in the repository" | acceptable only for material that supports no numbered claim |

### P2b — The pointer must match the object type

A **figure** can support: a distribution, a trend, a visible difference, a network topology.
A **table or workbook sheet** must support: an exact count, an exact p-value, a set membership.

Claude's Edits R7 and R8 both cited a figure for a count. Both were corrected to the table. The
general rule: **if the claim contains an integer or a p-value, the pointer goes to a table.**

### P2c — Load-bearing numbers stay in the main text

Edit R1 reverses a decision to move unadjusted Fisher p-values to a supplementary workbook. The
principle: a supplementary file adds detail; it never becomes the only home of a number the main
text argues from. If the text says a class is enriched, the printed table shows both the raw and
the adjusted p.

### P2d — Methods must cite their methods

Five statistical procedures (Spearman, Pearson, Bland-Altman, overlap coefficient + hypergeometric,
Kendall tau, label-shuffling permutation) went from zero citations to five (Edit M11). Five prior
studies went from zero to five (Edit M2). Software versions are pinned and cited. The AI tools are
pinned to a version and their prompt artefacts pointed at a public repository (Edit M14).

### P2e — Detection rule

```bash
# every sentence containing a number, FDR, p-value, OR percentage — does it carry a pointer?
grep -n -E "FDR = |p = |raw p|OR |[0-9]+ %|× 10-" manuscript.txt \
  | grep -v -E "\((Figure|File|Table)[^)]*\)"
```

Any line that survives both filters is a number with no source. Then, for each pointer, open the
object and confirm the number is in it.

---

## P3 — Sequential numbering by order of first citation

**Statement.** Figures, supplementary figures, tables, supplementary files, panels within a figure,
and references are numbered in the order in which the text first cites them. No exceptions, and the
rule applies recursively to panel letters.

**Why.** A reader following the text should never have to jump backwards. It is also the house rule
of essentially every journal, so violations are returned by the production editor.

**Evidence.** 22 citation-renumbering edits (Edit M1, D1). The whole reference list was renumbered
after the Methods-before-Results reorder moved 24 method citations ahead of everything else. Figure
S13 became Figure S1 because it is now cited first, in Methods (Edit M9); the length-distribution
panels that were S1A/S1B became S1D/S1E (Edit R3).

**Corollary — a structural move triggers a full renumbering pass.** Moving Materials and methods
ahead of Results is a one-line decision that invalidates the numbering of the entire reference list
and at least one supplementary figure. Any agent that performs a section move **must** re-run the
numbering audit and report the new order; it is not a follow-up task.

**Detection rule.**

```python
import re
pat = re.compile(r'(Figure S\d+[A-Z]?|Figure \d+[A-Z]?|File S\d+|Table \d+)')
order, seen = [], set()
for tok in pat.findall(open('manuscript.txt').read()):
    base = re.match(r'((?:Figure S|Figure |File S|Table )\d+)', tok).group(1)
    if base not in seen:
        seen.add(base); order.append(base)
print(order)   # must be S1, S2, S3, … and 1, 2, 3, … each in ascending order within its series
```

Run this on the *body text only*, **excluding the legends block**, or the legends mask the
violation — a legend header looks exactly like a citation and every object then appears cited.
The audit is implemented, with that exclusion and the legend-inventory cross-check, as

```bash
~/venvs/collagen_3_11/bin/python \
  ~/.claude/skills/scientific-review/scientific_review_tools.py numbering ms.txt
```

Main figures 1–10 and Tables 1–2 are in order. **Supplementary figures and files are not**, and
one supplementary figure is uncited — see Part IV.1.

---

## P4 — Scientific language: precise, plain, one subtheme per passage

**Statement.** No casual register, no over-built grammar. One idea per sentence, one subtheme per
paragraph, one internally consistent block of results per subsection. Terminology is single-valued
across the whole manuscript.

### P4a — Register

| Casual / imprecise | Corrected |
|---|---|
| "makes it possible to measure where TEs **sit** relative to…" | "enables comprehensive investigation of TE contributions to…" |
| "how close each TE group **sits** to…" | "how close each TE group **locates** to…" |
| "falling **just** outside the 0.05 threshold" | "falling outside the 0.05 threshold" |
| "**again** reflecting the likely random nature" | "reflecting the likely random nature" |
| "cell adhesion **is no longer represented at all**" | "cell adhesion **was not represented**" |

### P4b — Structure, not decoration

The M7 rewrite is the model:

* six sentences → one;
* `Reporting X would **therefore** …, **and** every statement is made on Y` → `**Since** reporting X
  would …, every statement is made on Y` — reason first, consequence second, both connectives gone;
* the definition stated once, in one vocabulary ("this enrichment score"), not twice in two.

### P4c — Terminology is single-valued

"class of gene" was replaced by "functional gene group" (Edit I2) because "class" is reserved for
TE classes. "the latter term" replaced by the term's name (Edit R10). "the choice" replaced by "the
10kb window choice" (Edit M3). **Rule: a noun phrase that has a technical meaning elsewhere in the
paper may not be reused with a general meaning, and a demonstrative may not reach back more than
one clause.**

### P4d — One subtheme per passage

Three structural moves in this pass:

1. "Statistical tests" moved from between the permutation background and the IFNA test, to after
   both. Result: the two permutation-based methods are adjacent, generic statistics follow (M10).
2. The concordance-methods sentence moved out of the sensitivity subsection into "Statistical
   tests", where the other test descriptions live (M11).
3. The "Sensitivity to the window and to the gene-set cut-off" pseudo-heading deleted and its
   content folded into the paragraph that introduces the window (M3).

The test for a subsection: **can you state, in one clause, what block of results it contains, with
no "and also"?** If not, split it or move the outlier.

### P4e — Tense and agreement

Results are past tense ("all **shared** a single core interferon gene set"). Comparative statements
carry "respectively" when two lists are paired. A bare integer in parentheses is written with its
unit ("16 terms") so it is not read as a citation.

---

## P5 — Every number is recalculated from the source; no rounded ranges

**Statement.** Before a number enters the manuscript, it is recomputed from the supplementary table
or the repository code that produced it. Approximations, ranges and "roughly N" are not acceptable
in a final version.

**Why.** Numbers propagate: a figure caption, a Results sentence, an abstract and a supplementary
sheet can drift apart across drafts. The only defence is to re-derive at the point of writing.

### P5a — Quote the true crossing point, not the nearest grid row

Edit M8 is the clearest case. `output/permutation_convergence_checkpoints.csv` has rows at
N = 50/100/200/250/300/400. Claude quoted N = 250 and N = 100 because those rows exist. Daniil
quoted N = 140, N = 250 and N = 315 — where the curves actually cross 1 %, 6 % and 10 %. **A
checkpoint table is a summary of a curve; quote the curve.**

### P5b — Quote the number in the unit the reader can verify

The same edit converted "0.06 standard deviations of drift" into "less than 1 % of the N = 500
value", because the y-axis of the panel being cited is a fraction of the N = 500 value. **The unit
in the text must be the unit on the axis of the figure it cites.**

### P5c — Nominal vs realised counts

"(2,872 genes)" → "(2,872 genes **at maximum**)" (Edit M5). The top-10 % cut is a ceiling; ties at
the boundary make the realised set smaller for some TE groups. **A cut-off is not a count.**

### P5d — Counts are stated, and therefore checkable

Edit R6 replaced "The metal claim is now zinc-only" with "6 GO terms with FDR between 0.05 and 0.1
were rejected". This is the principle working: stating the count made it auditable, and the audit
(Part IV.6) shows the count is wrong — there are 10 such terms for MIR. A vague sentence would have
hidden the error indefinitely.

### P5e — The recalculation protocol

For each number in the manuscript:

1. Identify the file that produced it (`output/*.csv`, `output/*.json`, `supplementary/*.xlsx`).
2. Recompute it in a fresh interpreter — do not trust a previously printed value, a `results.md`, or
   a `CLAUDE.md` "numbers worth not re-deriving" table. Those tables have been wrong before: this
   project's own `CLAUDE.md` recorded the length-bias correlation as R = 0.661 for months when the
   measured value is R = 0.985.
3. Confirm the rounding matches the manuscript's convention (3 significant figures for OR, 2 for
   percentages, exponent form for p < 0.001).
4. Confirm the same number in every other place it appears — abstract, results, caption, table,
   supplementary sheet.

**Detection rule.** Every number in the manuscript should appear at least once in
`output/` or `supplementary/`. Extract all numeric tokens from the manuscript and grep them against
the output tree; anything with no hit needs a derivation or must go.

---

## P6 — The manuscript reports the current state of the analysis, not its history

**Statement.** The published paper describes one analysis with one set of thresholds. It never
refers to a previous draft, a previous threshold, a term that "no longer" qualifies, or a claim that
is "no longer made".

**Why.** A reader of the published version has never seen the previous threshold. "Dong-R4 no longer
reaches significance" is meaningless to them and reads, to a referee, as an admission that results
were unstable across drafts. The revision history belongs in the response-to-reviewers letter.

**Evidence.** The largest single class of edit in this pass — see the frequency table under Edit R5.
Six deletions or rewrites; "no longer" went from 6 occurrences to 0.

**The one legitimate use.** A comparison to an alternative threshold, stated as a property of the
data rather than of the draft:

> "…(both FDR = 0.086), that could be deemed significant under a more relaxed FDR threshold of 0.1
> but are rejected here."

This tells the reader what a different analytical choice would have produced. That is a sensitivity
statement, not a changelog entry.

**Detection rule.**

```bash
grep -n -i -E "no longer|not retained|previous(ly)? threshold|at the tightened|\
in the (earlier|previous) (version|draft|analysis)|do(es)? not survive|we now |is now " manuscript.txt
```

**Rewrite recipe.** Ask: *would this sentence make sense to a reader who has never seen another
version of this paper?* If no — delete it, or restate it as a property of the data under a named
alternative threshold.

---

## P7 — No process scaffolding survives into the manuscript

**Statement.** Editor notes, TODOs, placeholders, instructions to the typesetter, and explanations
of how a script performed an edit do not belong in a document that will be submitted.

**Evidence.** Claude left nine `[EDITOR NOTE]` blocks. Daniil consumed two by replacing them with
the content they described (Edits I3 and R1). Seven remain and must be removed before submission
(Part IV.3).

The two he acted on are instructive: both notes were *describing content that should have been
there*. The correct resolution of an editor note is almost always **do the thing**, not **leave a
note about the thing**.

**Corollary for agents.** If you must leave an instruction for the human, put it in the run report,
the implementation plan or a separate markdown file — never in the `.docx`. If a note must
temporarily live in the document, it must be findable by one unambiguous string and the count must
be reported to the user at the end of the run.

**Detection rule.**

```bash
grep -n -E "\[EDITOR NOTE\]|\[TODO\]|\[ZENODO DOI\]|\(ref\)|XXX|TBD|FIXME|<placeholder>" manuscript.txt
```

Expected count at submission: **0**.

---

## P8 — Caveats are stated once, in the place the reader expects, and scoped forward

**Statement.** A limitation is stated in Limitations. A statement of scope goes at the end of the
Introduction. Neither is repeated defensively in the middle of a paragraph doing other work.

**Evidence.** Claude wrote the arms-race caveat twice: once mid-paragraph in the Introduction's
opening (Edit I1), once in Limitations. Daniil deleted the first, kept the second untouched, and
added a *third-position* sentence at the end of the Introduction (Edit I3) that is not a caveat but
a scope statement:

> Albeit the current analysis is purely correlative and cannot detect the ongoing evolutionary arms
> race, it points at genome loci and molecular processes that could be impacted by it for future
> investigation.

The grammar carries the distinction: the limitation is in a subordinate clause, the contribution in
the main clause. A defensive disclaimer inverts that and reads as an apology.

**Rewrite recipe for a caveat.** `Although [what the design cannot do], [what it does do and what
that enables].` Never the reverse, never twice, never mid-paragraph.

---

## P9 — House style is enforced mechanically

These are the conventions this manuscript follows. They are not negotiable per-sentence.

| Item | Convention | Sites in this pass |
|---|---|---|
| Spelling | US English | 6 |
| Exponents | plain text, `2.3 × 10-91` — never Unicode superscript | 8 |
| p-value | lower-case `p` everywhere, including `p-value`, `raw p`, `FDR` | 1 |
| Bare integers in parentheses | give the unit: `(16 terms)`, never `(16)` | 1 |
| Test names | `Mann-Whitney U raw p = …` — no comma between the test and `raw p` | 2 |
| Paired lists | close with `respectively` | 1 |
| Tables | placed at first citation, with both raw and adjusted p | 2 |

Unicode superscripts are the most damaging of these: they do not survive typesetting reliably and
they make the number unfindable by a plain-text search, which defeats P5's verification protocol.

---

## P10 — Priority order when principles conflict

1. **Correctness** (P5) — a wrong number outranks everything.
2. **Traceability** (P2) — a right number with no source is not usable.
3. **Numbering** (P3) — a right, sourced number the reader cannot navigate to is not usable either.
4. **No proofless statements** (P1) and **no history** (P6) — deletions, always safe to apply.
5. **Language** (P4) and **house style** (P9) — last, because they touch the most text for the least
   risk-reduction.

An agent working under time pressure should do 1–4 completely and 5 partially, never the reverse.

---

# Part III — Verification protocol

Run this before handing a manuscript back. Every step produces a number that goes in the run report.

### III.1 Structural

```bash
# 1. Editor notes / placeholders: expect 0
grep -c -E "\[EDITOR NOTE\]|\[TODO\]|\(ref\)|\[ZENODO DOI\]|TBD|FIXME" manuscript.txt

# 2. Orphan citation paragraphs: expect 0
grep -cE '^(\([0-9,–-]+\))+$' manuscript.txt

# 3. Citation order: must be ascending in the body
python refs_order.py manuscript.txt

# 4. Legend inventory matches body citations: every Figure SN cited has a legend and vice versa
```

### III.2 Numerical

For every number in the text:

```python
import pandas as pd
# example: the LTR family/class GO term counts
fam = pd.read_csv('output/GO_families_fdr005.csv')
cls = pd.read_csv('output/GO_classes_count_fdr005.csv')
assert fam[fam.class_name=='LTR']['Term ID'].nunique() == 33
assert cls[cls.class_name=='LTR']['Term ID'].nunique() == 22
assert fam.family_name.nunique() == 13          # families with >= 1 term at 0.05
```

Where a manuscript sentence enumerates a set ("6 GO terms with FDR between 0.05 and 0.1"), the
enumeration must be reconstructed by query, not by counting the items in the sentence:

```python
r = pd.read_csv('output/GO_families_fdr01_reference.csv')
band = r[(r.family_name=='MIR') & (r.FDR > 0.05) & (r.FDR <= 0.1)]
len(band)   # <-- this is the number the sentence must state
```

### III.3 Language

```bash
grep -n -i -E "rather than assumed|is what makes|measured rather than|robustness|reliab|\
no longer|not retained|previous threshold|tightened|notably|strikingly|remarkably" manuscript.txt
grep -n -P "[\x{2070}-\x{209F}]" manuscript.txt          # Unicode superscripts: expect none
grep -n -E "neighbourhood|characterised|signalling|randomised|colouri?s|visualis" manuscript.txt
```

### III.4 Document integrity (docx-specific)

```bash
# Mendeley content controls must survive every scripted edit
python -c "import zipfile; x=zipfile.ZipFile('m.docx').read('word/document.xml').decode();\
print(x.count('<w:sdt>'), x.count('MENDELEY_CITATION'))"
# no <w:t> inside a <w:del>; the count of <w:ins>/<w:del> before and after each script
```

---

# Part IV — Open items exposed by this comparison

These are unresolved as of 2026-08-09. Each is a concrete action.

**IV.1 — Supplementary figure and file numbering is not sequential, and the legends were not
renumbered with the body.**

Re-audited 2026-08-09 with legends correctly excluded (`scientific_review_tools.py numbering`),
which changes the picture from the first scan in Part 0 — that one counted the legend block as
body text and so credited three supplementary figures with citations they do not have:

| Series | Order of first citation in body prose | Verdict |
|---|---|---|
| Figure | 1 2 3 4 5 6 7 8 9 10 | OK |
| Table | 1 2 | OK |
| File S | 1 2 3 **5 6 4** | not ascending |
| Figure S | 1 2 3 4 5 6 7 **10 13 12** | not ascending |

Two further findings only the legend-excluded scan makes visible:

* **`Figure S8` and `Figure S9` are cited only from inside the legends of Figures 4 and 5**, never
  from body prose. That is not automatically wrong — their position in the sequence is the position
  of the figure whose legend cites them — but it must be confirmed deliberately rather than by
  accident. (`Figure S10` was in the same state in Claude's version; Daniil's Edit R9 gave it a
  body citation.)
* **`Figure S11` has a legend and is cited nowhere at all** — not in the body, not in any legend.
  An uncited supplementary figure is returned by the production editor.

Worse, the renumbering is half-applied. The Methods now cite `Figure S1A/B/C` for the convergence
curves and Results cites `Figure S1D/E` for the length distributions — but the legend block still
reads:

* `Figure S1. Ridge plots for length distribution comparison … (A) for all classes and (B) for
  individual classes.`
* `Figure S13. Stability of the gene sets across TSS window sizes, and convergence of the
  permutation background at 500 permutations. (A) … (B) … (C) … (D) …`

**Action.** Decide the intended S1 (merge convergence A–C with length D–E into one figure, or keep
them separate and renumber all of S1–S13); give `Figure S11` a citation or drop it; then renumber
the whole supplementary series by first citation, update every legend, update every panel letter,
and re-export the affected Figma frames. Note that the SVG filenames in `svg/` are deliberately
*not* renamed (the notebooks write them), so `svg/PLACEMENT.md` must carry the mapping.

**IV.2 — 28 paragraphs consist only of orphaned citation markers.**

Examples: `(9)(10)(11)(8)(5,6,18,42,45,46)`, `(50,51)(60)(59)(61)(62)(63)`,
`(104)(105,106)(107)(108)(109)(110)(111)(112)`. They sit where the extended-discussion passages were
deleted before being moved to File S6. The Mendeley content controls were not inside the tracked
deletions, so accepting the deletions left them behind.

**Action.** Delete these paragraphs in Word, then refresh the Mendeley bibliography — deleting them
will change the reference list, so the renumbering audit (IV.1, P3) must be re-run afterwards.
Confirm afterwards that every reference in the bibliography is still cited somewhere.

**IV.3 — Seven `[EDITOR NOTE]` blocks remain.** At body-text lines 195, 196, 259, 261, 285, 294, 297
of the extract (subsection-order notes, back-matter order, G2 legend placement, G14 acknowledgments,
the Zenodo DOI placeholder, and the G8 table-placement note). Each is either an instruction that has
already been carried out — delete — or one that has not — do it, then delete.

**IV.4 — Tables 1 and 2 now appear twice.** Once inline in Results (Daniil's insertion, with the raw
Fisher p column and mean/SD split) and once at the end of the main text (Claude's version, four
columns, `mean ± SD` combined). The two versions also disagree on column names and on the Helitrons
Fisher p (`0.41` vs `0.413`).

**Action.** Keep one. G3 house style puts tables at the end of the main text, but the reviewer asked
for the table to be readable; whichever is kept, delete the other and make the caption match
("gene TSS neighborhoods: counts and odds ratios" vs "gene TSS 10 kb neighbourhoods").

**IV.5 — Two incomplete insertions.**

* `Claude Code Opus 5.0 (ref)` — the placeholder needs a real citation.
* The Ethical statement ends mid-clause: *"All data analyzed were obtained from the T2T genome
  assembly and the RepeatMasker track"*. Complete the sentence and add the standard no-human-subjects
  / no-animal-subjects wording.

**IV.6 — Two enumerated GO-term counts do not match the source table.** Found by applying P5's
recalculation protocol to Daniil's own edits.

*MIR, borderline terms.* The text states **6** GO terms with FDR between 0.05 and 0.1 were rejected.
`output/GO_families_fdr01_reference.csv` gives **10** for MIR in that band:

| Term | FDR |
|---|---|
| macrophage activation | 0.053 |
| actin binding | 0.062 |
| phosphatidylinositol-4,5-bisphosphate binding | 0.062 |
| negative regulation of response to cytokine stimulus | 0.062 |
| cellular response to cadmium ion | 0.078 |
| detoxification of copper ion | 0.078 |
| complement component C1q complex | 0.078 |
| **neuromuscular junction** | **0.078** |
| **cell junction** | **0.098** |
| **actin filament binding** | **0.098** |

The four in bold are not in the sentence. Either state "10" and list all ten, or state "6 of the 10"
and say the listed ones are those with FDR ≤ 0.078, or restrict the stated band.

*LTR class, borderline terms.* The text names **two** LTR terms in the 0.05–0.1 band (glutamatergic
synapse, positive regulation of lipopolysaccharide-mediated signaling, both FDR = 0.086).
`output/GO_classes_count_fdr01_reference.csv` gives **three**: `nucleotide binding` (FDR = 0.076) is
missing. (LINE is correct: exactly one term, flavone metabolic process, FDR = 0.088.)

**IV.7 — The convergence numbers in Edit M8 need a source table cut to their definition.** The
manuscript now claims: running mean within 1 % of the N = 500 value by **N = 140**; running SD within
6 % by **N = 250**; least stable family's SD within 10 % after **N = 315**.
`output/permutation_convergence_checkpoints.csv` reports drift in *standard-deviation units* on a
grid of N = 50/100/200/250/300/400 and does not contain N = 140 or N = 315.

**Action.** Extend the convergence script to emit, per N, the running mean and running SD as a
fraction of the N = 500 value (class-level and family-level), write it to
`output/permutation_convergence_curves.csv`, and confirm the three crossing points. Cite that file
in the Methods sentence that says "the checkpoint values … are provided in the repository".

**IV.8 — Residual house-style hits.** `Visualisation` survives at one Results site (Edit R9's
sentence begins "Visualisation of up to five GO terms per family"); `neighbourhoods` survives in
both end-of-document table captions.

---

# Part V — Working checklist for a future review pass

Use in this order. Each line is a pass over the whole document.

### Before writing anything

- [ ] Read `revision_G3/CLAUDE.md` and `project_overview.md` for the frozen files, the decisions and
      the caveats.
- [ ] Extract the accept-all text of the current working document; that is the object under review.
- [ ] Build the list of every number the document quotes and the file each comes from.

### Pass 1 — correctness (P5)

- [ ] Recompute every number from `output/` or `supplementary/` in a fresh interpreter.
- [ ] For every enumerated set ("N terms", "N families"), reconstruct the set by query, not by
      counting the items written in the sentence.
- [ ] Check each number against every other place it appears (abstract, results, caption, table).
- [ ] No rounded ranges, no "roughly", no "approximately" in a final version.
- [ ] Cut-offs are labelled as cut-offs ("at maximum"), not as realised counts.

### Pass 2 — traceability (P2)

- [ ] Every claim containing a number carries a pointer to one object.
- [ ] Integers and p-values point at tables/sheets; distributions and trends point at panels.
- [ ] The pointed-at object actually contains the number — open it and check.
- [ ] Every method (statistical, computational, AI) is cited with a version.
- [ ] Every prior-work claim is cited.
- [ ] Nothing the main text argues from lives only in a supplementary workbook.

### Pass 3 — numbering (P3)

- [ ] Run the citation-order scan on body text with legends excluded.
- [ ] Figures, supplementary figures, tables, supplementary files and panel letters all ascend by
      first citation.
- [ ] Every legend exists, matches its number, and enumerates the panels the body cites.
- [ ] References renumbered after any structural move.

### Pass 4 — deletions (P1, P6, P7)

- [ ] Run the P1 grep. For each hit: replace with the measurement, or delete.
- [ ] Run the P6 grep. Delete every reference to a previous version or threshold, except sensitivity
      statements phrased as properties of the data.
- [ ] Run the P7 grep. Zero editor notes, TODOs, placeholders.
- [ ] Zero orphan citation-only paragraphs.

### Pass 5 — language and style (P4, P8, P9)

- [ ] Each paragraph states one subtheme; each subsection one internally consistent block.
- [ ] No demonstrative ("the latter", "this", "the choice") reaching back more than one clause.
- [ ] No technical noun reused with a general meaning.
- [ ] Caveats once, in Limitations or at the end of the Introduction, subordinate clause first.
- [ ] US spelling; plain-text exponents; lower-case `p`; units on parenthetical integers;
      `respectively` on paired lists.

### Pass 6 — document integrity

- [ ] Mendeley content-control count unchanged (or the change explained).
- [ ] No `<w:t>` inside a `<w:del>`.
- [ ] `<w:ins>` / `<w:del>` counts reported before and after every scripted edit.
- [ ] Frozen-file md5s still match.

---

## Appendix — the one-paragraph version

Write only what a reader can check. Every number points at the table or panel that holds it, and
was recomputed from that source before it was typed. Everything is numbered in the order the text
first mentions it. Say it once, plainly, in the section where it belongs. Never tell the reader what
the paper used to say, and never tell them how careful you were — show them the number that makes
the question moot.
