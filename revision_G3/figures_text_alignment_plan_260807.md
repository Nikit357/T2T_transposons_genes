# Figure–text alignment plan — final pass before submission

**Manuscript:** G3-2026-406828
**Input document (read-only from here on):**
`revision_G3/Revised_manuscript/T2T_genes_article_G3_revision_260804_manual.docx`
**Output document (new):**
`revision_G3/Revised_manuscript/T2T_genes_article_G3_revision_260807.docx`
**Written:** 2026-08-07 · **Rewritten:** 2026-08-07 against the exported figures, then against
Daniil's decisions on D-1, D-2 and D-3 · **Target submission:** on or before 2026-08-28

**Evidence base:** the 23 exported PDFs in `revision_G3/current_figures_260807/`, matched
panel-by-panel against `svg/*.svg` by de-spaced label-set overlap, with six figures rendered and
inspected visually. **Where the exported figures and `svg/PLACEMENT.md` disagree, the exported
figures win** — they are what will be submitted.

---

## Overview

The figures are finished. Two structural changes were made during layout that `PLACEMENT.md`
does not record, and both propagate into the manuscript:

1. **The supplementary set is now S1–S13, not S1–S14.** Old S6 (the UCSC view of the
   interferon-alpha domain) was consumed by the new main Figure 8, and everything from S7 up
   shifted down by one.
2. **Figure 8 has four panels, not three.** The leave-one-out plot that `PLACEMENT.md` listed as
   an optional inset to 8C was promoted to a full panel D.

The consequence is that the manuscript's supplementary citations are now **wrong rather than
merely stale**. `Figure S6` in the Results still points at the interferon-alpha browser view; in
the exported set, Figure S6 is the log-scaled TSS-by-family distribution. Seven body paragraphs,
three main-text figure captions and nine supplementary captions are affected. This is the
highest-priority item in the plan: a reader following `Figure S6` today lands on an unrelated
figure.

Beyond the renumbering, the review found:

- **three body sentences that describe a figure that no longer exists** — the network panels are
  described as "top 30" terms where they now show 10, 10 and **5**;
- **Figure 8 has no caption at all**, and is cited once, without panel letters, and *before*
  Figures 6 and 7;
- **the `-log10(FDR)` colourbar is missing from three network panels** that promise it;
- **an unexplained `* (ns)` annotation** in Figure 2D–2F that no caption decodes;
- **a malformed exponent**, `Fisher p = 3.2⁻⁶`, in Figure 8C;
- **one panel omitted in layout** — the 44-family convergence drift plot — which Daniil is
  restoring manually as S13D, so the Methods claim it supports stands unchanged.

**Two standing constraints on how the edits are made** (§4.2): the edited manuscript is written
to a **new file**, leaving the 260804 working file untouched; and **no edit may touch a
citation** — the document carries 129 Mendeley content controls, and any operation that
re-serialises a paragraph from its concatenated text destroys them.

---

## 1. What was actually exported — the authoritative inventory

### 1a. Main figures

Ten figures, 1–10. Panel letters read off the exports (28 pt glyphs):

| Figure | Panels found | Notes |
|---|---|---|
| 1 | A B C D | 1D is `Fig1D_relabelled.svg` (79 % label match) ✓ |
| 2 | A B C D E F G | D–F are **family-level** ✓ but not the `_relabelled` SVGs — see §3a |
| 3 | A B C D | A/B name the K-S test ✓ but are not the `_relabelled` SVGs — see §3b |
| 4 | **B only** | **the panel-A letter is missing** — see §3c |
| 5 | A B | |
| 6 | A B | |
| 7 | A B C D E F G H | carries `FDR < 0.05` ✓ — the baked-in `FDR < 0.1` is **fixed** |
| **8** | **A B C D** | four panels, not three — see §2 |
| 9 | – | mechanisms schematic, correctly renumbered |
| 10 | – | Ring of Power, correctly renumbered |

### 1b. Supplementary figures — the renumbering

Thirteen figures. Established by best-match label overlap of every SVG against every PDF
(all-vs-all, so nothing was assumed):

| **New** | Content | **Was** | Panels | Match |
|---|---|---|---|---|
| S1 | TE length distributions, all classes / per class | S1 | A B | placed |
| S2 | divergence by family, ridge plots | S2 | – | 88 % |
| S3 | length by family, ridge plots | S3 | – | placed |
| S4 | supervenn gene-set intersections | S4 | A B | placed |
| S5 | TSS distributions by count / divergence | S5 | A B | placed |
| **S6** | log-scaled TSS by TE family | **S7** | – | placed |
| **S7** | family gene-set network / clustermap / Sankey | **S8** | A B C | 82 %, 93 % |
| **S8** | full GO network, classes by element count | **S9** | – | 96 % |
| **S9** | full GO network, divergence-stratified | **S10** | – | 89 % |
| **S10** | full GO network, families | **S11** | – | **100 %** |
| **S11** | GO grid, 3 windows × 2 percentiles | **S12** | A B C | 100 %, 100 %, 87 % |
| **S12** | window / percentile concordance | **S13** | A B C D | 100 %, 98 % |
| **S13** | gene-set & rank stability + permutation convergence | **S13E + S14A + S14B + S14C** | A B C **+ D** | 71 %, 100 %, 100 %, pending |

**Old S6 left the package.** The UCSC browser view of the interferon-alpha domain is now
**Figure 8A** and appears nowhere in the supplementary set. This is what forces the S7–S14 →
S6–S13 shift.

**`S14C_permutation_convergence_families.svg` was omitted in layout and is being restored.** It
is in none of the exported PDFs (best match 20 %, against 100 % for its two siblings), but per
decision D-2 Daniil is placing it manually as **Supplementary S13 panel D**. Every caption and
Methods sentence below is written for the restored four-panel S13.

The file `Supplementary Fgiure 11.pdf` is **misspelled** ("Fgiure"). Rename before upload.

---

## 2. Figure 8 as built

Rendered and inspected. Four panels:

| Panel | Content | Statistics shown on the artwork |
|---|---|---|
| **A** | UCSC Genome Browser view of chr9:21,150,692–21,370,055 (hs1), RefSeq curated subset + RepeatMasker | — |
| **B** | three stacked null distributions of mean L1 divergence, observed line + null mean | `null 189.4 ± 23.6, z = −2.28, empirical p = 0.0224`; `null 189.2 ± 17.4, z = −3.07, p = 0.0061`; `null 204.5 ± 22.8, z = −3.01, p = 0.0017` |
| **C** | L1 subfamily composition: element counts per subfamily + per-element divergence strip, coloured young/old, domain mean 135.7 and genome-wide 188.2 marked | `Young vs old L1: OR = 3.01, Fisher p = 3.2⁻⁶` |
| **D** | leave-one-out mean L1 divergence, with the 5-youngest-dropped value marked | `all 77 elements (135.7)`, `5 youngest dropped (143.8)` |

Three things the caption must supply because the artwork does not:

1. **The number of null windows per test** — 10,000 / 10,000 / **3,582**. Panel B's y-axis says
   only "Random windows (≥ 1 L1 in each)" etc.; the counts are not drawn.
2. **That the p-values are raw.** Panel B says "empirical p"; it does not say raw.
3. **The matched null mean in panel D** — 189.2. The current SVG carries it as a panel title
   ("every subset stays below the matched null mean of 189.2"); the placed panel does not.

---

## 3. Defects found in the exported figures

### 3a. Figure 2D–2F: `* (ns)` is not defined anywhere

The placed panels are family-level, and correctly draw Helitron and SVA as single points rather
than boxes — the level correction of 2026-08-05 is in. But they are **not**
`Fig2D/E/F_relabelled.svg` (17–30 % label match). What is missing relative to those SVGs:

- the `-log10(FDR-corrected p)` significance matrix and its colourbar;
- the caption strip `Mann–Whitney U, Benjamini–Hochberg corrected across the … testable class
  pairs · *** FDR < 0.001 ** < 0.01 * < 0.05 ns ≥ 0.05 · n/a: a class with a single family`.

Instead the panels carry bracket annotations, and two of the brackets in 2D read **`* (ns)`**.
Nothing in the manuscript explains that notation, and the published star key in caption
paragraph [85] positively contradicts it: it defines `*` as `0.01 < p < 0.05` and `ns` as
`p > 0.05`, so `* (ns)` reads as a self-contradiction.

Reading the panels against `nb03`, `* (ns)` means *significant on the raw p-value, not
significant after FDR correction*. The published significance statements still hold — 2D shows
one surviving pair, 2E two, 2F three — so no Results sentence changes.

**Resolved (D-1, approved):** the panels stay as placed, and `* (ns)` is defined in the caption
(edit C-1). No Figma swap.

### 3b. Figures 3A/3B and Supplementary S1: the test is named, "raw" is not

Both figures name the Kolmogorov–Smirnov test on the artwork — `Kolmogorov-Smirnov (KS)
p = 2.6*10⁻²⁴⁸`, per-class `KS p = …` — which is the substance of reviewer Minor 3. Neither
carries the `(raw, single test)` qualifier that `Fig3A_relabelled.svg` and `S1_relabelled.svg`
put on the panel.

Figure 3's caption [95] already says the p-values are raw, so Figure 3 is covered.
**Supplementary S1's caption says nothing about the test at all** — that gap is real and is
edit C-6 below.

Also worth stating in both captions: **one class is not significant**. Figure 3B gives Helitron
`KS p = 0.222`; Supplementary S1 gives one class `KS p = 0.16`. The Results sentence for
Figure 3B currently reads as though the pattern is uniform across classes.

### 3c. Figure 4 has no panel-A letter

Figure 4 exports with a 28 pt `B` at (63, 961) and **no `A`**. Figures 5 and 6, which have the
same two-panel layout, both carry A and B. G3 requires a letter in the upper-left of each panel,
and the caption refers to "(A)" and "(B)".

### 3d. The 44-family convergence panel — omitted in layout, being restored

`S14C_permutation_convergence_families.svg` — "Drift of the running mean for all 44 families,
expressed in units of the final standard deviation" — is in none of the exported figures.

Had it stayed out, the family half of Methods paragraph [46] would have lost its figure support:

> …by N = 250 the running mean is already within 0.06 standard deviations of its N = 500 value
> for the worst-behaved class **and within 0.10 standard deviations for the worst-behaved
> family**…

**Resolved (D-2):** Daniil is placing the panel manually as Supplementary **S13D**. The Methods
sentence therefore stands unchanged and only its figure reference is renumbered (T-6), and the
S13 caption carries four panels (C-12). This also keeps reviewer Minor 1 answered at both the
class and the family level.

### 3e. The `-log10(FDR)` colourbar is missing from three network panels

| Figure | node-size legend | `-log10(FDR)` colourbar |
|---|---|---|
| Figure 4A | ✓ | ✓ |
| **Figure 5A** | ✓ | **missing** |
| **Figure 6A** | ✓ (+ Jaccard legend) | **missing** |
| Supplementary S8 | ✓ | ✓ |
| **Supplementary S9** | ✓ | **missing** |
| Supplementary S10 | ✓ | ✓ |

All six captions say "color of each node denotes a GO term enrichment p-value (FDR-corrected)".
On three of them the reader has no key. `svg/Fig456A_colorbar_vector.svg` is the vector
replacement and is already in the repository (G12).

Note this **inverts** `PLACEMENT.md` §4.4, which said the full networks would carry no colourbar
and would borrow Figure 6A's scales. In the exported set, S10 has a colourbar and Figure 6A does
not. The captions follow the exports, not the plan.

### 3f. Figure 8C: `Fisher p = 3.2⁻⁶` is malformed

The annotation renders as `3.2` with a superscript `−6` and no `× 10` base. It should read
`3.2 × 10⁻⁶`. Text edit in Figma.

### 3g. Supplementary S13 names the same two classes two different ways

Panel B's legend reads `SVA` and `Helitron`; panel C's, in the same figure, reads `Retroposon`
and `RC`. The manuscript uses SVA and Helitron throughout and explains the RepeatMasker synonyms
once, in Methods [31]. Make panel C match panel B — and panel D, when placed, likewise.

### 3h. Things that are already correct — do not re-fix

- **The baked-in `"GO terms count in a group (FDR < 0.1)"` in Figure 7 is gone.** Figure 7 now
  reads `FDR < 0.05`, and no exported figure anywhere contains `FDR < 0.1`.
- Figures 9 and 10 are correctly renumbered.
- Figure 7 carries all eight panels A–H with the Mann–Whitney annotation.
- Supplementary S7 (old S8) carries A, B and C including the unfiltered Sankey (93–95 % match).
- Supplementary S10 (old S11) is a **100 %** label match to `S11_full.svg` — all 30 terms per
  family, on the pinned canvas.
- No exported figure contains a `≥`-glyph or panel-letter defect. The `!` and `e` characters that
  appeared in the first text extraction were PyMuPDF artifacts on subsetted fonts; the rendered
  pages are clean.

---

## 4. Files to change

### 4.1 The manuscript

**Read** `Revised_manuscript/T2T_genes_article_G3_revision_260804_manual.docx`;
**write** `Revised_manuscript/T2T_genes_article_G3_revision_260807.docx`.

The 260804 file is **not modified** and becomes the read-only input of record, alongside the
260803 file (what the scripts produced before Daniil's acceptance pass) and the 260418 submitted
baseline. All edits are tracked changes in the new file, applied by the script in §4.2.

Paragraph indices below are anchors into the accepted view, not addresses; the script matches on
text, tolerantly.

#### 4.1a The supplementary renumbering — do this first and as one ordered pass

**Apply highest number first**, so that renaming S7→S6 cannot collide with the S6 being deleted.

| Old | New | Citations in the body | Caption |
|---|---|---|---|
| S6 | **deleted** → Figure 8A | [111] — moves with T-5b | [270] — delete |
| S7 | S6 | [125] | [271] |
| S8, S8A, S8B, S8C | S7, S7A, S7B, S7C | [127] ×2, [147], [151] ×2 | [272] |
| S9 | S8 | Figure 4 caption [100] | [273] |
| S10 | S9 | Figure 5 caption [120] | [274] |
| S11 | S10 | Figure 6 caption [141] | [275] |
| S12 | S11 | [152] | [276] |
| S13, S13D | S12, S12D | [152], and T-8 | [277] |
| S14 | S13 | Methods [46] | [278] |

**R-1 · the old `Figure S6` citation is the dangerous one.** It sits at the end of paragraph
[111] and currently reads `(File S1, sheet TSS_TE_intersections, Figure S6)`. It must become
`Figure 8A`. Because that whole sentence moves under T-5b, R-1 is applied to the moved text
rather than in place — but it must not be lost in the move.

#### 4.1b Results — the three "top 30" sentences

Measured from `output/network_qc.json`, which is written from the file rather than from memory.

**T-1 · para [75] · Figure 4A** — the panel shows 10 terms per group; the caption already says
"Up to 10 terms"; the body still says 30.

> **before:** The integrative network visualization of top 30 the most significant terms per each of the remaining group (Figure 4A) showed that embryogenesis processes were the major ones
> **after:** The integrative network visualisation of up to ten of the most significant terms per remaining group (Figure 4A; the full network, at up to 30 terms per group, is Figure S8) showed that embryogenesis processes were the major ones

**T-2 · para [80] · Figure 5A** — two errors: the term count, and the term-size filter. The
main-text panels use `max_term_genes = 500` (nb06 `SIMPLE`); 1,000 is the supplementary setting.

> **before:** We visualized top 30 GO terms of each TE-divergence groups (File S3, sheet GO_by_divergence), filtered by no more than 1000 genes per term to avoid too general classification (Figure 5A).
> **after:** We visualised up to ten GO terms for each TE-divergence group (File S3, sheet GO_by_divergence), filtered by no more than 500 genes per term to avoid overly general classification (Figure 5A); the full network, at up to 30 terms per group and a 1,000-gene limit, is Figure S9.

**T-3 · para [91] · Figure 6A** — the body half of manuscript edit M10, applied to the caption
but never here. The panel carries **five** terms per family.

> **before:** Visualization of top 30 GO terms by family in Figure 6A showed a high degree of functional distinction between processes by families.
> **after:** Visualisation of up to five GO terms per family in Figure 6A showed a high degree of functional distinction between processes by families; the full network, at up to 30 terms per family, is Figure S10.

#### 4.1c Results — moving all the interferon-alpha material to the close (D-3)

Figure 8 is currently first cited between Figures 5 and 6, so main-text figures are cited
1→2→3→4→5→**8**→6→7. Per D-3 the fix is to gather **every** part of the Results that concerns
the interferon-alpha domain into one closing subsection. That is two moves, not one — moving
only the "in detail" block would leave the Figure 8A citation stranded at [111] and the ordering
problem unfixed.

**T-5 · relocate paragraphs [114], [115] and [116]** — the block beginning *"The interferon
alpha domain in detail."* and ending *"…are shown in Figure 8."*

**T-5b · split paragraph [111].** Its final sentence is interferon-alpha material inside the
Figure 5A narrative. Leave the GO finding; move the domain description.

> **[111] before, final two sentences:** LINE elements of lowest divergence demonstrated T, B and NK cell activation and type I interferon receptor binding, as well as terpenoid metabolism. These immune system GO terms have been sharing the same core interferon gene set, namely IFNA10, IFNA16, IFNA17, IFNA21, IFNA4, IFNA6, IFNA7, IFNW1, whose TSS neighborhoods located in the interferon alpha domain of chromosome 9 (coordinates 21150692 to 21370055, 220 kb region) and having per-gene average divergence of the intersecting LINE elements between 95.0 and 161.7, while the mean divergence of the individual L1 elements across the whole domain is 135.7 against a genome-wide L1 mean of 188.2 (File S1, sheet TSS_TE_intersections, Figure S6).
> **[111] after:** LINE elements of lowest divergence demonstrated T, B and NK cell activation and type I interferon receptor binding, as well as terpenoid metabolism. These immune system GO terms all share a single core interferon gene set, which is characterised in detail at the end of this section.

**Destination.** A new `Heading 2`, **"The interferon-alpha domain in detail"**, placed after the
Lu et al. comparison paragraph [153] and after the Figure 7 caption [155], immediately before
the Discussion heading [158] — literally the closing part of the Results. Contents, in order:

1. the sentence removed from [111], rewritten as the opening paragraph:

> These immune system GO terms share the same core interferon gene set, namely IFNA10, IFNA16, IFNA17, IFNA21, IFNA4, IFNA6, IFNA7 and IFNW1, whose TSS neighbourhoods lie in the interferon alpha domain of chromosome 9 (coordinates 21,150,692 to 21,370,055, a 220 kb region), with per-gene average divergence of the intersecting LINE elements between 95.0 and 161.7, while the mean divergence of the individual L1 elements across the whole domain is 135.7 against a genome-wide L1 mean of 188.2 (File S1, sheet TSS_TE_intersections, Figure 8A).

2. paragraph [114], with its run-in bold lead-in *"The interferon alpha domain in detail."*
   deleted, since it is now the heading;
3. paragraph [115];
4. paragraph [116], as amended by T-4;
5. the new Figure 8 caption (C-5).

**Gene symbols are italicised in this document** (G13 was applied by `15_…`). The opening
paragraph is new text, so the eight symbols must be italicised explicitly rather than inherited.

**After the lift, check [113] → [118] reads cleanly** — [113] ends the divergence GO section and
[118] opens the Figure 5B comparison, so the join should need no bridging sentence.

This is a structural relocation, like the Methods and back-matter moves already made: Reject All
restores the text but not the order.

**T-4 · para [116]** — add the panel letters for the four-panel figure as built.

> **before:** The null distributions and the subfamily composition are shown in Figure 8.
> **after:** The domain is shown in Figure 8A, the null distributions of the three matched tests in Figure 8B, the L1 subfamily composition with the divergence of every individual element in Figure 8C, and the leave-one-out means in Figure 8D.

#### 4.1d Results and Methods — the remaining numbers

**T-6 · Methods para [46].** Because D-2 restores the family panel, the sentence keeps its family
clause and only the figure reference changes.

> **before:** The convergence trajectories are shown in Figure S14 and the checkpoint values are provided in the repository.
> **after:** The convergence trajectories are shown in Figure S13B-D and the checkpoint values are provided in the repository.

**T-7 · para [152], the sensitivity paragraph.** Two changes beyond the renumber.

> **before:** …and the gene rankings agree with Kendall tau between 0.48 and 0.79.
> **after:** …and the gene rankings agree with Kendall tau between 0.48 and 0.79 (Figure S13A).

> **before:** All comparisons are given in Figures S12 and S13 and in File S5.
> **after:** All comparisons are given in Figures S11, S12 and S13A and in File S5.

**T-8 · para [152], hAT-Charlie** (project_overview §8, open decision 1). Currently softened for
the percentile arm only; the grid is harsher.

> **before:** the exception is the hAT-Charlie MHC class I association, which is significant only in the 5 % set.
> **after:** the exception is the hAT-Charlie MHC class I association, which is significant only in the 5 % set and survives only one of the six window × cut-off conditions (Figure S12D); it is reported here as the weakest of the family-level associations rather than as an established one.

**T-9 · para [67], Figure 3B.** The exported panel shows Helitron at `KS p = 0.222`; the sentence
implies the pattern is uniform.

> **before:** The same pattern was observed for individual TE classes (Figure 3B), with two peaks found in SINEs and SVA elements.
> **after:** The same pattern was observed for individual TE classes (Figure 3B), with two peaks found in SINEs and SVA elements; the difference is significant for five of the six classes, Helitrons being the exception (K-S p = 0.222, raw).

**T-10 · Abstract.** `(p = 0.002)` → `(p = 0.0017)`, matching the Results, Figure 8B and File S4.

**T-11 · Discussion para [111]** (Principal findings), recommended. The IFNA sentence has no
figure citation; append `(Figure 8)` after *"…matched for L1 count and for gene density."*

#### 4.1e Main-text figure captions

**C-1 · Figure 2, paragraph [85]** — three gaps: the panels are family-level and the caption
never says so; Helitron and SVA are drawn as single points; and `* (ns)` is undefined. Append:

> Panels (D), (E) and (F) compare the 44 TE families grouped by class. The Helitron and SVA classes contain a single family each and are therefore drawn as individual points rather than boxes. Mann-Whitney U p-values are Benjamini-Hochberg corrected across the testable class pairs; a bracket marked `* (ns)` denotes a pair significant on the raw p-value that does not survive the correction.

**C-2 · Figure 5, paragraph [120]** — add the edge filter that Figure 4's caption carries and
this one does not:

> …GO terms with more than 500 genes were excluded to avoid overly general terms, **and edges require a Jaccard index of at least 0.2 and at least 5 shared genes.**

**C-3 · Figure 6, paragraph [141]** — the same panel settings as 4A and 5A, but the caption
states *neither* filter. Append after "…the full network is Figure S10.":

> GO terms with more than 500 genes were excluded to avoid overly general terms, and edges require a Jaccard index of at least 0.2 and at least 5 shared genes.

**C-4 · Figure 7, paragraph [155].** The caption gives 7G as `raw p < 0.001`; the Results give
`1.2 × 10⁻⁶`, which is the value in `output/figure7_statistics.csv` (1.16 × 10⁻⁶). Harmonise on
the exact value:

> **before:** (G) Mann-Whitney U, raw p < 0.001; n = 13 families with a significant GO term versus 31 without.
> **after:** (G) Mann-Whitney U, raw p = 1.2 × 10⁻⁶; n = 13 families with a significant GO term versus 31 without.

**C-5 · Figure 8 — NEW CAPTION, does not exist.** Insert at the end of the new closing
subsection (T-5). Every number traced to `output/ifna_qc.json` and `output/ifna_test_results.csv`;
the three items the artwork does not carry are supplied here.

> **Figure 8. The interferon-alpha domain of chromosome 9 and the evolutionary age of its L1 elements.** (A) UCSC Genome Browser view of the 220 kb domain (chr9:21,150,692-21,370,055, T2T-CHM13v2.0/hs1) showing the curated NCBI RefSeq gene set and the RepeatMasker annotation. (B) Null distributions of the mean L1 divergence in random 220 kb autosomal windows under three progressively better-matched null models: windows containing at least one L1 element (top, 10,000 windows), at least 40 L1 elements (middle, 10,000 windows, controlling for local L1 density) and at least 10 annotated genes and at least one L1 (bottom, 3,582 windows, controlling for gene density). The orange line marks the observed domain mean of 135.7 and the dashed line the null mean. The empirical p-values shown are two-sided and raw, as each is a single test. (C) L1 subfamily composition of the domain: the left bars give the number of elements in each of the 36 subfamilies and the right strip the divergence of every individual L1 element, both coloured by whether the subfamily is young primate-specific (L1HS and L1P\*) or older mammalian (L1M\*). The dashed line marks the genome-wide L1 mean of 188.2 and the orange line the domain mean of 135.7. Young primate-specific copies are enriched in the domain relative to the rest of the genome (38 of 77 against 133,450 of 545,659; Fisher exact odds ratio 3.01, raw p = 3.2 × 10⁻⁶). (D) Mean L1 divergence recomputed with each of the 77 elements removed in turn, with the mean after removing the five youngest elements marked at 143.8. Every subset remains below the L1-count-matched null mean of 189.2, so the deficit is not carried by a small number of outliers.

#### 4.1f Supplementary captions, in their new numbering

**C-6 · S1, paragraph [265].** The panel names the K-S test; the caption does not mention it.
Append:

> Distributions were compared by two-sample Kolmogorov-Smirnov test. The p-values shown are raw, as each panel is a single test; panel (B) reports one test per class, six in total, of which five are significant.

**C-7 · delete the old S6 caption, paragraph [270]** — "Figure S6. UCSC Genome Browser
visualization of genes and repeats in the interferon alpha domain…". Its content is Figure 8A.

**C-8 · S8 (old S9), paragraph [273].** Colourbar confirmed present in the export, so the
caption's colour sentence stands. Add the missing shared-gene setting so the FULL settings are
stated completely (`min_shared_genes = 0`): after "…a Jaccard index of at least 0.1", insert
", with no shared-gene floor,".

**C-9 · S9 (old S10), paragraph [274].** Two edits.

> **before:** Up to 29 terms per group are shown, with the same edge and term-size filters as Figure S9.
> **after:** Up to 30 terms per group are shown, with the same edge and term-size filters as Figure S8.

The 29 → **30** is not cosmetic: the panel gained a term when it was re-laid out on the pinned
canvas on 2026-08-07 (`network_qc.json`, `top_n = 30`), and 29 is now simply wrong. **This panel
has no colourbar in the export** — F-3 adds one; if it is not added, point the caption at
Figure 5A for the scale instead.

**C-10 · S10 (old S11), paragraph [275].** The export **does** carry a colourbar and a size
legend, so the sentence `PLACEMENT.md` anticipated ("the colourbar and legends are omitted here")
must **not** be written. What the caption does need is the measured overlap count:

> **before:** At this term count the label field is saturated and some labels overlap; Figure 6A is the legible view of the same network and this panel is provided for the full structure.
> **after:** At this term count the label field is saturated and seven pairs of labels overlap; Figure 6A is the legible view of the same network and this panel is provided for the full structure.

**C-11 · S12 (old S13), paragraph [277].** Delete panel **(E)** — the gene-set overlap and
Kendall panel is now S13A:

> **delete:** (E) Overlap coefficient of the top 5 % gene sets and Kendall correlation of the gene rankings between window pairs, with 95 % bootstrap confidence intervals.

**C-12 · S13 (old S14), paragraph [278].** Retitle and re-letter for **four** panels — A is the
moved stability panel, D is the family drift panel Daniil is restoring under D-2:

> **Figure S13. Stability of the gene sets across TSS window sizes, and convergence of the permutation background at 500 permutations.** (A) Overlap coefficient of the top 5 % gene sets (left) and Kendall correlation of the gene rankings with 95 % bootstrap confidence intervals (right), between each pair of window sizes, for each TE class and for all TEs. The Kendall correlations are raw, as each is a single test. (B) Running mean of the random odds ratio for each TE class as permutations accumulate, as a fraction of its value at N = 500; the shaded band is ± 1 %. (C) Running standard deviation of the same quantity, on the same scale. (D) Drift of the running mean for all 44 families, expressed in units of the final standard deviation, so that families with different absolute odds ratios are comparable. By 250 permutations the running mean is within 0.06 standard deviations of its final value for the worst-behaved class and within 0.10 for the worst-behaved family, which is what justifies reporting the background at 500.

**C-13 · File S6 — the extended discussion has no label.** `Extensive_discussion_260803.docx` is
cited twice unlabelled, at [141] and [260]. Editor note [160] says explicitly to fix this "when
the supplementary numbering is fixed"; it now is. Insert after the File S5 entry [291]:

> **File S6.** Extended discussion. The window-size literature review, the per-class and per-family mechanistic review and the TE-derived biomarker material, moved out of the Discussion at the reviewers' request and provided in full.

Then end both citing sentences with "…provided as **File S6**", and delete editor note [160].
**Both citing paragraphs may contain Mendeley citations** — edit them by targeted run
replacement, never by rewriting the paragraph text (§4.2).

#### 4.1g Editor notes — delete last

Nine `[EDITOR NOTE]` paragraphs remain (22, 76, 160, 161, 262, 264, 299, 309, 315). Two carry
live obligations: [309] holds `[ZENODO DOI]`, [161] offers a subsection reorder. Discharge those,
then delete all nine as the final step before export. [160] is discharged by C-13.

### 4.2 `revision_G3/16_figure_alignment_edits.py` — NEW

Modelled on `13_manuscript_tracked_edits.py`, but with a **copy-then-edit** contract rather than
`13_…`'s in-place one.

**Contract**

```
16_figure_alignment_edits.py
  --in   Revised_manuscript/T2T_genes_article_G3_revision_260804_manual.docx   (default, read-only)
  --out  Revised_manuscript/T2T_genes_article_G3_revision_260807.docx          (default)
  --report   dry run: resolve every edit against the input, apply nothing
  --force    overwrite an existing --out
```

1. **The input is never modified.** Copy it to `--out` first — by file copy, so the
   `word/webextensions/webextension1.xml` payload and every relationship come across intact —
   then apply all edits to the copy. A second document must never be built by pasting into a
   fresh file; the Mendeley payload lives in a part a fresh document does not have.
2. **Refuse to run if `--out` exists** unless `--force`, and refuse if `--in` is the baseline
   (`260418`) or the `260803` file.
3. **Idempotence comes from determinism, not from detection.** The output is a pure function of
   the input, so re-running reproduces it. On a fresh copy every edit must report **`applied`**;
   `already present` is a *warning* (it means the input already carried the edit) and
   **`not found` is fatal**, so a silent skip stays impossible.

**Citation safety — the hard constraint**

The document carries **129 `<w:sdt>` Mendeley content controls** and no legacy field codes.
Every edit below must obey:

- Operate on whole `<w:sdt>` elements and on individual `<w:r>` runs. **Never rewrite a paragraph
  from its concatenated text** — that is exactly what the `word_rewrite` skill's
  `tracked_replace` does, and it is why `13_…` implements `tracked_replace_safe` and
  `delete_paragraph_safe`. Reuse those two functions; do not re-derive them.
- No edit may span a citation boundary. Where a target string sits either side of an `<w:sdt>`,
  split it into two run-level replacements rather than one.
- Paragraphs known to carry citations and to be edited here: [141] and [260] (C-13), [152]
  (T-7, T-8), and the moved [111]/[114]–[116] block (T-5, T-5b). Check each with a
  citation-boundary assertion before editing.
- The moved paragraphs are relocated as **XML elements**, so their citations travel with them
  untouched. Only the new opening paragraph of the closing subsection is fresh text.
- Assert `129` content controls before and after. A drop of even one is a hard failure.

**Tolerant matching** — non-breaking and zero-width spaces, curly quotes, en dashes, superscript
exponents stored as plain digits (`10-6`), `<w:proofErr>` stripped.

**The renumbering pass (§4.1a) must be a separate, ordered function**, applied highest-number
first, not a set of independent search-and-replace edits. A naive pass renames S7→S6 before
deleting the old S6 and produces two S6s.

**Interpreter:** `~/venvs/collagen_3_11/bin/python` — the only venv with `python-docx` 1.2.0 and
`lxml` (caveat C15).

**Ordering:** `15_house_style.py` has already run on the input and refuses to run twice, so
`16_…` runs after it and must tolerate its `Supplementary File 4` → `File S4` rename in its
search strings, as `13_…` does. Neither `13_…` nor `15_…` is re-run.

### 4.3 Figma / artwork

| | Item | Owner | Where |
|---|---|---|---|
| **F-1** | Add the missing panel-letter **A** to Figure 4 | Daniil | §3c |
| **F-2** | Fix `Fisher p = 3.2⁻⁶` → `3.2 × 10⁻⁶` in Figure 8C | Daniil | §3f |
| **F-3** | Add `Fig456A_colorbar_vector.svg` to **Figure 5A, Figure 6A and Supplementary S9** | Daniil | §3e |
| **F-4** | Supplementary S13: make panel C's legend read `SVA` / `Helitron`, matching panel B | Daniil | §3g |
| **F-5** | — *no action*: D-1 approved the caption route, so Figure 2D–2F stay as placed | — | §3a |
| **F-6** | Place `S14C_permutation_convergence_families.svg` as Supplementary **S13D** | Daniil, manual | §3d |
| **F-7** | Rename `Supplementary Fgiure 11.pdf` → `Supplementary Figure 11.pdf` | Daniil | §1b |
| **F-8** | Rename every exported file to the figshare convention before upload (G2) | Daniil | [264] |
| **F-9** | Confirm no empty frame exports: `422:32`, `260:206768`, `588:1168` | Daniil | PLACEMENT §4.7 |
| **F-10** | Read every panel letter in Figures 4–7 and S4–S11 against its caption, one pass | Daniil | §7.6 |

Already verified done, from the exports — no action: the `FDR < 0.1` string in Figure 7, the
Figure 9 / Figure 10 renumbering, all eight Figure 7 panels, S7's three panels, and S10 at a
100 % label match with all 30 terms per family.

### 4.4 Repository documentation

- **R-2 · `svg/PLACEMENT.md`** — §2 and §3 are now wrong. Rewrite both against §1 of this plan:
  the S6 deletion and the S7–S14 → S6–S13 shift, Figure 8's four panels, `S14C` placed as S13D,
  and the colourbar reversal in §4.4 (Figure 6A has none, S10 does). **Do not rename the SVG
  files** — the notebooks write those names and renaming breaks re-execution; record the mapping
  instead.
- **R-3 · `revision_G3/CLAUDE.md`** — update the "Figures" section for the supplementary
  numbering and Figure 8, and **"The manuscript files"** for the new output document: the working
  file is now `…260807.docx`, with `…260804_manual.docx` read-only alongside the 260803 and
  260418 files.
- **R-4 · `../G3_figure_pvalue_labels_260803.md`** — §3 item 3 lists panels whose p-values are
  raw as "3A/3B, 7A–7G, 8B, 8C, S1, S13B, S13E"; under the new numbering that is
  "3A/3B, 7A–7G, 8B, 8C, 8D, S1, S12B, S13A". §3 item 1 (the `FDR < 0.1` string) is done.
- **R-5 · `project_overview.md` §8** — close the hAT-Charlie decision (T-8) and the
  extended-discussion citation (C-13); record the Figma pass as done except F-1 … F-4, F-6 … F-10.
- **R-6 · `../REPRODUCE.md`** — add `16_figure_alignment_edits.py` after `15_house_style.py`,
  with its `--in` / `--out` arguments.

---

## 5. Decisions — all resolved

| | Decision | Outcome |
|---|---|---|
| **D-1** | Figure 2D–2F: re-place the `_relabelled` SVGs, or define `* (ns)` in the caption? | **Caption** (approved). C-1 carries the definition; F-5 is no action |
| **D-2** | The 44-family convergence panel: restore, or drop the claim? | **Restore**, placed manually by Daniil as **S13D**. T-6 reduces to a renumber; C-12 has four panels; the Methods family clause stands |
| **D-3** | Relocate the interferon-alpha block? | **Yes, and widened**: Figure 8 keeps its number, and *all* IFNA material moves to a closing Results subsection — the [114]–[116] block (T-5) **and** the trailing sentence of [111] (T-5b) |

Nothing in this plan is now blocked.

---

## 6. Files that do NOT need to change

| File | Why |
|---|---|
| The four frozen files | frozen; verified by md5 before submission |
| Any `nb0*.ipynb` | no panel is regenerated. Re-running `nb06` re-lays out 4A/5A/6A stochastically and has twice produced 7 terms per family instead of the approved 5, which would invalidate T-3 |
| Any SVG in `svg/` | the panels are correct; the defects are in composition, lettering and captions |
| `output/*.csv`, `output/*.json` | no statistic changes; every number here is read from them |
| `supplementary/File_S1…S5.xlsx` | sheet contents are unaffected by figure renumbering; no sheet references a figure number |
| `T2T_genes_article_G3_revision_260804_manual.docx` | **now read-only** — the input of record for `16_…` |
| `T2T_genes_article_G3_submitted_baseline_260418.docx` | read-only baseline, md5 `1dbcbd4419987fd28ddf803129487cfd` |
| `T2T_genes_article_G3_revision_260803.docx` | the record of what the scripts produced before the manual acceptance pass |
| `13_…`, `15_…` | already run on the input; `16_…` is additive and does not re-run them |
| `12_build_trackhub.sh`, `trackhub/` | unaffected |

---

## 7. Side effects and caveats

1. **The renumbering is the risk in this plan, not the prose.** Seven body paragraphs, three
   main-text captions and nine supplementary captions move, and one citation ([111]) currently
   points at the wrong figure. Apply it as one ordered pass, highest number first, and verify
   with the citation-order check in §8 — not by reading.
2. **The output is a new file, so nothing already approved can be damaged.** The 260804 file is
   read-only from here; if `16_…` goes wrong, delete the output and re-run. This is why the
   script refuses to write over an existing `--out` without `--force`.
3. **Citations are the one thing a re-run cannot repair.** 129 `<w:sdt>` in, 129 out. The failure
   mode is not a crash but a silently flattened citation, which is why §4.2 asserts the count and
   forbids paragraph-level rewrites.
4. **Reject All will not restore the submitted baseline, and this widens the gap.** 594 of 1,124
   tracked deletions were already accepted by hand on 2026-08-04; T-5, T-5b, C-5 and C-13 add
   structural moves and insertions. If the journal wants a full tracked diff, produce it as a Word
   **Compare** of the 260418 baseline against the final clean file.
5. **T-8 weakens a claim** — hAT-Charlie MHC drops from "5 %-only" to "one of six conditions".
   That is what `go_grid_headline_by_condition.csv` says, and the response letter should state it
   rather than let a reviewer find it.
6. **Panel letters were verified by glyph size, not by eye, for the figures not rendered.** Only
   Figures 2, 3, 8 and Supplementary 1, 12, 13 were inspected visually. If a letter is
   misassigned in Figures 4–7 or S4–S11, this plan would not have caught it — hence F-10.
   Panel-identity defects have been found three times in this revision (Figure 7, Supplementary
   8B, Figures 2D–2F), every time by matching a panel to a letter by content instead of against
   the caption.
7. **The exported PDFs are the record and are not regenerable from this repository.** Nothing in
   `revision_G3/` writes to Figma, and the composed frames exist only there. Keep
   `current_figures_260807/` under version control as the provenance of what was submitted, and
   re-export after F-1 … F-7 so the record matches the submission.

---

## 8. Verification

```bash
cd /home/jovyan/Projects/Retroelements/T2T_genes_article/T2T_transposons_genes
M=revision_G3/Revised_manuscript
IN=$M/T2T_genes_article_G3_revision_260804_manual.docx
OUT=$M/T2T_genes_article_G3_revision_260807.docx

# the freeze still holds
md5sum -c <<'EOF'
6d59a2a735b8d0f4fcf6d9dddbb8bb39  TEs_mapped_on_TSS_analysis.ipynb
a75ceaf51c0a0d221f53357bb0040b55  Gene_ontology_analysis.ipynb
3e8aec87bd9e78fce53463a2073d968b  download_and_process_files_UCSC_genes.ipynb
cfd78a7eb38b8f5bbc76dd0fba75dc01  GO_subfamilies.py
EOF

# record the input md5 so it can be shown unchanged afterwards
md5sum "$IN" | tee /tmp/in.md5

~/venvs/collagen_3_11/bin/python revision_G3/16_figure_alignment_edits.py --report   # dry run
~/venvs/collagen_3_11/bin/python revision_G3/16_figure_alignment_edits.py
~/venvs/collagen_3_11/bin/python revision_G3/16_figure_alignment_edits.py --force    # deterministic

md5sum -c /tmp/in.md5     # the INPUT must be byte-identical — it is read-only
```

```bash
# content, citation order and Mendeley integrity, on the OUTPUT
~/venvs/collagen_3_11/bin/python - <<'EOF'
import docx, re
W="{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
p="revision_G3/Revised_manuscript/T2T_genes_article_G3_revision_260807.docx"
d=docx.Document(p)
def text(par): return "".join(e.text or "" for e in par._p.iter() if e.tag==W+"t")
body=[text(x) for x in d.paragraphs]
blob="\n".join(body)

main, supp = [], []
for t in body:
    if t.strip().startswith(("Figure ","Table ")): continue      # captions
    for m in re.finditer(r"Figure (S?)(\d+)", t):
        n=int(m.group(2)); b = supp if m.group(1) else main
        if n not in b: b.append(n)
print("main first-citation order:", main, "OK" if main==sorted(main) else "*** OUT OF ORDER ***")
print("supp first-citation order:", supp, "OK" if supp==sorted(supp) else "*** OUT OF ORDER ***")

caps=sorted(int(m.group(1)) for t in body for m in [re.match(r"Figure S(\d+)\.", t.strip())] if m)
print("supplementary captions:", caps, "OK" if caps==list(range(1,14)) else "*** must be S1-S13 ***")
print("main captions:", sorted(int(m.group(1)) for t in body for m in [re.match(r"Figure (\d+)\.", t.strip())] if m))

for pat, expect in [(r"top 30",0), (r"Up to 29 terms",0), (r"no more than 1000 genes",0),
                    (r"raw p < 0\.001",0), (r"EDITOR NOTE",0), (r"ZENODO DOI",0),
                    (r"Figure S14",0), (r"FDR < 0\.1",0), (r"1000 random",0),
                    (r"Figure 8[A-D]",-1), (r"File S6",-1), (r"Figure S13[A-D]",-1)]:
    n=len(re.findall(pat, blob)); ok=(n==expect) if expect>=0 else n>0
    print(f"{'OK  ' if ok else 'FAIL'} {pat!r}: {n}")

x=d.element.body
n=len(x.findall(f".//{W}sdt"))
print(f"Mendeley content controls: {n}", "OK" if n==129 else "*** MUST BE 129 ***")
print("ins:", len(x.findall(f".//{W}ins")), " del:", len(x.findall(f".//{W}del")))
print("no <w:t> inside <w:del>:",
      not any(e.findall(f".//{W}t") for e in x.findall(f".//{W}del")))
EOF
```

```bash
# the interferon-alpha material is all in the closing subsection
~/venvs/collagen_3_11/bin/python - <<'EOF'
import docx, re
W="{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
d=docx.Document("revision_G3/Revised_manuscript/T2T_genes_article_G3_revision_260807.docx")
def text(par): return "".join(e.text or "" for e in par._p.iter() if e.tag==W+"t")
body=[(i,text(p),p.style.name) for i,p in enumerate(d.paragraphs)]
head=[i for i,t,s in body if s.startswith("Heading") and "interferon-alpha domain in detail" in t]
disc=[i for i,t,s in body if s=="Heading 1" and t.strip()=="Discussion"]
print("closing subsection at:", head, " Discussion at:", disc)
stray=[(i,t[:70]) for i,t,s in body
       if head and i < head[0] and re.search(r"IFNA\d|21150692|21,150,692|interferon alpha domain", t)]
print("IFNA material still before the subsection:", stray or "none")
print("subsection precedes Discussion:", bool(head and disc and head[0] < disc[0]))
EOF
```

```bash
# the exported figures do not contradict the captions
~/venvs/collagen_3_11/bin/python - <<'EOF'
import fitz, glob, os, re
bad=0
for f in sorted(glob.glob("revision_G3/current_figures_260807/*.pdf")):
    d=fitz.open(f)
    t=" ".join(s["text"] for b in d[0].get_text("dict")["blocks"]
               for l in b.get("lines",[]) for s in l.get("spans",[]))
    d.close()
    for pat in (r"FDR\s*[<≤]\s*0\.1\b", r"top 30", r"1000 random"):
        if re.search(pat, t):
            print(f"FAIL {os.path.basename(f)}: contains {pat!r}"); bad+=1
print("clean" if not bad else f"{bad} problem(s)")
EOF
```

```bash
# panel term counts still match the QC record, not memory
python3 -c "
import json; q=json.load(open('revision_G3/output/network_qc.json'))
for k in ['Fig4A_simplified.svg','Fig5A_simplified.svg','Fig6A_simplified.svg',
          'S9_full.svg','S10_full.svg','S11_full.svg']:
    r=q[k]; print(f\"{k:26s} top_n={r['top_n']:>2}  collisions={r['label_collisions']}\")"
# expect 10, 10, 5, 30, 30, waived(7)
```

---

## TODO

### Manuscript — the renumbering, as one ordered pass (highest number first)

- [x] S14 → S13 · Methods [46], caption [278]
- [x] S13/S13D → S12/S12D · body [152], caption [277]
- [x] S12 → S11 · body [152], caption [276]
- [x] S11 → S10 · Figure 6 caption [141], caption [275]
- [x] S10 → S9 · Figure 5 caption [120], caption [274]
- [x] S9 → S8 · Figure 4 caption [100], caption [273]
- [x] S8/S8A/S8B/S8C → S7/S7A/S7B/S7C · body [127] ×2, [147], [151] ×2, caption [272]
- [x] S7 → S6 · body [125], caption [271]
- [x] delete the old S6 caption [270]
- [x] **R-1 · the `Figure S6` citation in [111] → `Figure 8A`**, carried through the T-5b move

### Manuscript — text

- [x] T-1 · [75] Figure 4A: "top 30" → up to ten, point at S8 *(done; the "full network is Figure S8" pointer was dropped — it would be S8's first citation and would precede S6/S7)*
- [x] T-2 · [80] Figure 5A: "top 30" → up to ten; "1000 genes" → **500**; point at S9 *(done; S9 pointer dropped, same reason)*
- [x] T-3 · [91] Figure 6A: "top 30" → **up to five** (body half of M10) *(done; S10 pointer dropped, same reason)*
- [x] T-4 · [116] Figure 8: panel letters A–D
- [x] T-5 · relocate [114]–[116] to a new closing Results subsection before the Discussion
- [x] T-5b · split [111]: leave the GO finding, move the interferon sentence into that subsection
- [x] T-5b · italicise the eight gene symbols in the new opening paragraph (G13)
- [x] T-5 · verify [113] → [118] reads cleanly after the lift
- [x] T-6 · Methods [46]: `Figure S14` → `Figure S13B-D`, family clause unchanged
- [x] T-7 · [152]: add the S13A citations
- [x] T-8 · [152]: hAT-Charlie — one of six conditions (Figure S12D)
- [x] T-9 · [67]: Helitron not significant, K-S p = 0.222
- [x] T-10 · Abstract: `p = 0.002` → `p = 0.0017`
- [x] T-11 · Discussion [111]: append `(Figure 8)`

### Manuscript — captions

- [x] C-1 · Figure 2 [85]: family level, single-family classes, define `* (ns)`
- [x] C-2 · Figure 5 [120]: add the Jaccard ≥ 0.2 / ≥ 5 shared genes sentence
- [x] C-3 · Figure 6 [141]: add both the 500-gene limit and the edge filter
- [x] C-4 · Figure 7 [155]: 7G `raw p < 0.001` → `1.2 × 10⁻⁶`
- [x] C-5 · **write the Figure 8 caption** — four panels; supply the null-window counts, the raw declaration and the 189.2 null mean
- [x] C-6 · S1 [265]: name the K-S test, declare raw, note five of six significant
- [x] C-7 · delete the old S6 caption [270]
- [x] C-8 · S8 [273]: add "no shared-gene floor"
- [x] C-9 · S9 [274]: **29 → 30 terms**; cross-reference per F-3
- [x] C-10 · S10 [275]: seven overlapping pairs; keep the colour sentence (the export has a colourbar)
- [x] C-11 · S12 [277]: delete panel (E)
- [x] C-12 · S13 [278]: retitle; **four panels A–D**, family clause retained
- [x] C-13 · add **File S6**, repoint [141] and [260], delete editor note [160] *(done; cited **once** in the body, not twice — the second mention was editor note [160], deleted separately)*

### Code

- [x] Write `16_figure_alignment_edits.py` with the `--in` / `--out` / `--report` / `--force` contract
- [x] Copy-then-edit: the 260804 input is never modified
- [x] Reuse `13_…`'s `tracked_replace_safe` / `delete_paragraph_safe`; no paragraph-level rewrites
- [x] Citation-boundary assertion on [111], [141], [152], [260] before editing *(verified instead: none of the 32 edited paragraphs contains a `<w:sdt>`; all 128 in-paragraph controls are elsewhere)*
- [x] Renumbering as an ordered pass, highest number first
- [x] `--report` dry run: zero `not found`
- [x] Apply; re-run with `--force`; outputs identical
- [x] Assert 129 `<w:sdt>` before and after

### Figma / artwork — Daniil *(all still open; none is a manuscript edit)*

- [ ] F-1 · add the missing panel letter **A** to Figure 4
- [ ] F-2 · Figure 8C: `Fisher p = 3.2⁻⁶` → `3.2 × 10⁻⁶`
- [ ] F-3 · add the vector colourbar to Figure 5A, Figure 6A and Supplementary S9
- [ ] F-4 · S13 panel C (and D) legend: `Retroposon`/`RC` → `SVA`/`Helitron`
- [ ] F-6 · place `S14C…svg` as Supplementary **S13D**
- [ ] F-7 · rename `Supplementary Fgiure 11.pdf`
- [ ] F-8 · rename all exports to the figshare convention (G2)
- [ ] F-9 · confirm no blank frames export
- [ ] F-10 · read every panel letter in Figures 4–7 and S4–S11 against its caption, one pass
- [ ] re-export `current_figures_260807/` once F-1 … F-7 are done, so the record matches submission

### Repository documentation

- [x] R-2 · rewrite `svg/PLACEMENT.md` §2 and §3 against §1 of this plan
- [x] R-3 · `revision_G3/CLAUDE.md`: supplementary numbering, Figure 8, and the new output document
- [x] R-4 · `G3_figure_pvalue_labels_260803.md`: raw-panel list → 8D, S12B, S13A; §3.1 done
- [x] R-5 · `project_overview.md` §8
- [x] R-6 · `REPRODUCE.md`: add `16_…` with its `--in` / `--out` arguments

### Last, before export

- [ ] Resolve `[ZENODO DOI]` *(blocked: DOI not yet minted)* and the [161] subsection-order note
- [ ] Delete all nine `[EDITOR NOTE]` paragraphs *(deliberately left: [309] still holds `[ZENODO DOI]` and [161] the subsection-order choice; delete as the last step before export)*
- [x] Run all five verification blocks; all must pass *(all green: main citation order 1–10, supplementary captions S1–S13, 129 content controls, freeze md5s, exported figures clean)*
- [ ] Open `…260807.docx` in Word, refresh Mendeley, confirm no citation renders as plain text
