# G3 revision implementation plan — manuscript G3-2026-406828

**Manuscript:** "Evolutionary arms race between transposable elements and human genes:
telomere-to-telomere genome comprehensive analysis identifies young L1 clusters in the
interferon-alpha domain"
**Journal:** G3: Genes|Genomes|Genetics — **conditionally accepted**, minor revisions, 30-day window
(letter dated July 29, 2026; plan revised 2026-08-03, so target submission **on or before 2026-08-28**).
**Manuscript file under review:** originally `T2T_genes_subfamilies_article_figures/T2T_genes_article_for_plos_one.docx`; now `revision_G3/Revised_manuscript/T2T_genes_article_G3_revision_260803.docx`, with the read-only submitted baseline and the extended-discussion file beside it
(14.9 MB, last modified Apr 18 2026 — the filename is stale, it is the G3 submission, not a PLOS ONE one).
**Reviewer report:** `G3_reviewer_report_260802.md`
**Status:** all decisions answered by Daniil (D1–D9 and D4a, §2). No open questions.

---

## 1. Overview

The paper is conditionally accepted. Nothing in the report challenges the core result; the requests
are (a) reframing causal language, (b) adding statistical rigour around the headline IFNA–L1 finding,
(c) addressing orthogonal validation, (d) sensitivity analyses, (e) a shorter Discussion, and (f) a set
of cosmetic/consistency fixes.

Four decisions shape the whole revision:

1. **The permutation background stays at N = 500** and the Methods text is corrected down from its
   erroneous "1,000" (D1). N = 500 is what was actually run, it is what Results, Discussion and
   Figure 1D already state, and it is independently confirmed by Table 1's empirical p = 0.004 =
   2/501. Nothing in the paper is re-run for permutation reasons, which removes a whole class of
   numeric churn. Justification for 500 is *demonstrated* rather than asserted, via a convergence
   figure, and rests on the correct argument: the permutation background exists to **correct the odds
   ratio for element-length bias**, not to serve as the primary significance test (that is the Fisher
   exact p). The 11 GB legacy permutation store is compacted losslessly to ~110 MB (**WP1b**), which is
   what makes the rest of the plan fit on a filesystem with 2.4 GB free.
2. **GO FDR tightens to 0.05 everywhere** (D2), with no "suggestive 0.05–0.1 band" anywhere — main
   text, supplementary figures and supplementary tables are all filtered at 0.05. Three items leave the
   paper: *flavone metabolism* (FDR 0.088), *copper/cadmium detoxification* (0.078), and 2 of 3 MHC-I
   terms (0.051). Every GO count in every figure, table and caption is recomputed.
3. **No new epigenomic or genome-wide orthogonal analysis** (D4). Major 3 is answered in the
   Discussion, with a new subsection making a positive methodological argument: this paper supplies the
   **proximity-based null model** that the TE-epigenomic literature needs — which signals arise merely
   because certain TE families sit closer to certain genes, at T2T resolution — and the requested
   epigenomic integration exists as the companion 2026b preprint.
4. **Figures are produced only up to SVG.** All revision plotting happens in **Jupyter notebooks under
   `revision_G3/`** so Daniil can inspect and adjust each subfigure before it goes anywhere; the
   notebooks write SVGs to `revision_G3/svg/`, and Daniil places them into the Figma frames by hand.
   No script writes to Figma. §5b is therefore a **read-only reference map** (frame IDs, raster
   exposure, baked-in parameter strings) for that manual step, not a set of instructions to a tool.

Also decided: neutral title (D3); full sensitivity analysis at 5/10/20 kb windows × 5/10 % percentiles
with a dedicated robustness/intersection comparison against the original 10 kb + 5 % results (D5);
Lu et al. comparison by gene-set overlap only, with the hg38 re-run and the clustering reimplementation
both declared out of scope in the response letter (D6); Discussion cut to ~2,200 words but **keeping**
subsection *3.5 Connection of TE enrichment with cancer*, with the excised material rebuilt as a
standalone, internally coherent **"Extensive discussion"** docx (D7); GitHub Pages track hub (D8);
Zenodo DOI minted by Daniil at the end (D9).

Every manuscript edit is made as **Word tracked changes** using the existing
`~/.claude/skills/word_rewrite.md` + `word_rewrite_trackchanges.py` helpers, and must preserve the
manuscript's **128 Mendeley Cite citation content controls** — see §3.5, which documents exactly how
they are stored and why the naive approach to moving text into a second document breaks them.

Work is organised as 16 work packages (WP1–WP16 incl. WP1b) keyed to the reviewer comments, in a new
`revision_G3/` subdirectory following the modularisation rule in the project `CLAUDE.md`.

**The existing notebooks are frozen.** `TEs_mapped_on_TSS_analysis.ipynb`,
`Gene_ontology_analysis.ipynb` and `download_and_process_files_UCSC_genes.ipynb` are **not edited at
all** — not a cell, not a label, not a threshold. Every change this revision needs is implemented in
**completely new notebooks** under `revision_G3/`, and those notebooks cover **only** what this plan
requires: nothing is re-derived or re-plotted for its own sake. The consequence is that the helper
functions currently defined inside the old notebooks' cells must be **copied** into a new shared module
`revision_G3/revision_lib.py` and modified there (§7.2, §7.3, WP16, caveat C19). §7 therefore lists the
frozen files and what replaces each of them, rather than cell-level edits.

---

## 2. Decisions — resolved

| # | Decision | Daniil's answer | Consequence for the plan |
|---|---|---|---|
| D1 | Permutation count (Minor 1) | **(b)** keep N = 500 and correct the Methods, **plus** compact `permutation_results/` | WP1 becomes a text-correction + justification package; the streaming script is retained but is now used only for the *new* window sizes (WP5). New **WP1b** = lossless compaction 11 GB → ~110 MB. Empirical p floor stays **0.004**, not 0.002. No enrichment CSV, Table or Figure is regenerated for permutation reasons. |
| D2 | GO FDR threshold (Minor 2) | **0.05 everywhere**; rewrite the text and redo the plots; no suggestive band | Option (c) from the original plan is dropped entirely. Supplementary figures *and* tables are filtered at 0.05 too. flavone / copper / cadmium / 2 MHC-I terms are **removed**, not demoted. |
| D3 | Title (Major 1) | **(b)** neutral title | *"Telomere-to-telomere co-mapping of transposable elements and human genes identifies a cluster of young L1 elements in the interferon-alpha domain"*. Retitle the companion subfamilies manuscript for consistency. |
| D4 | Orthogonal data (Major 3) | **No additional epigenomic-file analysis.** Answer in the Discussion: correlative nature of the results; epigenomics is our declared future focus; cite the 2026 preprint; and make the argument that this paper is the **methodological background correction** for the epigenomic TE literature — which signals arise simply because certain TEs are closer to certain genes, at T2T scale | WP3 becomes **text-only**. `03b_encode_ifna_chromatin.py` is cut. Figure 8 loses its planned chromatin panel D. |
| D4a | GTEx tissue-specificity τ (the non-epigenomic half of the original WP3) | **Drop it.** | `03a_gtex_tissue_specificity.py` is cut. Major 3 is answered entirely in the Discussion, with **no** orthogonal dataset added — which is the coherent position: a single τ analysis would have invited "why expression but not chromatin or eQTL?" without settling the comment. No GTEx download, no supplementary figure, no Results paragraph. |
| D5 | Sensitivity scope (Major 5) | **(a)** full, with the streaming permutation script: **5/10 % percentiles + 5/10/20 kb windows**, plus a new notebook section comparing the original results (hand-written scripts, 10 kb, 5 %) against the alternatives, with intersection, robustness analysis and statistical tests | WP5 grows a dedicated comparison/robustness deliverable (`nb05_sensitivity_robustness.ipynb`). New permutation runs are needed only for 5 kb and 20 kb, at **N = 500** for comparability with the retained 10 kb background. |
| D6 | Lu et al. comparison (Major 4) | **(a) only**; state in the response letter that a full reimplementation of their clustering was out of scope | `04b_hg38_rerun.sh` and `04c_hg38_vs_t2t_compare.py` are **cut** — no hg38 re-run, no hg38 download, ~6–10 h of compute and 4 person-days saved. Also removes the planned hg38-vs-T2T supplementary table and figure. See caveat C6 for how to answer the "assembly or methodology?" question without that experiment. |
| D7 | Discussion cuts (Major 6) | Cut as proposed **but keep** *3.5 Connection of TE enrichment with cancer*; compose a separate docx **"Extensive discussion"** from the cut parts, written as integral well-connected text, not a dump. **Copy/delete in-text references as whole objects** to preserve Mendeley rendering. Do all changes as **tracked changes** | Target rises from ~2,000 to **~2,200 w** with six subsections. New deliverable: `Extensive_discussion_260803.docx`. Mendeley mechanics documented in §3.5 and caveat C10; tracked-changes mechanics in §3.6. |
| D8 | Track hub hosting (Editor) | **(a)** GitHub Pages | Unchanged from the recommendation: `gh-pages` branch, bigBed split per TE class. |
| D9 | Archival DOI | **Daniil will mint the Zenodo snapshot himself at the end** — remind him | Zenodo items in the checklist are marked `DANIIL:`. Three explicit reminders are placed: Phase 6, Phase 8, and §13 "Reminders for Daniil". |

### The frozen-notebook rule (applies to every work package)

Because the existing notebooks may not be touched, three structural rules hold throughout:

1. **New notebooks only, minimal scope.** Five notebooks in `revision_G3/` (WP16) produce every figure
   this revision needs — and nothing else. A figure that does not change is not regenerated.
2. **Helpers are copied, not edited in place.** `run_goatools_enrichment`,
   `run_goatools_ordered_enrichment`, `save_go_network_svg` (old `Gene_ontology_analysis.ipynb` cell 6),
   `visualize_go_class_network` (cell 36) and `save_go_network_svg_families_by_classes` (cell 175) are
   copied into `revision_G3/revision_lib.py`, where the FDR threshold becomes 0.05 and the new
   `min_shared_genes` / `max_term_genes` / collision-check behaviour is added. The originals keep working
   unchanged at FDR 0.1 (caveat C19 covers the drift risk).
3. **Nothing regenerates data the frozen notebooks own.** `enrichment_*.csv` and the class/family
   enrichment tables stand as published (D1), so the new notebooks *read* them rather than recomputing
   them.

---

## 3. Background / reference data

### 3.0 Canonical input files (use these, and only these, for TE ↔ TSS intersections)

Daniil copied two files into the working directory on 2026-08-03. **Every test involving TE–TSS
intersection must use them** — this supersedes the broken `../../T2T_article/T2T_repeat_masker_processed.csv`
path and the missing `T2T_genes_sorted.bed` documented in `CLAUDE.md`.

| File | Size | Rows | Columns | Verified |
|---|---|---|---|---|
| `T2T_repeat_masker_processed_sorted.bed` | 162,325,096 B | **3,709,429** | `chrom, start, end, score(divergence), subfamily, family, class` | **Byte-size-identical to `../epigenomic_files/repeats_all.bed`** and identical first records — it is the same data the legacy permutations used, so `cut -f4,5,6,7` remains the correct projection. |
| `T2T_genes.bed` | 1,203,205 B | **38,704** | `chrom, start, end, gene` | Its interval set is **exactly identical** (all 38,704 intervals, `diff`-clean) to `../epigenomic_files/BE2C.H3K9me3.chm13v2.0.mapped_on_TSS.bedGraph`, which the legacy permutation runs used as the `-b` file. The bedGraph merely carries an extra signal column. |

**Why this matters:** because both files are provably interchangeable with what the legacy pipeline
consumed, the new 5 kb / 20 kb runs and the retained 10 kb N = 500 background are directly comparable —
there is no hidden geometry change to explain in the Methods.

### 3.1 Facts verified in this repository (use these numbers, not the manuscript's)

| Item | Verified value | Source of truth | Consequence |
|---|---|---|---|
| Permutations actually run | **500** (`repeats_intersected_with_TSS_random_1..500.bed`, 501 files incl. `consolidated_random_data.csv`) | `../epigenomic_files/permutation_results/` | Methods (line 277) saying 1,000 is wrong and is corrected **down to 500**. Results / Fig. 1D already say 500. Table 1's empirical p = 0.004 = 2/(500+1) confirms it independently. |
| Empirical p-value floor | **2/501 = 0.003992 ≈ 0.004** | arithmetic; matches Table 1 | Stays as published. The IFNA test (WP2) draws its own 10,000 random windows, so its floor is 2/10,001 = 0.0002 and is unaffected by N. |
| Permutation storage | **11 GB** (500 BEDs ≈ 5.1 GB + `consolidated_random_data.csv` **6.37 GB**); free disk **2.4 GB** | `du -sh`, `ls -l`, `df -h` | WP1b compacts this losslessly; new window sizes use streaming counts (WP1/WP5). |
| Per-seed BED shape | 10,729,853 B, **545,099 rows**, 4 columns (`score, subfamily, family, class`), **70,040 unique tuples** | `repeats_intersected_with_TSS_random_1.bed` | Row order is meaningless (an unordered multiset of intersections), so run-length encoding by tuple is **lossless** for every downstream statistic. Basis of WP1b. |
| Generator script setting | `NUM_PERMUTATIONS = 1000` | `download_and_process_files_UCSC_genes.ipynb` cell 34 | **This is the source of the inconsistency** — fix the generator down to 500, not the consumer up. |
| Consumer setting | `NUM_PERMUTATIONS = 500` | `TEs_mapped_on_TSS_analysis.ipynb` cell 40; cell 36 of the download notebook | Correct as-is. No change beyond a clarifying comment and the new counts reader. |
| Cores / bedtools | 8 cores; `bedtools` 2.31.1 at `/usr/local/bin/bedtools` | `nproc`, `which` | ~8-way `xargs -P` parallelism. |
| GitHub URL in Data Availability | `…/T2T_genes_evolution` (**broken**) | manuscript line 371 | Actual remote: `https://github.com/Nikit357/T2T_transposons_genes`. This is the editor's complaint. |
| `flavone metabolic process` | **Real GO term, GO:0051552**, 5 overlapping genes (LINE class), FDR = 0.088 | `go-basic.obo`; `GO_top_5_perc_genes_by_class_number_with_all.csv` | Minor 5 is answered by citing the accession — *not* a typo for flavonoid. It is then removed by the 0.05 threshold. |
| Discussion length | 3,970 words (vs. Results 4,889; Methods 954; 123 references) | word count of the converted docx | Target ≈ 2,200 w with the cancer subsection retained. |

### 3.2 Effect of tightening GO FDR from 0.1 → 0.05 (measured, not estimated)

| Analysis level | Terms at FDR < 0.1 | Terms at FDR < 0.05 | Groups that lose all terms |
|---|---|---|---|
| Classes by count | 504 | 425 (−15.7 %) | none |
| Classes by divergence | 516 | 414 (−19.8 %) | none |
| Families by count | 196 | 160 (−18.4 %) | Dong-R4 (its single term, which the manuscript already calls "likely random") |

Headline claims, with measured FDR values. Under D2 there is **no suggestive band** — anything above
0.05 simply leaves the paper:

| Claim (abstract / title-level) | Term | FDR | Survives 0.05? |
|---|---|---|---|
| Young L1 / IFNA immune link | type I interferon receptor binding | 5.0 × 10⁻⁴ | **yes** |
| " | B / NK / T cell activation | 8.2 × 10⁻⁴ / 1.3 × 10⁻³ / 9.1 × 10⁻³ | **yes** |
| SVA → transcription termination | termination of RNA polymerase II transcription | 0.020 | **yes** |
| Alu → RNA processing/splicing | (53 → 31 terms retained) | — | **yes** |
| LINE → olfactory perception | olfactory receptor activity | 6.9 × 10⁻¹⁰ | **yes** |
| LINE → lipid metabolism | negative regulation of fatty acid metabolic process | 9.9 × 10⁻⁴ | **yes** |
| LINE → flavone metabolism | flavone metabolic process (GO:0051552) | 0.088 | **no — remove** from Abstract and Results |
| MIR → zinc | cellular response to zinc ion / Zn homeostasis | 0.006 / 0.014 | **yes** |
| MIR → copper, cadmium | detoxification of copper ion / response to cadmium | 0.078 / 0.078 | **no — remove**; MIR claim becomes zinc-only |
| hAT-Charlie → MHC class I | MHC class I protein complex | 0.025 | **yes** (1 of 3 terms) |
| " | β2-microglobulin binding; antigen processing MHC-I | 0.051 / 0.051 | **no — remove**; "three GO terms" → "one" |
| LTR → potassium channels | (LINE highest-divergence voltage-gated K⁺ complex 9.7 × 10⁻⁴) | — | **yes** |

### 3.3 The IFNA domain, pre-computed (chr9:21,150,692–21,370,055)

| Quantity | Value |
|---|---|
| All TEs in the 220 kb window | 175 |
| **L1 elements** | **77** (44 % of all TEs in the window) |
| Other content | Alu 33, L2 15, MIR 13, hAT-Charlie 11, ERV1 10, ERVL-MaLR 6, others 10 |
| L1 density in window | 350 L1 / Mb vs. **182 L1 / Mb genome-wide** (565,459 L1 / 3.1 Gb) → ~1.9× |
| Mean L1 divergence in window | **135.7** vs. genome-wide mean **188.2** (median 197) |
| Youngest elements | L1P3 (div 0 — needs QC, see caveat C7), L1PA3 (19), L1PA2 (24), L1P1 (25), L1PA4 (28), L1PA5 (37), L1PA3 (39), L1PA4 (39) |
| Subfamilies present | primate-specific young: L1PA2/3/4/5/6/7/8/10/14/15, L1P1/3/4e/5, L1PB, L1PREC2; older mammalian: L1MA1/2/3/5A/6/9, L1MB2, L1MC/3/4a, L1MD3, L1M2c/3/4, L1ME3A/4b/4c, L1MEf, L1MCb |

**This directly answers Major 2's "not driven by a few outliers":** 77 elements, spanning ≥ 30
subfamilies, with a clear excess of the young primate-specific L1PA/L1P clades.

**These descriptive numbers go into the Results text**, not only into the response letter — a new
paragraph in the divergence section reporting the element count, the class breakdown, the subfamily
inventory, the window-vs-genome density ratio and the two divergence means. The formal tests that
accompany them are WP2. Checklist item in Phase 5.

### 3.4 Reviewer comment → work package map

| Reviewer item | WP | Type | Effort |
|---|---|---|---|
| Minor 1 — 500 vs 1,000 permutations | WP1 | text + justification figure | 0.5 d |
| (D1 companion) compaction of the 11 GB permutation store | **WP1b** | new script | 0.5 d + ~1–2 h compute |
| Major 2 — IFNA L1 permutation test + counts/subfamilies | WP2 | new analysis + new figure | 2 d |
| Major 3 — orthogonal data (ENCODE/GTEx) | WP3 | **text only** (D4) | 0.5 d |
| Major 4 — direct comparison with Lu et al. (ref 14) | WP4 | gene-set overlap only (D6) | 1 d |
| Major 5 — window (5/10/20 kb) and percentile (5/10 %) sensitivity + robustness comparison | WP5 | new analysis + supp. tables + notebook | 3.5 d |
| Minor 2 — FDR 0.1 → 0.05 | WP6 | re-run GO + text | 1 d |
| Major 1 — reframe "arms race"/causality | WP7 | text only | 0.5 d |
| Major 6 — shorten/refocus Discussion + Extensive discussion docx | WP8 | text + new deliverable | 2 d |
| Minor 3 — p-value provenance in all axis labels/legends + labelling reference doc | WP9 | figure/caption edits + new doc | 1.3 d |
| Minor 4 — Table 1 too wide | WP10 | table restructure | 0.5 d |
| Minor 5 — "flavone metabolism" | WP11 | text only | 0.1 d |
| Minor 6 — simplify Figures 4A/5A/6A, no overlapping text | WP12 | re-plot + supplementary | 2 d |
| Editor — data/code links | WP13a | repo + text | 0.5 d |
| Editor — browser/track instance | WP13b | new deliverable | 2 d |
| Response letter, tracked changes, submission package | WP14 | manuscript production | 2 d |
| Final formatting to G3 house style (`G3_article_guidelines.md`, gaps G1–G18) | WP15 | manuscript production | 2 d (G1 is Daniil's) |
| Notebook-based figure production surface + `revision_lib.py` (cross-cutting; existing notebooks frozen) | WP16 | infrastructure | 1.5 d |

**≈ 22 person-days plus ~6–10 h of compute** against a window that started July 29. The long poles are
WP5, WP8 and WP12. Compute is no longer a critical path (D1 and D6 removed the two big re-runs).

### 3.5 How the manuscript stores its citations (verified — read before any text edit)

Daniil's D7 instruction to "copy/delete in-text references as a whole to preserve Mendeley citation
rendering" is exactly right, and the mechanism is more specific than field codes:

| Fact | Value |
|---|---|
| Citation mechanism | **Mendeley Cite** (the Word *web add-in*), not the legacy COM plugin |
| In-text citations | **128** `<w:sdt>` content controls in `word/document.xml`, each tagged `MENDELEY_CITATION…` |
| Bibliography | **1** `<w:sdt>` tagged `MENDELEY_BIBLIOGRAPHY` |
| Citation payload | a single `MENDELEY_CITATIONS` property inside `word/webextensions/webextension1.xml` (**2.5 MB**) |
| Style | `MENDELEY_CITATIONS_STYLE` = *NLM/Vancouver: Citing Medicine 2nd edition (citation-sequence)*, `format: numeric` |
| Locale | `MENDELEY_CITATIONS_LOCALE_CODE` = `en-GB` (consistent with the British spellings G3 will want changed) |
| Bibliography state | `MENDELEY_BIBLIOGRAPHY_IS_DIRTY = true` — it already needs a refresh in Word |
| Legacy field codes | **none** — zero `ADDIN CSL_CITATION`, `instrText`, `fldSimple`, `fldChar` anywhere in the archive |

Three consequences that change how the work is done:

1. **Moving or deleting a citation means moving or deleting the entire `<w:sdt>` subtree.** Touching
   only the visible run inside it leaves an orphaned content control or, worse, plain text that Mendeley
   will not re-render on refresh. All tracked edits operate on `<w:sdt>` elements as units.
2. **The "Extensive discussion" docx must be built by copying the manuscript file and deleting
   everything else** — never by creating a fresh document and pasting text. The citation payload lives
   in a webextension part that a fresh document does not have; paste-based construction produces dead
   plain-text citation numbers. Procedure in WP8.
3. **G1 (numbered → author–year) is a Mendeley style switch**, not a text edit: Daniil changes
   `MENDELEY_CITATIONS_STYLE` from Vancouver to the G3/CSE author–year style in the Mendeley Cite pane
   and refreshes. That is why D7/WP15 assign G1 to him. What *we* must still do by hand is the
   **in-sentence grammar** that numbered citations allowed and author–year does not
   (*"the previous landmark study by (14)"* → *"the previous landmark study by Lu et al. (2020)"*),
   because that text sits outside the content controls.

### 3.6 Tracked-changes mechanics

Daniil asked for all changes as tracked changes, "like in the skill /nar-review". There is no
`/nar-review` skill in this environment; the mechanism he means is
**`~/.claude/skills/word_rewrite.md`** with its helper library
**`~/.claude/skills/word_rewrite_trackchanges.py`**, whose guiding principles were written from the NAR
rewrite. It emits genuine `<w:ins>`/`<w:del>` revision markup, so *Accept All* yields the clean
revision and *Reject All* restores the submitted baseline exactly — which is also precisely what the
G3 letter asks for ("a clean version" plus "a highlighted or tracked version").

Two operational notes:

- The manuscript currently contains **zero** `<w:ins>`/`<w:del>` elements, i.e. a clean baseline. Good:
  the tracked diff will be unambiguous.
- The helper library needs `python-docx` + `lxml`. They are in **`~/venvs/collagen_3_11`**, not in
  `~/venvs/Retroelements_3_11`. Either run the docx work in the collagen venv or `pip install
  python-docx` into the Retroelements venv — decide once and record it in `CLAUDE.md` (caveat C15).

---

## 4. Work packages

### WP1 — Justify and correct the permutation count at N = 500 (Minor 1)

**Reviewer:** "The number of random permutations is stated as 500 in the Results and Discussion but as
1000 in the Methods. Please unify this and justify the chosen number."

**Response strategy (D1).** N = 500 is what was run; the Methods number was a drafting error. Correct
the Methods down to 500, and justify 500 on the grounds that are actually true:

1. **The permutation background is a bias correction, not the significance test.** Raw Fisher odds
   ratios scale with element length (random OR vs. mean element length, **Pearson R = 0.985** across
   the 44 families, n = 44, p = 8.4 × 10⁻³⁴ — re-measured 2026-08-04; this line and `CLAUDE.md`
   previously said 0.661, which was wrong), so SINEs
   (~300 bp) would be systematically under-called and LINEs/LTRs (~6 kb) over-called. The permutation
   distribution exists to produce the observed/random OR ratio; significance calls come from the
   FDR-adjusted Fisher exact test. What N must deliver is therefore a **stable mean and SD of the random
   OR**, not a small p-value floor.
2. **That stability is reached well before 500.** Report the observed random-OR SDs (0.005–0.163 in
   Table 1) and add a **convergence figure**: running mean and running SD of the random OR per class and
   per family against iteration index, showing the curves flat long before 500. This demonstrates the
   justification instead of asserting it.
3. **State the floor honestly.** The empirical p floor is 2/(N+1) = **0.004** at N = 500. Say so, and
   say that the empirical p is reported as a supporting statistic alongside the Fisher exact p rather
   than as the primary test — which is why a floor of 0.004 is not limiting. Note explicitly that for
   families/classes whose empirical p sits at the floor, BH correction across 44 families still reaches
   significance because many hit the floor simultaneously; a single isolated floor hit would not, and we
   do not rely on one.
4. **The IFNA test is unaffected.** WP2 draws its own 10,000 matched random windows (floor 0.0002), so
   nothing about the headline finding depends on N = 500.

**Scripts:**

- `revision_G3/nb01_permutation_convergence.ipynb` — reads the compacted store (WP1b) and emits the
  convergence panels to `revision_G3/svg/S_permutation_convergence_*.svg`. New **Supplementary
  Figure** for the Methods justification.
- `revision_G3/01_permutations_stream.sh` (+ `01a_consolidate_counts.py`) — the streaming permutation
  runner. **Retained, but its purpose has changed:** it is no longer used to re-run the 10 kb
  background; it generates the 5 kb and 20 kb backgrounds for WP5, at N = 500, in the compact format.

```bash
# revision_G3/01_permutations_stream.sh  (N=500, one window size per invocation)
#   INPUT_REPEATS = T2T_repeat_masker_processed_sorted.bed   (§3.0; == epigenomic_files/repeats_all.bed)
#   TSS_FILE      = windows_5kb.bed | T2T_genes.bed | windows_20kb.bed
run_one () {
  SEED=$1
  bedtools shuffle -i "${INPUT_REPEATS}" -g "${GENOME_FILE}" -seed "${SEED}" \
    | bedtools intersect -a - -b "${TSS_FILE}" -wa \
    | cut -f4,5,6,7 \
    | sort | uniq -c \
    | awk -v s="${SEED}" '{print s"\t"$2"\t"$3"\t"$4"\t"$5"\t"$1}' \
    > "${OUTPUT_DIR}/counts_seed_${SEED}.tsv"
}
export -f run_one
seq 1 500 | xargs -P "${N_JOBS}" -I{} bash -c 'run_one {}'
```

- Output columns are `seed, score, subfamily_name, family_name, class_name, n` — **the same schema
  WP1b produces**, so compaction and the new runs converge on one format and one reader.
- Identical seeds 1–500 and identical shuffle/intersect semantics as the legacy run, so
  `01a_consolidate_counts.py --check-legacy` can prove the pipeline reproduces
  `consolidated_random_data.csv` before anything is deleted.
- Per-window output ≈ 110 MB uncompressed counts, ≈ 20 MB compressed.

**Text edits:** Methods line 277 `1000` → `500` plus the justification paragraph; verify the Results
(line 37, 155) and Figure 1D caption (line 33) already say 500; keep "empirical p-value = 0.004"
everywhere it appears — **it does not change.**

**Downstream regeneration: none.** Table 1/2 values, Figure 1D bars, `enrichment_families_with_random.csv`
and `enrichment_subfamilies_with_random*.csv` all stay as published. This is the single largest
simplification D1 buys.

---

### WP1b — Lossless compaction of `permutation_results/` (D1 companion)

**Goal (Daniil):** compress the 11 GB as much as possible while preserving the permuted results.

**Why it is losslessly compressible by ~100×.** Each per-seed file is an *unordered multiset* of
4-tuples `(score, subfamily, family, class)` — one row per intersecting TE, in arbitrary order. Row
order carries no information, so run-length encoding by tuple loses nothing that any downstream
statistic (count, mean, distribution, quantile) can detect. Measured on
`repeats_intersected_with_TSS_random_1.bed`:

| Representation | Size | Ratio |
|---|---|---|
| Raw BED, 545,099 rows | 10,729,853 B | 1× |
| Same content, `gzip -9` | ~1.7 MB | ~6× |
| **Counts TSV, 70,040 unique tuples** | 1,693,541 B | **6.3×** |
| **Counts TSV + `gzip -9`** | 275,224 B | **39×** |
| **Counts TSV + `zstd -19`** | **217,587 B** | **49×** |

Across the store: 500 × 217,587 ≈ **109 MB**, and `consolidated_random_data.csv` (**6.37 GB**) becomes
redundant — it is only the concatenation of the per-seed files with a `permutation_id` column, which the
compact store carries as `seed`. **Total: 11 GB → ~110 MB, ~100×, lossless.**

**New script:** `revision_G3/01b_compact_permutation_results.py`

Procedure, designed to be safe on a filesystem with 2.4 GB free:

1. Process **one seed at a time**: read `repeats_intersected_with_TSS_random_${SEED}.bed`, group by
   `(score, subfamily, family, class)`, write
   `revision_G3/output/permutation_counts_10kb/counts_seed_${SEED}.tsv.zst`.
2. **Verify each file before deleting anything**: total row count must equal the source line count, the
   per-class totals must match, and the reconstructed multiset must be identical (compare
   `sort | uniq -c` of the source against the counts file). Only then delete the source BED. Peak extra
   disk use is one file (~11 MB), not the whole store.
3. After all 500: run the aggregate check against `consolidated_random_data.csv` — per
   `(permutation_id, class)` and per `(permutation_id, family)` totals, and the divergence distribution
   per class (mean, SD, deciles) — then delete the 6.37 GB CSV.
4. Write `revision_G3/output/permutation_counts_10kb/MANIFEST.json`: source file list, row counts,
   checksums, compression settings, script version, date. This is the provenance record that makes the
   deletion defensible.
5. Provide `revision_G3/01c_expand_counts.py --seed N` to reconstruct an exact legacy BED on demand, so
   the old format is never truly lost.

**Reader.** The frozen `TEs_mapped_on_TSS_analysis.ipynb` reads the 6.37 GB CSV and will simply stop
working once it is deleted — that is expected and acceptable, because its outputs are already final
under D1. The replacement reader lives in `revision_G3/revision_lib.py` as
`load_permutation_counts(window="10kb")`, returning a DataFrame with `n` as a weight column. It is what
`nb01`/`nb03`/`nb05` use, and it is fast and low-memory where the old path was neither. Note the
consequence in `REPRODUCE.md`: anyone re-running the old notebook must first call
`01c_expand_counts.py` to rebuild the legacy files.

**Order of operations:** WP1b runs **first**, before anything else needs disk — it is what frees the
space for WP5's two new window backgrounds.

**Fallback if a check fails:** stop, keep the source, and fall back to plain `zstd -19` on the BEDs
in place (~6× → 11 GB becomes ~1.8 GB). Still enough to proceed. Never delete a source file whose
verification did not pass.

---

### WP2 — Statistical test for the IFNA domain (Major 2)

**Reviewer:** "…do not compare this to genome-wide background or to matched random regions. Please
perform a statistical test (e.g., permutation test)… Also, please report the number of L1 elements in
this region and their specific subfamilies."

**New script:** `revision_G3/02_ifna_domain_test.py` — **inputs are the §3.0 files**
(`T2T_repeat_masker_processed_sorted.bed`, `T2T_genes.bed`), not the broken CSV path.

Four tests:

1. **Divergence, unmatched genome-wide background.** Mean divergence of the 77 window L1s (135.7) vs.
   10,000 random 220 kb autosomal windows (`bedtools shuffle -incl` a mappable-region file, excluding
   centromeric/acrocentric arrays, requiring ≥ 1 L1). Empirical two-sided p (floor 0.0002).
2. **Divergence, L1-count-matched background.** Same, restricted to windows containing ≥ 40 L1
   elements — controls for the trivial possibility that low mean divergence follows from high L1
   density. This is the test the reviewer's phrase "matched random regions" is asking for.
3. **Divergence, gene-density-matched background.** Restricted to windows containing ≥ 10 genes from
   `T2T_genes.bed` (the window has 12).
4. **Subfamily composition test.** 2×2 Fisher: young primate-specific L1 (L1HS/L1P*/L1PA*/L1PB/L1PREC*)
   vs. older L1M* in the window, against the genome-wide L1 composition. Plus a leave-one-out and a
   trimmed mean (drop the 5 youngest elements) to show the signal is not carried by them.

**Persistence and review (Daniil's instruction).** The script writes every result into
`revision_G3/output/` as CSV/TSV — `ifna_window_elements.csv` (all 175 TEs), `ifna_test_results.csv`
(the four tests with observed value, null mean/SD, empirical p), `ifna_null_distributions.csv` (the
10,000 null values per test) — and **`revision_G3/nb02_ifna_domain.ipynb` reads those files and renders
them** so Daniil can see the numbers and the figures before anything is placed. The notebook writes
panel SVGs to `revision_G3/svg/`; it does not touch Figma.

**Also fix the framing — two distinct quantities, both stated in the text.** The manuscript currently
reports "average divergence of intersecting LINE elements at the level of 95–161.7", which is the
**per-gene average over the 8 IFNA TSS windows**, not the divergence of the elements themselves. Both
go into the paper, separately labelled:

- *per-gene mean LINE divergence across the 8 IFNA TSS windows: 95 – 161.7*
- *mean divergence of the 77 L1 elements in the 220 kb domain: 135.7 (genome-wide L1 mean 188.2)*

**Outputs:**
- New **main Figure 8 "The interferon-alpha domain"** (promotes the current Supplementary Figure 6):
  (A) UCSC browser view of chr9:21,150,692–21,370,055 (existing vector panel); (B) null distributions
  for tests 1–3 with the observed value marked; (C) L1 subfamily composition and per-element divergence
  strip/lollipop plot. *There is no chromatin panel D — D4 removed it.* Current Figures 8 and 9
  (schematics) become Figures 9 and 10.
- New Supplementary File: per-element table of the 175 TEs in the window.
- Results: the descriptive paragraph from §3.3 **plus** the four test results, and a rewritten
  Discussion paragraph (line 239).

**Risk:** medium-low. The signal is large (135.7 vs 188.2) and n = 77, so tests 1 and 3 will almost
certainly be significant. Test 2 (count-matched) is the one that could come back weaker — if it does,
report it as such and soften the claim to "the domain combines above-average L1 density with a
significant excess of young primate-specific L1 subfamilies", which test 4 still supports.

---

### WP3 — Answer Major 3 in the Discussion, without new analysis (D4)

**Reviewer:** "The study relies solely on a single genome assembly (T2T-CHM13) and does not integrate
orthogonal data (e.g., epigenetic marks, expression data, eQTL) to support functional relevance. It
would be helpful to add such analyses using public resources (e.g., ENCODE, GTEx, TCGA)."

**Decision (D4 + D4a): no new epigenomic-file analysis and no GTEx.** The answer is a new, short
Discussion subsection plus a paragraph in the response letter. This is a defensible answer because the
reviewer wrote "it would be helpful", not "this is required", and because the requested integration
already exists as our companion study.

**New Discussion subsection — *Proximity as a null model for TE–epigenomic association studies*
(~250 w).** The four points Daniil specified, in this order:

1. **The design is correlative and static.** This is a co-localisation map of TE positions relative to
   TSS neighbourhoods in one assembly. It measures where elements are, not what they do; functional
   relevance requires perturbation.
2. **The positive methodological argument.** A large literature reports that some TE family carries
   some epigenomic mark near genes of some function, and reads that as evidence of regulatory
   recruitment or domestication. Any such claim requires a **positional baseline**: a share of the
   signal follows from nothing more than certain TE families being physically closer to certain gene
   classes than chance allows. This paper supplies exactly that baseline — genome-wide, at
   telomere-to-telomere resolution, with a length-bias-corrected permutation background across 6
   classes and 44 families — and it is the reference against which mark-based enrichments should be
   normalised. Framed this way, the absence of epigenomic data here is the point rather than a gap:
   the proximity layer has to be measured on its own before it can be conditioned on.
3. **Epigenomics is our declared next step**, and it exists: cite the 2026b companion preprint
   (DOI 10.64898/2026.03.19.712972 — 7 marks × 12 cell lines × 1,122 subfamilies on T2T), which is
   precisely the orthogonal integration the reviewer asks for, and state that the present proximity map
   is its background correction.
4. **Named limitations and future tests** (this doubles as Major 6's item 4): no eQTL/TCGA integration;
   single assembly; GO annotation bias; multiple-TSS bias (already flagged at line 217); TE annotation
   uncertainty in newly resolved regions. Concrete follow-ups: CRISPRi of individual IFNA-domain L1s
   under IFN stimulation; SVA_B deletion at the SSU72L cluster; conditioning published mark-based
   enrichments on this paper's proximity baseline.

**Response-letter framing:** state plainly that we did not add the analysis, why the proximity baseline
is a prerequisite rather than a substitute, and where the epigenomic integration lives. Disclose the
companion preprint in the cover letter so the editor does not read it as undisclosed overlap (C4).

**No scripts. No new figures. No new supplementary files.** `03a_gtex_tissue_specificity.py` and
`03b_encode_ifna_chromatin.py` are both removed from the deliverable list.

---

### WP4 — Comparison with Lu et al. 2020, gene-set level only (Major 4, D6)

**Reviewer:** "The author notes methodological differences… but does not directly compare results on
the same dataset. Is this difference due to methodology or the updated genome assembly?"

Reference 14 = Lu JY, Shao W, Chang L et al. *Genomic Repeats Categorize Genes with Distinct Functions
for Orchestrated Regulation.* Cell Rep 2020;30(10):3296.

**New script:** `revision_G3/04a_lu2020_geneset_overlap.py` (option (a) only, per D6)

Pull their published gene categories from the Cell Reports supplementary tables, map to current HUGO
symbols, and cross-tabulate against our top-5 % sets per class (Fisher + Jaccard + a supervenn panel,
reusing the existing supervenn code path in `Gene_ontology_analysis.ipynb` cells 32/90). Deliverable: a
matrix showing which of their repeat-based gene groups our SINE/LINE/LTR/DNA-enriched sets correspond
to, and where they diverge.

**Answering "assembly or methodology?" without the hg38 re-run.** D6 drops the assembly-controlled
experiment, so the answer must be argued from evidence we do have, and the response letter must be
candid that it is an argument rather than a controlled test:

- Quantify the **assembly-attributable** component descriptively: how many of our TEs and TSS windows
  fall in regions absent or misassembled in hg38 (newly resolved acrocentric arms, centromeric
  regions, segmental duplications). This bounds how much of the difference the assembly *could*
  explain.
- Attribute the rest to methodology, and name the specific differences already in the manuscript
  (TSS-window proximity vs. region-binning; per-element vs. per-bin normalisation; length-bias
  correction via permutation).
- State explicitly in the response: *"a full reimplementation of their region-binning and clustering
  was out of scope for this revision"*, and note the gene-set overlap as the direct comparison we can
  make on shared quantities.

**Outputs:** one new supplementary table (category × class overlap with Fisher p and Jaccard), one new
supplementary figure (supervenn + overlap heatmap), and a rewritten Discussion paragraph replacing
line 217's speculative "the comparison of both methods… could be helpful" with the actual overlap
result and the bounded assembly argument.

**Cut from the original plan:** `04b_hg38_rerun.sh`, `04c_hg38_vs_t2t_compare.py`, the hg38
RepeatMasker/RefSeq downloads, 1,000 hg38 permutations, and the hg38-vs-T2T table and figure.

---

### WP5 — Window and percentile sensitivity, with a robustness comparison (Major 5, D5)

**Reviewer:** "…sensitivity analyses for key findings (e.g., the IFNA-L1 association, the
SVA-termination association) using alternative window sizes (5 kb and 20 kb) and alternative
percentiles (e.g., 10%)."

**Inputs:** the §3.0 files throughout.

**New scripts / notebook:** `revision_G3/05a_build_windows.sh`, `05b_window_sensitivity.py`,
`05c_percentile_sensitivity.py`, `nb05_sensitivity_robustness.ipynb`

- **Windows.** Build TSS neighbourhoods at ±2.5 kb (5 kb) and ±10 kb (20 kb) from the same TSS
  definition that produced `T2T_genes.bed`, re-map TEs, and run **500** streaming permutations per
  window (WP1 machinery, ≈ 20 MB compressed each). 10 kb reuses the retained legacy background — this is
  legitimate precisely because `T2T_genes.bed` and the legacy TSS file have identical interval sets
  (§3.0).
- **Full enrichment recomputation** at 5/10/20 kb for 6 classes and 44 families → one supplementary
  table with observed OR, random OR, obs/random, and both p-values per window.
- **Percentiles.** Re-run GO for the top/bottom **10 %** gene sets (2,872 genes) alongside the original
  **5 %** (1,436), for all classes, divergence groups and families, at FDR 0.05.
- **Key-finding spot checks**, named explicitly in the response letter: (i) IFNA — re-run all four WP2
  tests at 5 and 20 kb *and* confirm the low-divergence-LINE GO gene set still contains the IFNA core
  genes; (ii) SVA — confirm "termination of RNA polymerase II transcription" persists at 5 and 20 kb.

**The robustness comparison Daniil asked for** — a dedicated section of
`nb05_sensitivity_robustness.ipynb` comparing the **original results (hand-written scripts, 10 kb,
top/bottom 5 %)** against every alternative:

| Comparison | Method |
|---|---|
| Enrichment agreement across windows | Spearman and Pearson correlation of obs/random OR across the 6 classes and 44 families, pairwise between 5/10/20 kb; Bland–Altman plot of the differences; count and name every significance flip |
| Gene-set stability | Jaccard and overlap coefficient between the top/bottom gene sets at each window and percentile; hypergeometric p for each overlap against chance; supervenn of the three window sets per class |
| GO-term stability | Jaccard between the FDR < 0.05 term sets; fraction of original terms preserved; every term gained or lost, listed |
| Rank stability | Kendall τ between gene rankings by TE count across windows; bootstrap CI on τ |
| Headline-claim table | one row per abstract-level claim × 6 conditions (3 windows × 2 percentiles), marked survives / weakens / lost |
| Statistical test of concordance | permutation test on the observed correlation against label-shuffled null, so "concordant" is a measured statement rather than an eyeball judgement |

Expect rank correlation ≳ 0.95 between windows; if not, that is itself a reportable result about the
spatial scale of the associations.

**Outputs:** 2 new supplementary tables, 1 new supplementary figure (OR concordance across windows +
GO-term overlap across percentiles), the robustness section in the notebook, 1 new Results paragraph,
1 Methods paragraph. Panel SVGs to `revision_G3/svg/`.

---

### WP6 — GO FDR 0.05 everywhere (Minor 2, D2)

**Reviewer:** "The GO analysis uses an FDR threshold of 0.1. Please justify this relatively lenient
threshold or, ideally, tighten it to 0.05."

**Decision (D2): 0.05 only.** No suggestive band, in main text, supplementary figures **or**
supplementary tables. Every figure and table that currently reports GO results at 0.1 is redone.
The threshold change is made in the **copied** helpers in `revision_G3/revision_lib.py`
(`fdr_threshold=0.05`), not in the frozen notebooks or in `GO_subfamilies.py`; `06_go_rerun_fdr005.py`
drives them. Measured consequences are in §3.2. Text edits:

- **Abstract:** drop "flavone metabolism" from the LINE claim; change the MIR claim from "zinc, copper
  and cadmium detoxification" to "zinc ion response and homeostasis".
- **Results line 93:** remove flavone from the LINE sentence.
- **Results line 119:** "14 families" → recount at 0.05 (Dong-R4 loses its only term → 13 families);
  recompute the "3 of these families had only 1 significant GO term" sentence.
- **Results line 121:** drop flavone from the L1 sentence.
- **Results line 127 / Discussion line 229:** hAT-Charlie "three GO terms" → **one** (MHC class I
  protein complex, FDR = 0.025). The other two are removed, not demoted.
- **Results line 125 / Discussion line 223:** MIR metals wording → zinc only.
- **Methods line 287:** threshold 0.1 → 0.05. No sentence about a suggestive band.
- **Every GO term count** in Figures 4B, 5B, 6B, 7A–7H, Supplementary Figure 8B and their captions is
  recomputed — these are counts of significant terms, so they all shift.
- **Supplementary tables** are re-emitted with the 0.05–0.1 rows removed. The CSVs still store the full
  FDR column internally, but nothing published shows a term above 0.05.
- **Inside the Figma figures (§5b) — Daniil's manual step:** Figure 7 (frame `861:34`) contains the
  literal axis label `"GO terms count in a group (FDR < 0.1)"`, which must become `0.05`. Grep every
  frame for `FDR < 0.1`, `0.1` and `p-value` once WP6 lands; parameters baked into figure text are the
  easiest way to end up contradicting our own Methods (C11).

**Script:** `revision_G3/06_go_rerun_fdr005.py` — re-runs all GO levels at FDR 0.05 and writes to
`revision_G3/output/`. Plotting happens in `nb06_go_networks_fdr005.ipynb` (WP12).

---

### WP7 — Reframe causality and the "arms race" narrative (Major 1, D3)

**Reviewer:** "…the data are purely correlative and static… Please revise throughout the manuscript to
reflect the descriptive, hypothesis-generating nature of the findings."

Systematic language pass (line numbers are in the converted text; map to the docx by string search):

| Location | Current | Revision |
|---|---|---|
| Title | "Evolutionary arms race between transposable elements and human genes…" | per D3: *"Telomere-to-telomere co-mapping of transposable elements and human genes identifies a cluster of young L1 elements in the interferon-alpha domain"* |
| Abstract (line 13) | "suggesting a recent evolutionary arms race influencing innate immune responses" | "consistent with recent L1 activity in this domain; whether these elements influence innate immune gene regulation remains to be tested" |
| Abstract (line 13) | "LTR elements were potentially associated with…" | keep "associated"; remove all "influence/impact" verbs applied to our own data |
| Intro (line 21) | "indicating the recent example of evolutionary arms race shaping innate immune response" | "identifying a candidate locus of recent L1 activity for functional follow-up" |
| Results (lines 103, 107) | "indicating possible molecular or evolutionary mechanism" | "a pattern whose mechanism is not addressed by this design" |
| Discussion (line 239) | "suggests an ongoing evolutionary arms race" | "is compatible with recent L1 activity; we cannot distinguish insertion bias, reduced purifying selection, and local L1-attracts-L1 feed-forward expansion with these data" |
| Throughout | "impact", "influence", "drive", "shaping" applied to our results | "co-localise with", "are enriched near", "are associated with" |
| New | — | the Limitations subsection (WP8 §3) and the proximity-null subsection (WP3) |

Keep "arms race" only where it refers to *cited literature*, not to our own inference, and say so once
explicitly. Note that the title change also implies retitling the companion subfamilies manuscript.

---

### WP8 — Shorten and refocus the Discussion; build the Extensive discussion docx (Major 6, D7)

**Reviewer:** "…substantially shorten and refocus the Discussion on: (1) what the author found; (2) how
it compares to prior work; (3) limitations; and (4) specific hypotheses for future testing."

Restructure 3,970 w → **~2,200 w** in six subsections. Per D7, the cancer subsection is **kept**:

| New subsection | Target | Sourced from | Action |
|---|---|---|---|
| 1. Principal findings | ~400 w | lines 149, 195, 197, 199, 205, 207, 237 | compress; keep the 99.11 %-of-TSS framing, the class/family enrichment summary, the IFNA result |
| 2. Comparison with prior work | ~450 w | lines 153, 217, 219, 221, 225 + WP4's overlap result | replaces line 217's speculation with the measured overlap and the bounded assembly argument |
| 3. Limitations | ~350 w | new (WP7) + lines 191, 217's caveat, 231 | new subsection |
| 4. Hypotheses for future testing | ~300 w | new; derived from WP3 | CRISPRi of IFNA-domain L1s under IFN stimulation; SVA_B deletion at the SSU72L cluster; conditioning published mark-based enrichments on this proximity baseline |
| **3.5 Connection of TE enrichment with cancer** | **~150 w** | lines 243, 245 | **KEPT per D7** — condensed to the core claim with a much reduced citation set; the biomarker/gene-signature material (refs 86–97) moves to the Extensive discussion |
| 5. Proximity as a null model for TE–epigenomic studies | ~250 w | new (WP3/D4) | the methodological-baseline argument + 2026b citation |
| 6. Mechanistic framework (Figure 9/10) | ~300 w | lines 161–169 | keep the 4-mechanism list + schematic, trim the supporting prose |

Cuts (~2,100 w moved out, not deleted):
- Window-size literature review (lines 153, 155; 346 w) → ~100 w; the rest moves to Methods (this also
  serves WP5's Methods paragraph).
- Per-class/per-family mechanism review (lines 175–191, ~1,000 w) → ~350 w; keep only what bears on our
  own enrichment values.
- Cancer biomarker/signature material (line 245, 179 w, 12 self-citations) → out of the main text,
  into the Extensive discussion; the *3.5* subsection survives as a condensed ~150 w statement.

**The "Extensive discussion" deliverable (D7).** `Extensive_discussion_260803.docx` — the excised
material rewritten as an integral, well-connected text, not a pile of orphaned paragraphs: its own
short introduction stating what it extends and why it was moved, then thematic sections
(window-size choice in the literature; per-class mechanistic review; TE enrichment and cancer in
depth), with connective sentences so it reads as a standalone essay. Submitted as a supplementary
file (`File Sn`), cited once from the Discussion.

**Build procedure — this is the part that can go wrong (§3.5).** Do **not** create a new document and
paste. Instead:

1. Copy `T2T_genes_article_G3_revision_260803.docx` → `Extensive_discussion_260803.docx`. This carries
   `word/webextensions/webextension1.xml`, and with it the `MENDELEY_CITATIONS` payload, so the
   citations in the moved text stay live.
2. In the copy, delete everything except the retained extended-discussion text, deleting each
   in-text citation's **entire `<w:sdt>` content control** with its paragraph rather than the visible
   run alone.
3. Keep the `MENDELEY_BIBLIOGRAPHY` content control so Daniil can refresh in Word and get a
   bibliography containing only the references still cited in that file.
4. In the main manuscript, delete the same passages as **tracked deletions** (§3.6), again moving whole
   `<w:sdt>` elements.
5. Daniil opens both files in Word, refreshes Mendeley in each, and confirms no citation renders as
   plain text and no `[EDITOR NOTE]` placeholder remains.

Reducing the main reference list from 123 toward ~100 is a side benefit and is appropriate for a
research article.

---

### WP9 — p-value provenance in every figure (Minor 3)

**Reviewer:** "…clearly indicate in axis labels or legends whether p-values are raw or FDR-adjusted…
In Figure 1D, for example, the third vertical bar plot labeled '-log10(p-value)' should read
'-log10(FDR-corrected p-value)'."

Figure 1D is produced by the frozen `TEs_mapped_on_TSS_analysis.ipynb` **cell 63**, and the plotted
column is already `Enrichment_p_value_adjusted` — the label was simply wrong. Because that notebook may
not be edited, the panel is **re-created in `revision_G3/nb03_relabelled_figures.ipynb`**: copy cell 63's
plotting code into the new notebook, read the unchanged `enrichment_families_with_random.csv`, and fix
the two labels there. The data is identical, so the panel is visually identical apart from the labels.

```python
# The wrong labels, as they stand in the frozen notebook (cell 63)
ax_p.set_title("Significance\n-log10(adj p)", fontsize=10)
ax_p.set_xlabel("-log10(p)", fontsize=9)

# What nb03_relabelled_figures.ipynb emits instead
ax_p.set_title("Significance\n-log10(FDR-corrected p)", fontsize=10)
ax_p.set_xlabel("-log10(FDR-corrected Fisher exact p-value)", fontsize=9)
```

Then audit every remaining figure and caption for the same problem. Everything in the "Producing cell"
column below refers to the **frozen** notebooks and is the code to copy into `nb03` (or `nb06` for the
GO panels), never to edit in place:

| Figure | Producing cell | Statistic | Action |
|---|---|---|---|
| 1D bars 1 & 2 | `TEs…ipynb` 63 | OR; obs/random fold change | add "empirical p (FDR-corrected)" to the "ns" legend text |
| 2D–2F | `TEs…ipynb` 105 | Mann-Whitney via statannotations | caption already says FDR-corrected — make the axis/legend say it too |
| 3A/3B | `TEs…ipynb` 171–175 | Kolmogorov-Smirnov | state raw vs FDR in caption; currently unstated |
| 4B, 5B, 6B, Supp 8B | `Gene_ontology_analysis.ipynb` 70, 118, 191, 211, 318 | Fisher, star annotations | stars are FDR-corrected — label the colourbar/legend accordingly |
| 4A, 5A, 6A node colours | GO nb 6, 36, 175 | goatools `p_fdr_bh` | legend must read "FDR-corrected GO enrichment p-value" |
| 7A–7G | GO nb 196–202 | Pearson r / Mann-Whitney | state raw (these are single tests) explicitly |
| Supp 1–3 | `TEs…ipynb` 181, 194 | Mann-Whitney FDR | already stated; verify |
| 8B (new) | `nb02_ifna_domain.ipynb` | empirical permutation p | label "empirical p (10,000 matched windows)"; raw, single test per panel |

**New deliverable (Daniil's request): `G3_figure_pvalue_labels_260803.md`** — a standalone reference
document for the manual Figma pass, listing **every statistic in every figure and panel** with:

| Column | Content |
|---|---|
| Figure / panel | e.g. `Figure 4B` |
| Figma frame | node ID from §5b, e.g. `859:25` |
| Statistic | the test that produced the number |
| Source column | the exact CSV column plotted (`Enrichment_p_value_adjusted`, `FDR`, `p_raw_empirical`, …) |
| **Correct label** | the verbatim string to put on the axis/legend/colourbar |
| Raw or corrected | one word, so there is nothing to infer |
| Where it appears | axis label / legend / colourbar / caption / in-panel annotation |

Written so Daniil can work through the Figma frames without re-deriving anything.

Add one sentence to Methods: "Unless stated otherwise, all reported p-values are Benjamini–Hochberg
FDR-adjusted; raw p-values are labelled as such."

**This is journal policy, not only a reviewer preference.** G3's statistics guidance requires that it be
clear whether reported p-values are raw or corrected, that the type of correction be stated, **and that
raw p-values be available in the supporting materials** so readers can apply their own correction
(`G3_article_guidelines.md` §5). Our enrichment CSVs already carry both (`Enrichment_p_value` /
`Enrichment_p_value_adjusted`, `p_raw_empirical` / `p_adjusted_empirical_bh`) and the GO tables carry
both `P-value` and `FDR` — so the obligation is met by *saying so*: extend the Methods sentence to name
the supplementary files that carry the raw values, and repeat it in the Data availability statement.

---

### WP10 — Restructure Table 1 (Minor 4)

**Reviewer:** "Table 1 appears too wide… consider splitting into two tables, reducing column widths,
or moving less critical columns to supplementary materials."

Current: 11 columns × 6 rows. Proposed split:

**Table 1. Enrichment of TE classes in gene TSS 10 kb neighbourhoods.** (5 columns)
| Class | TEs in TSS windows | TEs total | Observed OR | Observed/random OR |

**Table 2. Statistical support for TE class enrichment.** (5 columns)
| Class | Adjusted Fisher p | Random OR (mean ± SD) | Empirical p | Adjusted empirical p |

- Merge "Mean of random OR" and "SD of random OR" into one `mean ± SD` column.
- Drop the raw (unadjusted) Fisher p from the main tables → supplementary (redundant with the adjusted
  one at these magnitudes).
- Renumber the existing Table 1 references in text; there is currently only one other table reference
  to check.
- Keep scientific notation consistent (`<10⁻²⁰⁰` is fine; G3 will typeset it). Fix the malformed
  `9.3*10⁻¹³³` style (WP15/G8).
- **Values do not change** — D1 keeps N = 500, so this is a reformat, not a regeneration.
  `revision_G3/10_tables.py` reads the existing `enrichment_families_with_random.csv` and emits
  `revision_G3/output/Table1.csv` / `Table2.csv`, so the split is reproducible rather than hand-typed.

---

### WP11 — "flavone metabolism" (Minor 5)

Verified: `flavone metabolic process` is **GO:0051552**, a real term, with 5 overlapping LINE-set
genes — not a typo for flavonoid. Two-part response: (i) confirm the accession in the response letter,
(ii) note that the term is at FDR = 0.088 and is therefore removed under the new 0.05 threshold (WP6),
so the sentence disappears from the Abstract and Results. Report the 5 genes in the response letter
(likely CYP/UGT family — confirm from `GO_top_5_perc_genes_by_class_number_with_all.csv`) so the
reviewer can see it was not an annotation error.

---

### WP12 — Simplified network figures, with no overlapping text (Minor 6)

**Reviewer:** "Figures 4A, 5A, and 6A are dense network visualizations that are difficult to interpret
in print. Please provide simplified versions in the main text… and move the full versions to
supplementary materials."

The existing plotting functions already accept most of what is needed — `save_go_network_svg(df, ...,
jaccard_threshold=0.1, top_n=5)` in `Gene_ontology_analysis.ipynb` cell 6, `visualize_go_class_network`
in cell 36, `save_go_network_svg_families_by_classes` in cell 175. Current calls use top 30 terms per
group (≈180 nodes for Figure 4A). **Those functions are copied into `revision_G3/revision_lib.py`** and
extended there with `min_shared_genes`, `max_term_genes` and the collision check; the originals stay
untouched (frozen-notebook rule, §2).

Main-text (simplified) settings — Daniil has agreed to this filtering:
- `top_n = 10` (from 30)
- edge filter: Jaccard ≥ 0.2 **and** ≥ 5 shared genes (add a `min_shared_genes` argument — currently
  only the Jaccard threshold exists)
- exclude terms with > 500 genes (currently > 1,000) — removes "protein binding" (468 overlapping
  genes), "cytoplasm", "nucleus", "cytosol", "nucleoplasm", which carry no interpretation and dominate
  the layout
- label all retained nodes (legible at print size with top_n = 10)
- **FDR < 0.05 only** (WP6)

**No overlapping text — a hard requirement.** Daniil has been fixing label collisions by hand in Figma;
the notebook output must arrive clean. Implement in `nb06_go_networks_fdr005.ipynb`:
- `adjustText` with `force_text`/`force_points` tuned per figure and `expand_text` margins, iterated
  until zero collisions;
- a **programmatic collision check** before writing the SVG — compute every label's rendered bounding
  box (`text.get_window_extent(renderer)`), test all pairs for intersection, and assert zero overlaps
  and zero label–node overlaps. Raise if any remain, so a bad figure cannot be saved silently;
- if collisions persist after adjustment, fall back in this order: shorten GO labels to a curated short
  form (keep a `long_label → short_label` map in `revision_G3/output/go_label_shortnames.csv`), then
  leader lines from node to offset label, then reduce `top_n` for that panel;
- reuse the same collision check for Figures 7A–7H and the WP5 concordance panels.

**Supplementary versions:** keep the current layout settings (top 30, Jaccard 0.1, ≤ 1,000 genes) so no
structural information is lost, **but filter at FDR < 0.05** like everything else (D2) — the
supplementary networks and tables must not carry 0.05–0.1 terms. New Figures S9, S10, S11 (full 4A, 5A,
6A), renumbered after WP2's promotion of the old Supplementary Figure 6.

**Figma side (§5b) — Daniil's manual step.** The dense networks are panel A of frames `859:25`
(Figure 4), `861:28` (Figure 5) and `861:33` (Figure 6). Their descendant counts — 16,761 / 12,755 /
**32,007** — are the quantitative confirmation of the reviewer's complaint, and also why the file is
slow and the figure zip is 53 MB. The notebooks write `svg/Fig4A_simplified.svg`,
`svg/Fig5A_simplified.svg`, `svg/Fig6A_simplified.svg` plus the three full versions; Daniil swaps
panel A inside the existing frames (keeping frame IDs stable) and creates new frames for S9–S11 so the
originals remain available.

`genes_subfamilies_network.py` (`JACCARD_THRESHOLD = 0.025`) belongs to the companion subfamilies
manuscript; apply the same treatment there only when that paper is refreshed (C3).

---

### WP13 — Data availability and a genome browser instance (Editor)

**13a — fix the links.**
- Manuscript line 371: `https://github.com/Nikit357/T2T_genes_evolution` →
  `https://github.com/Nikit357/T2T_transposons_genes` (verified as the actual remote).
- Verify every large deliverable is actually *in* the public repo: `Supplementary File 1.csv` (15 MB),
  `TEs_on_genes.csv` (23 MB), `TEs_on_genes_counts_subfamilies.csv` (89 MB — near GitHub's 100 MB hard
  limit; ship it gzipped or split). `goa_human.gaf` stays gitignored (190 MB) with download
  instructions.
- Add the compacted permutation store (`permutation_counts_10kb/`, ~110 MB after WP1b) plus its
  `MANIFEST.json` and the `01c_expand_counts.py` reconstructor — the permutation background becomes
  publicly reproducible for the first time, which is worth stating in the response letter.
- **DANIIL:** Zenodo archival snapshot + DOI (D9), cited alongside the repo URL.
- Add `revision_G3/README.md` and a top-level `REPRODUCE.md` giving the exact run order
  (notebooks → scripts → revision notebooks) and resolved paths, since the existing relative paths are
  inconsistent (documented in `CLAUDE.md`).
- Rewrite the Data Availability statement to name each deliverable and where it lives.

**13b — UCSC track hub** (`revision_G3/12_build_trackhub.sh`, D8)

```
trackhub/
  hub.txt                # hub name, shortLabel, longLabel, email, descriptionUrl
  genomes.txt            # genome hs1  /  trackDb hs1/trackDb.txt
  hs1/
    trackDb.txt          # composite "T2T_TEs" + derived tracks
    TEs_LINE.bb  TEs_SINE.bb  TEs_LTR.bb  TEs_DNA.bb  TEs_SVA.bb  TEs_RC.bb
    TSS_10kb_windows.bb
    genes_TE_top.bb  genes_TE_bottom.bb
    IFNA_domain.bb
```

- Build from `T2T_repeat_masker_processed_sorted.bed` (§3.0): `name` = subfamily, `score` = divergence
  (0–1000, already in that scale), `itemRgb` = the project class palette (`LINE #cc660b`,
  `LTR #70453c`, `SINE #ab1f20`, `DNA #195f90`, `Retroposon #765297`, `RC #238023`) so the hub matches
  the figures. `TSS_10kb_windows.bb` from `T2T_genes.bed`.
- `bedToBigBed` + `hubCheck` + `fetchChromSizes` from
  `https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/` (verified reachable); chrom sizes from the
  existing `chm13.genome`.
- **Split per class** — one combined bigBed of 3.7 M elements would risk exceeding GitHub's 100 MB
  per-file limit; the largest class (SINE, 1.7 M) should land around 35 MB.
- Host on a `gh-pages` branch (GitHub Pages serves HTTP range requests, which UCSC requires; Zenodo
  does not).
- Validate with `hubCheck -checkSettings` and by loading the one-click URL:
  `https://genome.ucsc.edu/cgi-bin/hgTracks?db=hs1&hubUrl=https://nikit357.github.io/T2T_transposons_genes/trackhub/hub.txt&position=chr9:21150692-21370055`
  (defaulting to the IFNA domain is a nice touch for reviewers).
- Add the hub URL and the one-click link to Data Availability and to the repo README.

---

### WP14 — Response letter, tracked changes, submission package

1. **Rename the working file:** `T2T_genes_article_for_plos_one.docx` →
   `T2T_genes_article_G3_revision_260803.docx` (the current name is misleading and risks the wrong file
   being uploaded). Keep the original untouched as the submitted baseline for the tracked-changes diff.
2. **All manuscript edits are tracked changes** (§3.6), so the clean and tracked versions are the same
   file accepted/rejected rather than two separately maintained documents.
3. **Three files for submission**, as the letter requires: clean version; tracked/highlighted version
   linking each response to the changed text; separate point-by-point response document.
4. Draft the response as `G3_response_to_reviewers_260803.md` → export to docx. One heading per
   reviewer item, in their numbering, each with: *the comment (verbatim)* → *what we did* → *where it is
   in the revised manuscript (section, figure, table, line)*. Quote the new numbers. Three items need
   particular care because they decline part of what was asked:
   - **Major 3** — the proximity-null argument plus the 2026b companion (WP3);
   - **Major 4** — gene-set overlap plus the explicit "full reimplementation of their clustering was
     out of scope" sentence (WP4/D6);
   - **Minor 1** — N = 500 with the bias-correction justification and the convergence figure (WP1).
5. **Final formatting pass to G3 house style — WP15.**
6. Cover letter: note the companion preprints (2026a preprint of this work, DOI 10.32942/X2FM2M; the
   epigenomic companion 2026b, DOI 10.64898/2026.03.19.712972) and the scope boundary between them
   (C4). State that a Reagent Table is not applicable (purely computational study) so its absence does
   not read as an omission.
7. **Renumbering sweep** after WP2/WP12: Figures 1–7, **new 8 (IFNA)**, 9, 10 (schematics);
   supplementary items renamed to G3 convention (`Figure S1`…, `File S1`…, see WP15/G2):
   Figures S1–S5, **S6 (was Supp. 6, now folded into main Fig. 8A — confirm whether to keep a
   standalone version)**, S7, S8, **new S9–S11 (full networks)**, **new S12 (Lu et al. overlap),
   S13 (window/percentile concordance), S14 (permutation convergence)**.
   Do this once, at the end, from a written map — it is the most error-prone step of the revision.
   The sweep runs in **four places**: the manuscript text, the figure captions, the on-disk filenames in
   `manuscript_figures_supplementary_260803/`, and **the Figma frame names (§5b)**. Scope the Figma part
   **by node ID, not by name** — the same canvas carries the subfamilies paper's `Figure 9`–`Figure 14`
   labels, and `Figure 9` already exists twice (frame `861:35` vs. a loose label at 14084, 2019).

---

### WP15 — Final formatting to G3 house style

**Driven by:** `G3_article_guidelines.md` (written 2026-08-03), which is the distilled requirement set
plus an 18-item gap list (G1–G18). Read it before touching the manuscript file.

**Provenance caveat you need to know about.** The live guidelines page
(`https://academic.oup.com/g3journal/pages/author-guidelines?login=false#section-11`) is **not
retrievable from this environment** — OUP returns HTTP 403 to every access route tried (plain curl,
browser and Googlebot user-agents, reader proxies), and the Internet Archive's captures of that URL are
OUP's crawl-prevention placeholder rather than page content. `G3_article_guidelines.md` is therefore
built from the archived GSA "Preparing Manuscripts for Submission" page, the current house style
extracted programmatically from six G3 articles published in 2026 (PMC13365844,
PMC13334165/167/169/170/186), and the decision letter — with every statement tagged by source and a
§13 list of 10 items **Daniil must confirm in a browser**, because the letter warns the guidelines were
recently updated.

**Daniil implements all of WP15 except G1 (citation style), which he will do in Mendeley.**

| Gap | Why it matters | Owner | Effort |
|---|---|---|---|
| **G1 — citation style** | G3 uses **author–year** ("Lu et al. 2020"), not numbered. The manuscript is Vancouver numeric via Mendeley Cite (§3.5). **DANIIL:** switch `MENDELEY_CITATIONS_STYLE` to the G3/CSE author–year style in the Mendeley Cite pane and refresh (the bibliography is already flagged `IS_DIRTY`). | **Daniil** | — |
| **G1b — in-sentence citation grammar** | Several sentences are ungrammatical under author–year: *"the previous landmark study by (14)"*, *"A landmark study by (10) that showed…"*. This text sits **outside** the Mendeley content controls, so it is ours to fix — as tracked changes, after Daniil's style switch. | us | 0.5 d |
| **G2 — supplementary naming** | G3 convention is `Figure S1`, `Table S1`, `File S1`; our draft and on-disk filenames use "Supplementary Figure 1" / "Supplementary File 1". Rename in body text, captions, on-disk files, Figma frame names (by node ID, C9), `CLAUDE.md`, response letter. Do together with WP14.7. | us | 0.5 d |
| **G7 — section order and headings** | Body: Abstract → Keywords → Introduction → *Materials and methods* (sentence case) → Results → Discussion → Supplementary material. Back matter: Acknowledgments → **Data availability** → Funding → Conflicts of interest → (Code availability) → **Literature cited** (not "References"). Our draft has all-caps headings, `REFERENCES`, an `ETHICAL STATEMENT` section G3 does not use, and British `ACKNOWLEDGEMENTS` (consistent with the `en-GB` Mendeley locale — switch that too). | us | 0.5 d |
| **G4/G5/G6 — missing title-page elements** | No keywords (need 3–10), no ~35-character running title, corresponding-author block has only an email (needs institutional address, phone, ORCID 0000-0003-1029-1174). | us | 0.2 d |
| **G11/G12 — figure production** | Figure titles must live in the manuscript legend, not in the image: our matplotlib `fig.suptitle(...)` calls bake titles into the panels (e.g. "Detailed TE Family Enrichment Analysis" on Figure 1D) — remove them. Charts export as vector PDF from the Figma frames (§5b); `.jpg`/`.docx` figures are rejected outright. **Raster check is narrower than first assumed:** the UCSC panels (Figure 1C, Supplementary Figure 6) are vector, so the ≥ 350 dpi rule bites only on the pasted colorbars in Figures 4/5/6 + Supplementary Figure 8 and on the schematics (Figure 8: 3 bitmaps; Figure 9: 10 bitmaps, two reused from the PhD thesis). | us + Daniil | 0.5 d |
| **G17 — baked-in FDR text** | `"GO terms count in a group (FDR < 0.1)"` in Figure 7 (frame `861:34`) must become `0.05`; grep all frames for `0.1` / `500` / `p-value` (C11). | Daniil (Figma) | 0.1 d |
| **G18 — mixed fonts and label sizes** | Inter throughout except Helvetica in the vector UCSC panels; axis-label sizes vary (10 / 13.3 / 14.7 / 16 px). Unify while the frames are open. | Daniil (Figma) | 0.2 d |

**Convenient alignments** (nothing to do): Investigations have no length limit; our abstract is 205 words
against a 250-word cap; `GLOBAL_FONT_SIZE = 10` already matches the sans-serif 10 pt requirement for
figure labels; our enrichment/GO tables already carry raw *and* FDR-adjusted p-values, which is what the
journal's statistics policy demands (see WP9).

**Also note:** G3's statistics policy independently requires the raw-vs-adjusted labelling that
Reviewer 1 asks for in minor comment 3 — cite that alignment in the response letter, it makes the answer
stronger.

**Sequencing.** WP15 runs **last**, after every number has settled (WP5, WP6) and after the figure
renumbering map exists (WP14.7). G1b (in-sentence citation grammar) must come after Daniil's Mendeley
style switch and after the Discussion trim (WP8), or the sweep gets done twice: WP8 → WP14.7 →
**Daniil: Mendeley style switch** → G1b → rest of WP15 → tracked-changes verification.

---

### WP16 — Notebook-based figure production surface (cross-cutting)

**Daniil's constraints, combined:** *"create new notebooks in a new `revision_G3/` subdirectory that
return exactly those plots that reviewer asked to modify or add"*; *"Do not write scripts that will
export figures to figma, just save them as SVG and I will add them manually myself"*; *"implement
figures generation in jupyter notebooks for me to clearly see and edit the subfigures before saving to
figma"*; *"create completely new notebooks in the `revision_G3` folder, do not change the currently
existing ones. The new notebooks should cover only the changes requested in this plan"*.

**Consequences for the whole plan:**

1. **The existing notebooks are frozen — zero edits.** `TEs_mapped_on_TSS_analysis.ipynb`,
   `Gene_ontology_analysis.ipynb` and `download_and_process_files_UCSC_genes.ipynb` are read-only
   references from here on. Record their current MD5s in `revision_G3/README.md` as the freeze baseline
   (all three currently have uncommitted working-tree changes, so **commit them first** — otherwise
   "unchanged" is unverifiable):

   ```
   6d59a2a735b8d0f4fcf6d9dddbb8bb39  TEs_mapped_on_TSS_analysis.ipynb
   a75ceaf51c0a0d221f53357bb0040b55  Gene_ontology_analysis.ipynb
   3e8aec87bd9e78fce53463a2073d968b  download_and_process_files_UCSC_genes.ipynb
   cfd78a7eb38b8f5bbc76dd0fba75dc01  GO_subfamilies.py
   ```

2. **Every revision plot lives in a new notebook**, not a script. Computation stays in `.py`/`.sh` (they
   are cacheable and re-runnable); the notebook reads the computed CSVs from `revision_G3/output/` and
   renders. This satisfies both the modularisation rule in `CLAUDE.md` and Daniil's need to inspect each
   subfigure inline.
3. **Shared code lives in `revision_G3/revision_lib.py`**, populated by **copying** the helpers out of
   the frozen notebooks' cells and modifying the copies: the two `run_goatools_*` functions
   (FDR 0.05), the three network plotters (plus `min_shared_genes`, `max_term_genes`, the collision
   check), `load_permutation_counts()` for the compacted store, the TE class palette, and the font
   settings. The module carries a header comment naming the notebook and cell each function came from,
   so the provenance of every copy is explicit (caveat C19).
4. **Minimal scope: a figure that does not change is not regenerated.** Figures 1A–1C, 2A–2C, and
   Supplementary Figures 4–7 are untouched — they keep their current SVGs. Only panels that the reviewer
   asked to change, that carry a wrong statistics label, that report GO results at 0.1, or that are new
   appear in the notebooks.
5. **Notebooks emit SVG only**, into `revision_G3/svg/`, one file per panel, named for its destination:
   `Fig1D_relabelled.svg`, `Fig4A_simplified.svg`, `Fig8B_null_distributions.svg`,
   `S13_window_concordance.svg`, … `revision_G3/svg/PLACEMENT.md` maps every SVG to its target Figma
   frame ID and panel position, so the manual step needs no guesswork.
6. **Nothing writes to Figma.** No `use_figma` calls, no automated frame edits, no exports. §5b is a
   read-only reference for Daniil's manual assembly.

**Five notebooks, and only five:**

| Notebook | Covers | Produces |
|---|---|---|
| `nb01_permutation_convergence.ipynb` | WP1 | running mean/SD of the random OR per class and family → S14 (the N = 500 justification figure) |
| `nb02_ifna_domain.ipynb` | WP2 | Figure 8B (null distributions), 8C (subfamily composition + per-element divergence); displays all four test results inline for review |
| `nb03_relabelled_figures.ipynb` | WP9 | re-creates the panels whose only defect is a statistics label, from unchanged data: Figure 1D (copied from frozen cell 63), 2D–2F, 3A/3B, Supplementary 1–3. **No numbers change** — this notebook exists purely to satisfy Minor 3 without editing the frozen notebooks |
| `nb05_sensitivity_robustness.ipynb` | WP5 | window/percentile concordance panels → S13; the full robustness/intersection comparison against the original 10 kb + 5 % results |
| `nb06_go_networks_fdr005.ipynb` | WP6, WP12 | simplified 4A/5A/6A + full S9–S11, all at FDR 0.05, all collision-checked; regenerated 4B/5B/6B, 7A–7H, S8B with recomputed GO counts; vector colorbars |

Each notebook: `plt.rcParams["svg.fonttype"] = "none"`, `GLOBAL_FONT_SIZE = 10`, the shared TE class
palette (all from `revision_lib.py`), and a final cell listing every SVG written with its target frame —
so the handoff to Daniil is explicit.

---

## 5. New files to create

```
revision_G3/
  README.md                          # run order, inputs, outputs, runtimes + frozen-notebook MD5 baseline
  revision_lib.py                    # helpers COPIED out of the frozen notebooks and modified here:
                                     #   run_goatools_enrichment / _ordered (fdr=0.05)
                                     #   save_go_network_svg, visualize_go_class_network,
                                     #   save_go_network_svg_families_by_classes
                                     #     (+ min_shared_genes, max_term_genes, collision check)
                                     #   load_permutation_counts(window=...)  -> compacted store reader
                                     #   TE class palette, font settings, assert_no_label_collisions()
  # --- compute (scripts: cacheable, re-runnable) ---
  01_permutations_stream.sh          # WP1/WP5: N=500 streaming permutations per window (5 kb, 20 kb)
  01a_consolidate_counts.py          # WP1: consolidated counts + legacy N=500 regression check
  01b_compact_permutation_results.py # WP1b: 11 GB -> ~110 MB lossless, verified per seed
  01c_expand_counts.py               # WP1b: reconstruct a legacy BED from counts on demand
  02_ifna_domain_test.py             # WP2: 4 tests -> output/ifna_*.csv
  04a_lu2020_geneset_overlap.py      # WP4: overlap with Lu et al. 2020 categories
  04b_newly_resolved_regions.py      # WP4/C6: upper bound on the assembly-attributable component
                                     #   (added during implementation; NOT the cut 04b_hg38_rerun.sh)
  05a_build_windows.sh               # WP5: 5 kb / 20 kb TSS neighbourhoods + TE mapping
  05b_window_sensitivity.py          # WP5: enrichment at 5/10/20 kb, 6 classes + 44 families
  05c_percentile_sensitivity.py      # WP5: GO at top/bottom 10 % vs 5 %
  06_go_rerun_fdr005.py              # WP6: all GO levels at FDR 0.05
  10_tables.py                       # WP10: Table1.csv + Table2.csv (reformat, values unchanged)
  11_results_numbers.py              # Phase 5: re-derives every number the revised text quotes
                                     #   -> output/results_numbers.{json,txt} (added during
                                     #   implementation, so nothing is typed from memory)
  12_build_trackhub.sh               # WP13b: bigBeds + hub.txt/genomes.txt/trackDb.txt + hubCheck
                                     #   (renumbered from 13_ to keep the manuscript scripts last)
  # --- manuscript (Phase 5; run in the collagen venv, which has python-docx) ---
  13_manuscript_tracked_edits.py     # Phase 5 D-K: baseline docx -> revision docx, all edits
                                     #   tracked; idempotent because it always starts from the
                                     #   baseline. Citation-safe: never rewrites a <w:sdt>
  14_build_extensive_discussion.py   # WP8/D7: Extensive_discussion_260803.docx by copy-and-delete
  15_house_style.py                  # Phase 7: G3 house style (G2-G16); MUST run after 13_
                                     #   from the baseline, so the moved citations stay live
  # --- figures (NEW notebooks only; the existing ones are frozen) ---
  nb01_permutation_convergence.ipynb # WP1 justification figure
  nb02_ifna_domain.ipynb             # WP2 figures + test results for review
  nb03_relabelled_figures.ipynb      # WP9 label-only re-creations: 1D, 2D-2F, 3A/3B, S1-S3
  nb05_sensitivity_robustness.ipynb  # WP5 figures + robustness/intersection comparison
  nb06_go_networks_fdr005.ipynb      # WP6/WP12 networks and heatmaps at 0.05, collision-checked
  external/                          # WP4 third-party inputs; the only files here we do not produce
    PROVENANCE.md                    #   URLs, sizes, md5s, and the download routes that no longer work
    lu2020/mmc2..7.xlsx              #   Lu et al. 2020 supplementary tables (Table S1 = mmc2)
    HOM_MouseHumanSequence.rpt       #   MGI mouse-human homology (gitignored, 15 MB, re-downloadable)
    hs1ToHg38.over.chain.gz          #   UCSC T2T -> GRCh38 liftOver chains
  output/                            # all derived tables (gitignored if > 50 MB)
    permutation_counts_10kb/         # WP1b compacted store + MANIFEST.json
    legacy_fdr01_n500/               # C3 snapshot of pre-revision CSVs
  svg/                               # every panel SVG for Daniil's manual Figma placement
    PLACEMENT.md                     # SVG -> Figma frame ID + panel position
```

Repo root:

```
G3_response_to_reviewers_260803.md   # WP14
G3_figure_pvalue_labels_260803.md    # WP9: which p-values are raw vs FDR, per figure/panel/frame
revision_G3/Revised_manuscript/       # the three manuscript docx files, kept together
  T2T_genes_article_G3_submitted_baseline_260418.docx   # read-only; what the journal received
  T2T_genes_article_G3_revision_260803.docx             # the revision, all edits tracked
  Extensive_discussion_260803.docx     # WP8/D7: built by copying the manuscript, not by pasting
REPRODUCE.md                         # WP13a
trackhub/                            # WP13b deliverable (published via gh-pages)
```

Already written (companion documents, repo root):

```
G3_article_guidelines.md             # WP15: G3 house style + gap list G1-G18 + browser-verify list
G3_revision_implementation_plan_260803.md   # this file
G3_reviewer_report_260802.md         # the decision letter (input)
```

Naming/style: `plt.rcParams["svg.fonttype"] = "none"`, `GLOBAL_FONT_SIZE = 10`, the shared TE class
palette, SVG output — per the project `CLAUDE.md`. Every script writes its intermediate DataFrame to
`revision_G3/output/` (persistence rule).

**Cut from the original plan** (D4, D4a, D6, D1): `03a_gtex_tissue_specificity.py`,
`03b_encode_ifna_chromatin.py`, `04b_hg38_rerun.sh`, `04c_hg38_vs_t2t_compare.py`,
`01b_permutation_convergence.py` (became `nb01`), `12_networks_simplified.py` (became `nb06`).

**No file outside `revision_G3/` and the repo-root list above is created or modified**, except the
manuscript, the three documentation files in §7.5–7.7, and the `trackhub/` deliverable. In particular no
existing notebook or existing `.py` script is touched (§7.2–7.4, §7.8).

---

## 5b. Figure production surface — the Figma file (read-only reference)

Panels are exported from the notebooks as SVG; **Daniil composes, letters and exports the final figures
in Figma by hand** (D-level instruction). This section exists so that manual step needs no
investigation: it records where each figure lives, what is raster, and what text has analysis
parameters baked into it. **Nothing in this plan writes to Figma.**

**File:** `Figures T2T genes` — `https://www.figma.com/design/WRNeTzKZObdmAQ8QG1EZlq/Figures-T2T-genes`
**File key:** `WRNeTzKZObdmAQ8QG1EZlq` · **single page:** `0:1 "Page 1"`
Deep-link any figure with `…/Figures-T2T-genes?node-id=<id with ':' replaced by '-'>`
(e.g. Figure 1 → `?node-id=856-7`).

### 5b.1 Node-ID map (verified 2026-08-03 by read-only inspection)

Every manuscript figure exists as a named top-level `FRAME`. These are the authoritative artboards —
export from these, not from the loose groups around them.

| Manuscript figure | Figma frame | Size (px) | Descendants | Notes |
|---|---|---|---|---|
| Figure 1 | `856:7` | 1210 × 1854 | 4,901 | contains the CIB3 UCSC panel (1C) — **vector**, Helvetica text |
| Figure 2 | `856:8` | 1204 × 1363 | 2,978 | |
| Figure 3 | `856:9` | 1194 × 915 | 5,411 | KS-test annotations live as text nodes here |
| Figure 4 | `859:25` | 1284 × 1768 | 16,761 | GO network 4A — **simplify (WP12)**; 2 raster colorbars |
| Figure 5 | `861:28` | 1222 × 1749 | 12,755 | GO network 5A — **simplify (WP12)**; 2 raster colorbars |
| Figure 6 | `861:33` | 1238 × 1956 | **32,007** | GO network 6A — **simplify (WP12)**; 2 raster colorbars; already carries a "-log (FDR adjusted p-value)" label |
| Figure 7 | `861:34` | 1216 × 1097 | 1,487 | **contains the literal string "GO terms count in a group (FDR < 0.1)" → must become 0.05 (WP6/G17)** |
| Figure 8 (schematic) | `766:1208` | 630 × 591 | 367 | mechanisms schematic; 3 raster images. **Becomes Figure 9** |
| Figure 9 (schematic) | `861:35` | 733 × 636 | 96 | "integrated map"/Ring of Power; **10 raster images**, two named "New Figures for Disser 1" and "Danya Nikitin PhD theme overview 1" (reused from the PhD thesis). **Becomes Figure 10** |
| Supplementary Figure 1 | `856:10` | 1141 × 610 | 3,914 | |
| Supplementary Figure 2 | `856:11` | 1198 × 1902 | 24,552 | |
| Supplementary Figure 3 | `856:12` | 1198 × 1902 | 24,918 | |
| Supplementary Figure 4 | `860:26` | 1182 × 1003 | 5,552 | |
| Supplementary Figure 5 | `860:27` | 1120 × 1122 | 6,338 | |
| Supplementary Figure 6 | `861:29` | 1155 × 691 | 866 | IFNA-domain UCSC view — **vector, no raster fills**; **promoted to main Figure 8 panel A by WP2** |
| Supplementary Figure 7 | `861:30` | 861 × 1828 | 3,540 | |
| Supplementary Figure 8 | `861:32` | 1229 × 1871 | 6,516 | 1 raster colorbar |

### 5b.2 What the inspection established

1. **The dpi worry was wrong for the browser panels.** Figure 1C and Supplementary Figure 6 contain
   **no raster fills** — they are vector UCSC exports. The real raster exposure is elsewhere: small
   pasted **colorbars** (23 × 58 to 30 × 175 px, `scaleMode=CROP`) inside Figures 4, 5, 6 and
   Supplementary Figure 8, and the **schematics** — Figure 8 (3 images) and especially Figure 9 (10
   images). Those are what must clear the ≥ 350 dpi rule; the colorbars should ideally be re-emitted as
   vector from matplotlib rather than pasted as bitmaps (the WP12 notebook can do that).
2. **Figure text encodes analysis parameters** and will silently contradict the revised Methods unless
   edited: `"GO terms count in a group (FDR < 0.1)"` in Figure 7, and any sibling "FDR < 0.1" wording.
   Grep the frames for `0.1`, `500`, `p-value` after WP6 lands (C11, G17).
3. **The canvas is shared with the subfamilies manuscript.** Below the framed figure row
   (y ≈ −2,400 … +400) sits a working band (y ≈ +1,700 … +4,700) of 154 loose top-level groups plus
   label texts reading `Figure 9`–`Figure 14` and `Supplementary Figure 9`–`Supplementary Figure 12`.
   Those belong to the companion paper. **`Figure 9` exists twice** — as the G3 schematic frame
   (`861:35`) and as a loose label at (14084, 2019) in the subfamilies band. Any renumbering must be
   scoped to the framed row by node ID, never by name search (C9).
4. **Font consistency:** Inter throughout, except **Helvetica** in Figure 1 and Supplementary Figure 6
   (inherited from the UCSC exports). Axis-label sizes vary across figures (10 / 13.3 / 14.7 / 16 px).
   Neither is fatal, but both are visible in print and easy to fix while the figures are open (G18).
5. **Housekeeping:** three empty leftover frames (`422:32` 9 × 17, `260:206768` 32 × 43, `588:1168`
   0 × 0) sit at top level and would export as blank artboards in a bulk export. Delete or skip
   explicitly.

### 5b.3 Manual figure workflow (Daniil)

For each changed figure: **(1)** the notebook writes the panel SVG to `revision_G3/svg/` and records
its target in `svg/PLACEMENT.md` → **(2)** place that SVG into the Figma frame named in the map, keeping
the frame ID stable so links and this table stay valid → **(3)** check panel letters (G3 requires a
letter in the upper-left of each panel) → **(4)** export the frame as **PDF** (vector) plus PNG for
preview → **(5)** copy into `manuscript_figures_supplementary_260803/` → **(6)** re-check the caption
text in the manuscript.

Do **not** re-create frames from scratch: the node IDs in 5b.1 are referenced by this plan, by
`G3_article_guidelines.md`, by `G3_figure_pvalue_labels_260803.md`, and by the response letter's figure
references.

---

## 6. Compute and storage plan (read before starting)

| Job | Output size | Wall clock (8 cores) | Notes |
|---|---|---|---|
| **WP1b compaction of `permutation_results/`** | 11 GB → **~110 MB** | ~1–2 h | **run first**; per-seed verify before each delete; frees ~10.9 GB |
| WP5 5 kb permutations, N = 500 | ~20 MB compressed | 2–4 h | streaming counts; background |
| WP5 20 kb permutations, N = 500 | ~20 MB compressed | 2–4 h | streaming counts; background |
| WP2 IFNA null windows (10,000 × 3 matching schemes) | < 50 MB | ~20 min | `bedtools shuffle -incl` |
| WP6 GO re-run at FDR 0.05 (all levels) | ~200 MB (`GO_tables/` 230 files) | ~1–2 h | overwrite after the C3 snapshot |
| WP13b bigBeds | ~100–150 MB total | ~1 h | split per class, ≤ 100 MB per file |

**Total new compute ≈ 6–10 h**, all backgroundable, and none of it on the critical path. D1 removed the
1,000-permutation re-run; D6 removed the hg38 re-run (which was 6–10 h on its own plus a large
download); D4 removed the GTEx download.

**Free disk is 2.4 GB, and WP1b is what fixes that** — after compaction there is ~13 GB free, which
comfortably covers everything above. Do not start WP5 before WP1b finishes.

---

## 7. Existing files to change

### 7.1 `T2T_genes_subfamilies_article_figures/T2T_genes_article_for_plos_one.docx` (→ renamed and relocated)

Now `revision_G3/Revised_manuscript/T2T_genes_article_G3_revision_260803.docx`. All three manuscript
docx files — the read-only submitted baseline, the tracked revision and the extended-discussion
supplementary file — sit together in `revision_G3/Revised_manuscript/`, away from the subfamilies
paper they were originally filed with. `13_manuscript_tracked_edits.py` and
`14_build_extensive_discussion.py` resolve that directory relative to themselves, so the move needs
no configuration.

The manuscript itself. All edits are **tracked changes** (§3.6) and must move Mendeley `<w:sdt>`
content controls as whole units (§3.5).

#### 7.1a Order of manuscript edits (Daniil's request: longest-lead first)

Compute and analysis come first because they set every number the text quotes. Nothing in the manuscript
is touched until stage C.

| Stage | Work | Why here |
|---|---|---|
| **A. Free the disk** | WP1b compaction | Blocks everything else; nothing fits until it is done |
| **B. Long compute, in parallel with nothing** | WP5 5 kb + 20 kb permutations (2–4 h each); WP6 GO re-run at 0.05; WP2 IFNA null windows | The longest jobs, and every downstream number depends on them. Launch A → B on day 1 |
| **C. Analysis and figures** | WP5 sensitivity + robustness notebook; WP6 outputs; WP2 tests; WP4a overlap; WP12 networks; WP1 convergence figure | Produces the numbers and panels the text quotes |
| **D. Methods first** | N = 500 justification; FDR 0.05; window rationale moved from Discussion; 5/20 kb + 10 % sensitivity methods; the "all p-values are FDR-adjusted unless stated" sentence; §3.0 input files named | Everything else references Methods; editing it first prevents contradictions |
| **E. Tables** | Table 1 → Tables 1 + 2 (WP10) | Referenced from Results |
| **F. Results — numbers** | Recomputed GO counts (lines 93–141); "14 families" → 13; hAT-Charlie three → one; MIR zinc-only; flavone removed (lines 93, 121) | Needs C complete |
| **G. Results — new text** | IFNA descriptive paragraph (§3.3) + four test results + the two divergence quantities; sensitivity/robustness paragraph; Lu et al. overlap paragraph | Needs C complete |
| **H. Abstract and title** | D3 title; causal language; flavone out; MIR wording; new IFNA statistics | Written last among the front matter, once Results are final |
| **I. Introduction** | Remove "arms race" as a claim about our data (WP7) | Short, independent |
| **J. Discussion** | Restructure to 6 subsections ~2,200 w; keep *3.5 cancer*; new Limitations; new Future testing; new proximity-null subsection; build `Extensive_discussion_260803.docx` | Largest text job; must come after F/G so it discusses final numbers |
| **K. Data availability** | Corrected repo URL, compacted permutation store, track hub + one-click link, per-deliverable listing, raw-p-value sentence | Needs WP13 deliverables to exist |
| **L. Renumbering sweep** | Figures/supplements in text, captions, filenames, Figma frame names — once, from a written map | Must follow every figure addition |
| **M. Daniil: Mendeley style switch** | Vancouver numeric → G3/CSE author–year; refresh bibliography | Blocks G1b |
| **N. House style** | G1b in-sentence citation grammar, then the rest of WP15 (G2, G4–G7, G11/G12, G17, G18) | Last, per C12 |
| **O. Package** | Response letter, tracked/clean versions, cover letter, upload | Final |

Regenerate the plain-text view after each editing round:

```bash
pandoc -t plain --wrap=none revision_G3/Revised_manuscript/T2T_genes_article_G3_revision_260803.docx -o /tmp/ms.txt
```

Highest-risk text edits (numbers that must be right): Methods "1000 random permutations" (line 277) →
**500**; keep "empirical p-value = 0.004" (line 37) **unchanged**; Methods GO threshold (line 287) →
0.05; every GO term count in Results lines 93–141; the family count "14 families" (line 119) → 13;
hAT-Charlie "three GO terms" (line 127) → one; MIR metals (lines 125, 223) → zinc only; flavone
(lines 93, 121, 217) removed; the IFNA divergence range "95 – 161.7" (line 103, clarify what it
averages and add the 135.7 element-level mean); Data Availability (line 371).

### 7.2 `TEs_mapped_on_TSS_analysis.ipynb` — **FROZEN, no edits**

Not modified. What each previously-planned cell edit becomes instead:

| Frozen cell | What it does | Replacement |
|---|---|---|
| 40 (permutation config) | `NUM_PERMUTATIONS = 500`, reads the 6.37 GB consolidated CSV | Already correct at 500 (D1), so nothing to fix. The compacted-store reader becomes `load_permutation_counts()` in `revision_lib.py`; the frozen cell stops working once the CSV is deleted, which is fine because its outputs are final — `REPRODUCE.md` documents `01c_expand_counts.py` for anyone who needs to re-run it |
| 63 (Figure 1D) | wrong axis/title labels on an already-FDR-corrected column | code copied into `nb03_relabelled_figures.ipynb`, labels fixed there (WP9) |
| 55, 60, 82 (enrichment) | produce `enrichment_*_with_random*.csv` | **nothing to do** — D1 keeps N = 500, so these outputs stand as published |
| 105, 171–175, 181, 194 | Figures 2D–2F, 3A/3B, Supplementary 1–3 | code copied into `nb03`, statistics provenance labels added (WP9) |
| the `../T2T_article/…csv` path | broken as written | not fixed in place; the new notebooks and scripts use the local `T2T_repeat_masker_processed_sorted.bed` (§3.0). `CLAUDE.md` records that the frozen notebook still contains the broken path |

### 7.3 `Gene_ontology_analysis.ipynb` — **FROZEN, no edits**

Not modified. Its helper functions are **copied** into `revision_G3/revision_lib.py` and changed there:

| Frozen cell | Function / output | Replacement in `revision_lib.py` / notebooks |
|---|---|---|
| 6 | `run_goatools_enrichment`, `run_goatools_ordered_enrichment` (`fdr_threshold=0.1`), `save_go_network_svg` | copied with `fdr_threshold=0.05` (D2) and the new `min_shared_genes` / `max_term_genes` / collision-check arguments |
| 36, 175 | `visualize_go_class_network`, `save_go_network_svg_families_by_classes` | copied and extended the same way |
| 38, 112, 176 | the network call sites | re-created in `nb06` with `top_n=10, jaccard_threshold=0.2, min_shared_genes=5, max_term_genes=500` for main text; current layout settings but FDR 0.05 for supplementary (WP12) |
| 70, 118, 191, 211, 318, 320, 322 | functional-group heatmaps | re-created in `nb06` at FDR 0.05 with corrected colourbar/legend wording (WP9) |
| 103, 104 | two duplicate throwaway IFNA one-liners | superseded by `02_ifna_domain_test.py` + `nb02`; the frozen cells are simply left as they are |
| 196–202, 367–369 | GO-count correlations | re-created in `nb06`; counts change under FDR 0.05 |

### 7.4 `GO_subfamilies.py` — **not modified**

`fdr_threshold=0.1` at line 47 stays. This script belongs to the companion subfamilies manuscript, whose
refresh is deferred (C3), so changing it now would invalidate that paper's current figures without
regenerating them. Subfamily-level GO at 0.05, where the G3 revision needs it, comes from
`revision_G3/06_go_rerun_fdr005.py --level subfamilies`, which writes to `revision_G3/output/` and leaves
the shared `GO_tables/` directory untouched. Note this in `CLAUDE.md` so the discrepancy is documented
rather than surprising.

### 7.5 `CLAUDE.md` (this directory)

- **Record the notebook freeze.** State that `TEs_mapped_on_TSS_analysis.ipynb`,
  `Gene_ontology_analysis.ipynb`, `download_and_process_files_UCSC_genes.ipynb` and `GO_subfamilies.py`
  are **frozen for the G3 revision**, that `revision_G3/` is authoritative for everything the revision
  changes, and that `revision_G3/revision_lib.py` holds modified **copies** of the old notebooks' helper
  functions (list which function came from which cell — caveat C19). Include the freeze MD5s.
- **Resolve the permutation-count inconsistency:** state that **N = 500 is authoritative**, that the
  generator's `NUM_PERMUTATIONS = 1000` (download notebook cell 34) was the error, that it **remains in
  the frozen notebook and was never executed**, and that the manuscript Methods were corrected down to
  500 — not the other way round (caveat C20).
- **Record that `GO_subfamilies.py` stays at FDR 0.1** while the G3 revision uses 0.05 via
  `revision_G3/06_go_rerun_fdr005.py`, and that the companion paper's refresh is deferred (C3).
- **Record that the frozen `TEs_mapped_on_TSS_analysis.ipynb` still contains the broken
  `../T2T_article/…` path** and that the working replacement is the local
  `T2T_repeat_masker_processed_sorted.bed` (§3.0) — the row in "Working Directory Gotchas" is therefore
  *superseded*, not *fixed*.
- Document the **compacted permutation store**: location, schema, ~100× ratio, `MANIFEST.json`, and
  `01c_expand_counts.py` for reconstruction. Mark `permutation_results/` as removed after verification.
- **Add §3.0's canonical input files** to the "Working Directory Gotchas" table: `T2T_genes.bed` and
  `T2T_repeat_masker_processed_sorted.bed` are present locally and are interval-identical /
  data-identical to the previously referenced files. Mark the `../T2T_article/` path row as **resolved**.
- Add a `revision_G3/` section to "Running the Analysis" and "Key Files", distinguishing compute
  scripts from figure notebooks.
- Update "Key Parameters to Tune": GO FDR 0.1 → **0.05 everywhere, no suggestive band**; the new
  network parameters (`min_shared_genes`, `max_term_genes`, collision check); the 5/10/20 kb window set;
  the 5/10 % percentile set.
- Add the track hub as a deliverable.
- **Record the Figma figure source** (file key `WRNeTzKZObdmAQ8QG1EZlq`, single page `0:1`,
  frame-per-figure map) and state that Python/notebooks produce **panels only, as SVG** — the repository
  currently documents the plotting scripts as if they produced final figures. Same note belongs in
  `T2T_genes_subfamilies_article_figures/CLAUDE.md`.
- **Record the docx/Mendeley facts** (§3.5) and the tracked-changes toolchain (§3.6), including the
  `python-docx` venv gotcha (C15).

### 7.6 `README.md` (repo front page)

- Add the G3 acceptance/citation, the corrected repository URL, the Zenodo DOI (once Daniil mints it),
  the track hub one-click link, and a pointer to `REPRODUCE.md`.

### 7.7 `T2T_genes_subfamilies_article_figures/CLAUDE.md`

- Note that the G3 manuscript is `T2T_genes_article_G3_revision_260803.docx` (renamed from
  `…for_plos_one.docx`) and that the figure inventory in that file describes the **subfamilies**
  paper's figures, not the G3 paper's — currently ambiguous.
- Note the shared Figma canvas and the `Figure 9` name collision (C9).

### 7.8 `download_and_process_files_UCSC_genes.ipynb` — **FROZEN, no edits**

Not modified. Cell 34's `NUM_PERMUTATIONS = 1000` is the original source of the manuscript's
inconsistency and it **stays as it is**.

**This has a consequence you should accept deliberately (caveat C20):** the repository is cited in the
Data availability statement, so a reader can open this notebook and find `1000` while the corrected paper
says `500`. Since the notebook may not be edited, the discrepancy is handled by documenting it in three
places — `revision_G3/README.md`, `CLAUDE.md`, and `REPRODUCE.md` — each stating plainly that N = 500 was
the number actually run (501 files in `permutation_results/`, empirical p floor 2/501 = 0.004), that the
generator's `1000` was never executed, and that `revision_G3/01_permutations_stream.sh` supersedes the
BED-writing path. If you would rather have the notebook itself corrected, that is a one-line change and
the only exception I would argue for — say so and I will make it.

---

## 8. Files that do NOT need to change

| File | Why |
|---|---|
| **`TEs_mapped_on_TSS_analysis.ipynb`** | **Frozen — zero edits** (§7.2). Every planned cell edit is re-implemented in `revision_G3/nb03_relabelled_figures.ipynb` or `revision_lib.py`. Verify with the MD5 in §10. |
| **`Gene_ontology_analysis.ipynb`** | **Frozen — zero edits** (§7.3). Its helper functions are copied into `revision_lib.py` and modified there; its call sites are re-created in `nb06`. |
| **`download_and_process_files_UCSC_genes.ipynb`** | **Frozen — zero edits** (§7.8). Cell 34's `NUM_PERMUTATIONS = 1000` stays; the discrepancy is documented instead (C20). |
| **`GO_subfamilies.py`** | **Not modified** (§7.4). Stays at FDR 0.1 for the companion paper; the revision gets 0.05 from `06_go_rerun_fdr005.py --level subfamilies`. |
| `GO_tables/` (230 files) | Left at FDR 0.1, since `GO_subfamilies.py` is not re-run. The revision's subfamily GO output goes to `revision_G3/output/` instead. |
| `enrichment_families_with_random.csv`, `enrichment_subfamilies_with_random*.csv` | **Unchanged under D1.** The background stays at N = 500, so every OR, random OR, obs/random ratio and empirical p stays as published. Table 1/2 and Figure 1D are a reformat and a relabel, not a regeneration. |
| `T2T_genes_subfamilies_article.docx` | The companion subfamilies manuscript is not under review at G3. Its *GO numbers* change (shared FDR threshold), so it needs a data refresh before its own submission — tracked as caveat C3, not part of this revision. |
| `genes_subfamilies_network.py`, `genes_subfamilies_network_clusters.py` | Produce Figures 6 and 7 of the *subfamilies* paper only. Not cited by the G3 manuscript. |
| `draw_length_divergence_corr.py` | Side analysis, not in the G3 manuscript. |
| `go-basic.obo`, `goa_human.gaf` | Ontology/annotation snapshots — must stay frozen at the Dec 31 2025 version cited in Methods, otherwise the FDR comparison in §3.2 is not reproducible. **Do not update them.** |
| `TEs_mapped_on_TSS_analysis-Copy1.ipynb` | Stale 29 MB duplicate. Delete only with your confirmation; irrelevant to the revision. |
| `individuals_by_classes_TE.csv`, `families_by_classes_TE.csv`, `sequence_report.tsv` | Static mapping/metadata inputs. |
| Supplementary Files 1–8 (root) | Structure unchanged. **Values change only where GO FDR bites** (the GO-derived files), not from permutations. The 5/20 kb and 10 % results become *new* supplementary files, not replacements. |
| `../epigenomic_files/*.mapped_on_TSS_with_TEs.bed` (54 files) | D4 removed the epigenomic analysis; these are untouched by this revision. |

---

## 9. Side effects and caveats

- **C1 — Numeric churn is smaller than first planned, and localised.** D1 means **no** number changes
  for permutation reasons: Table 1/2 values, Figure 1D bars, empirical p = 0.004 and all enrichment CSVs
  stand. The only numbers that move are **GO-derived** (term counts, family counts, the specific terms
  in §3.2) and the **new** analyses (WP2, WP5). Budget a numeric consistency pass over Results,
  Discussion and captions after WP5/WP6 — but the blast radius is the GO sections, not the whole paper.
- **C2 — Disk, and the one irreversible step.** WP1b deletes 11 GB after verification. The
  per-seed check must pass before each delete and the aggregate check before `consolidated_random_data.csv`
  goes; keep `MANIFEST.json` and `01c_expand_counts.py` so the legacy format is reconstructable. If any
  check fails, stop and fall back to `zstd -19` in place (~6×). **Never delete a source whose
  verification did not pass.**
- **C3 — Cross-manuscript coupling.** `GO_tables/`, `top_5_perc_genes_by_*.csv` and the GO summary CSVs
  are shared with the subfamilies manuscript. Tightening to FDR 0.05 silently invalidates its current
  figures and stated counts. Snapshot everything into `revision_G3/output/legacy_fdr01_n500/` **before**
  overwriting, and plan a data refresh for that paper before its own submission.
- **C4 — Self-overlap disclosure.** The 2026b preprint is the epigenomic companion and is now load-bearing
  in our answer to Major 3. Cite it in the Discussion and disclose the relationship in the cover letter,
  so the editor reads it as a declared companion rather than undisclosed overlap.
- **C5 — GO annotation freeze.** Keep `goa_human.gaf`/`go-basic.obo` at the December 31 2025 snapshot; a
  silent re-download would change term memberships and make §3.2 wrong.
- **C6 — Major 4 is answered without the controlled experiment (D6).** The reviewer asked a causal
  question ("methodology or assembly?") and the gene-set overlap does not settle it. Mitigation: bound
  the assembly contribution descriptively (TEs and TSS windows in newly resolved regions), attribute the
  remainder to the named methodological differences, and say explicitly that the reimplementation was out
  of scope. **Risk:** a reviewer could ask for the hg38 re-run in a second round. It is ~6–10 h of
  compute and now affordable on disk after WP1b, so it can be added if requested.
- **C7 — Data QC in the IFNA window.** One element is annotated `L1P3` with divergence **0**, implausible
  for that clade and possibly a RepeatMasker artefact. Check before it appears in a figure; if it is an
  artefact, report the trimmed statistics (WP2 already includes a leave-one-out).
- **C8 — 5 kb window vs. the SVA finding.** The SSU72L1–L5 + POLR2A cluster spans 116 kb, so the
  SVA → termination association may weaken at ±2.5 kb. Report it honestly if so; a scale-dependent result
  is a finding, and pre-committing to that framing now avoids the temptation later.
- **C9 — The Figma canvas is shared between both manuscripts (§5b).** One page holds the G3 figure frames
  *and* the subfamilies paper's working band, with colliding names (`Figure 9` twice). Every rename or
  bulk export must be scoped by node ID. A name-based sweep will silently corrupt the companion paper's
  figures, and there is no branch protecting them.
- **C10 — Mendeley content controls are fragile (§3.5).** 128 in-text citations are `<w:sdt>` content
  controls whose payload lives in a webextension part. Editing the visible run without the control, or
  building the Extensive discussion docx by pasting into a fresh document, yields dead plain-text
  citations. Move whole `<w:sdt>` elements; build the second document by **copying** the manuscript.
  Verify in Word (refresh Mendeley, check no citation renders as plain text) before submission.
- **C11 — Analysis parameters are baked into figure text.** Figure 7 hard-codes "FDR < 0.1" as an axis
  label; other frames may hard-code "500 permutations" or p-value wording. These do not update when the
  CSVs do, so they are the most likely source of a figure that contradicts the revised Methods. Grep the
  frames after WP6 (G17).
- **C12 — House style cannot be parallelised with the text work.** G1b (in-sentence citation grammar)
  touches every sentence that cites anything, so it must come after the Discussion trim (WP8), after all
  recomputed numbers, and after Daniil's Mendeley style switch — otherwise the sweep gets done twice.
  Reserve the last 3 days for WP15 + the tracked-changes verification, and do not start it early.
- **C13 — The guidelines document is second-hand.** `G3_article_guidelines.md` is assembled from an
  archived GSA page, six 2026 G3 articles, and the decision letter, because the live OUP page returns 403
  to every automated route. Its §13 lists 10 items to confirm in a browser (data-availability placement,
  figshare vs submission-system supplements, current reference format, figure resolutions, AI-disclosure
  requirement and placement, ORCID, graphical/plain-language abstract, Reagent Table applicability,
  length caps, Code availability section). Confirm those before the final upload.
- **C14 — Figure 9's schematic reuses thesis artwork.** Frame `861:35` contains 10 bitmaps, two named
  "New Figures for Disser 1" and "Danya Nikitin PhD theme overview 1" — carried over from the PhD thesis
  (Zenodo, DOI 10.5281/zenodo.19052416). Self-reuse is fine, but (i) the resolution must still clear
  ≥ 350 dpi, and (ii) if any panel is a substantially unchanged reproduction, say so in the legend and
  cite the thesis.
- **C15 — `python-docx` is in the wrong venv.** The tracked-changes helpers need `python-docx` + `lxml`,
  which are in `~/venvs/collagen_3_11`, not `~/venvs/Retroelements_3_11`. Decide once (use the collagen
  venv for docx work, or install into Retroelements) and record it in `CLAUDE.md`.
- **C16 — Label-collision checking must be enforced, not attempted.** Daniil has been fixing overlaps by
  hand. The WP12 notebook asserts zero label collisions before writing an SVG, so a colliding figure
  cannot be saved silently. If the assertion fires, use the documented fallback ladder rather than
  loosening the check.
- **C17 — Timeline.** ~22 person-days plus 6–10 h of compute against a 30-day window that started
  July 29. The long poles are WP5 (windows/percentiles + robustness), WP8 (Discussion restructure +
  Extensive discussion docx) and WP12 (collision-free networks). **Start WP1b today**, then launch the
  5 kb and 20 kb permutations and the GO re-run so compute proceeds while text work does. If the schedule
  slips, request the extension the editor explicitly offered.
- **C18 — The reviewer never asked for new biology.** Resist scope creep: every new analysis here exists
  to answer a specific numbered comment. Anything else belongs in the subfamilies paper.
- **C19 — Copied helpers will drift from their originals.** The frozen-notebook rule means
  `revision_lib.py` holds *copies* of `run_goatools_enrichment`, `run_goatools_ordered_enrichment`,
  `save_go_network_svg`, `visualize_go_class_network` and `save_go_network_svg_families_by_classes`, with
  a different FDR threshold and extra arguments. Two versions of the same logic will now exist, at
  different thresholds, and the older one still runs. Mitigations: (i) a header comment in
  `revision_lib.py` naming the source notebook and cell for every copied function; (ii) a statement in
  `CLAUDE.md` that for anything G3-related `revision_lib.py` is authoritative; (iii) when the subfamilies
  paper is refreshed (C3), retire the notebook copies in favour of the module rather than editing them.
  **Do not** import from the frozen notebooks (e.g. via `nbimporter`) — that reintroduces the coupling
  the freeze was meant to remove and would silently pick up FDR 0.1.
- **C20 — The repository keeps a visible contradiction, deliberately.** The frozen download notebook
  still sets `NUM_PERMUTATIONS = 1000` while the corrected paper says 500, and the repo is cited in Data
  availability. This is a documented, not fixed, discrepancy (§7.8): the explanation goes in
  `revision_G3/README.md`, `CLAUDE.md` and `REPRODUCE.md`. It is the one place where the freeze costs us
  something a reader could notice, so the wording must be unambiguous — "500 permutations were run
  (501 files, empirical p floor 2/501 = 0.004); the generator's 1000 was never executed."

---

## 10. Verification commands

```bash
source ~/venvs/Retroelements_3_11/bin/activate
cd /home/jovyan/Projects/Retroelements/T2T_genes_article/T2T_transposons_genes

# --- the frozen notebooks are actually untouched (run this at the END, before submission) ---
md5sum -c <<'SUMS'
6d59a2a735b8d0f4fcf6d9dddbb8bb39  TEs_mapped_on_TSS_analysis.ipynb
a75ceaf51c0a0d221f53357bb0040b55  Gene_ontology_analysis.ipynb
3e8aec87bd9e78fce53463a2073d968b  download_and_process_files_UCSC_genes.ipynb
cfd78a7eb38b8f5bbc76dd0fba75dc01  GO_subfamilies.py
SUMS
# expect 4x OK. If they were committed in Phase 0 (recommended), this is equivalent:
git status --short -- TEs_mapped_on_TSS_analysis.ipynb Gene_ontology_analysis.ipynb \
  download_and_process_files_UCSC_genes.ipynb GO_subfamilies.py    # expect empty output

# --- the new notebooks do not import from the frozen ones (C19) ---
grep -rn "nbimporter\|import_ipynb\|Gene_ontology_analysis\|TEs_mapped_on_TSS_analysis" \
  revision_G3/*.ipynb revision_G3/*.py | grep -v "^.*#"   # expect only comments/provenance notes

# --- §3.0: the canonical inputs are what we think they are ---
wc -l T2T_repeat_masker_processed_sorted.bed          # expect 3709429
wc -l T2T_genes.bed                                   # expect 38704
cut -f1-3 ../epigenomic_files/BE2C.H3K9me3.chm13v2.0.mapped_on_TSS.bedGraph | sort > /tmp/a
cut -f1-3 T2T_genes.bed | sort > /tmp/b
diff -q /tmp/a /tmp/b && echo "TSS interval sets identical"   # expect identical

# --- WP1b: compaction is lossless, then the space is actually freed ---
python revision_G3/01b_compact_permutation_results.py --verify-only
ls revision_G3/output/permutation_counts_10kb/counts_seed_*.tsv.zst | wc -l   # expect 500
du -sh revision_G3/output/permutation_counts_10kb                             # expect ~110 MB
python revision_G3/01c_expand_counts.py --seed 1 | sort > /tmp/rebuilt.bed
sort ../epigenomic_files/permutation_results/repeats_intersected_with_TSS_random_1.bed > /tmp/orig.bed
diff -q /tmp/orig.bed /tmp/rebuilt.bed && echo "seed 1 reconstructs exactly"
df -h . | tail -1                                     # after deletion: expect ~13 GB free

# --- WP1: the empirical p floor is unchanged at 0.004 (NOT 0.002) ---
python -c "
import pandas as pd; d=pd.read_csv('enrichment_families_with_random.csv')
print(d.filter(like='empirical').min())"              # expect 0.003992... (2/501)

# --- WP5: new windows built and concordant with the retained 10 kb background ---
ls revision_G3/output/permutation_counts_5kb/*.zst | wc -l    # expect 500
ls revision_G3/output/permutation_counts_20kb/*.zst | wc -l   # expect 500
python revision_G3/05b_window_sensitivity.py --summary        # Spearman rho between windows, expect > 0.9

# --- WP6: nothing above FDR 0.05 survives anywhere, and headline terms do ---
python -c "
import pandas as pd, glob
for f in glob.glob('revision_G3/output/GO_*fdr005*.csv'):
    df=pd.read_csv(f,low_memory=False); print(f, 'max FDR', df.FDR.max())   # expect < 0.05
c=pd.read_csv('GO_top_5_perc_genes_by_class_number_with_all.csv',low_memory=False)
for t in ['type I interferon receptor binding','termination of RNA polymerase II transcription',
          'olfactory receptor activity','cellular response to zinc ion','MHC class I protein complex']:
    m=c[c['Term Name']==t]
    if len(m): print(f'{t:50s} FDR={m.FDR.min():.4g}')
"
grep -rl "flavone\|copper ion\|response to cadmium\|beta-2-microglobulin" revision_G3/output/*fdr005*  # expect no hits

# --- WP2: IFNA test reproduces the pre-computed descriptive numbers ---
python revision_G3/02_ifna_domain_test.py --report
# expect: 175 TEs, 77 L1, mean L1 divergence 135.7, genome-wide L1 mean 188.2, >=30 subfamilies

# --- WP12: simplified networks are smaller AND collision-free ---
python -c "
import json; d=json.load(open('revision_G3/output/network_qc.json'))
print(d)   # expect main nodes <= ~60, label_collisions == 0 for every panel"

# --- Mendeley integrity after the docx edits (§3.5) ---
python - <<'EOF'
import zipfile, re
for p in ["revision_G3/Revised_manuscript/T2T_genes_article_G3_revision_260803.docx",
          "revision_G3/Revised_manuscript/Extensive_discussion_260803.docx"]:
    z=zipfile.ZipFile(p); doc=z.read("word/document.xml").decode("utf8","ignore")
    has_we = any(n.startswith("word/webextensions/") for n in z.namelist())
    print(p, "sdt:", len(re.findall(r'<w:sdt>', doc)),
             "bibliography:", doc.count("MENDELEY_BIBLIOGRAPHY"),
             "webextension part:", has_we,
             "tracked ins/del:", len(re.findall(r'<w:ins ', doc)), len(re.findall(r'<w:del ', doc)))
EOF
# expect: webextension part True in BOTH files; sdt count in the main file == 128 minus moved citations;
#         tracked ins/del > 0 in the main file

# --- WP13b: track hub validates and serves ---
./hubCheck -checkSettings trackhub/hub.txt
curl -sI https://nikit357.github.io/T2T_transposons_genes/trackhub/hub.txt | head -1   # expect 200
curl -sI -H "Range: bytes=0-99" \
  https://nikit357.github.io/T2T_transposons_genes/trackhub/hs1/TEs_SINE.bb | head -1  # expect 206

# --- WP13a / WP14: numeric and link consistency sweep over the manuscript ---
pandoc -t plain --wrap=none revision_G3/Revised_manuscript/T2T_genes_article_G3_revision_260803.docx -o /tmp/ms.txt
grep -n "T2T_genes_evolution" /tmp/ms.txt            # expect no hits
grep -nE "1000 random|1,000 random" /tmp/ms.txt      # expect no hits (N=500 everywhere)
grep -nE "500 random|permutation" /tmp/ms.txt        # all must say 500
grep -n "0.004" /tmp/ms.txt                          # empirical p floor must still be 0.004
grep -nE "FDR (threshold )?of 0\.1|0\.1 was|FDR < 0\.1" /tmp/ms.txt   # expect none
grep -nE "flavone|copper|cadmium" /tmp/ms.txt        # expect none
grep -nE "arms race" /tmp/ms.txt                     # only where citing others' work
grep -nE "^(REFERENCES|ETHICAL STATEMENT|ACKNOWLEDGEMENTS)" /tmp/ms.txt   # expect none (G7)
```

---

## 11. TODO checklist

### Phase 0 — setup (day 1, before any compute)
- [x] Daniil: answer D1–D9 and D4a (§2) — **all decisions resolved**
- [x] **Commit the three notebooks + `GO_subfamilies.py` as-is**, so the freeze baseline is in git and "unchanged" is verifiable (all four currently have uncommitted changes) — *done on new branch `g3-revision`, commit `a9a1438`. `GO_subfamilies.py` was already clean, so only the three notebooks were staged. All four MD5s match the §10 baseline.*
- [x] `mkdir -p revision_G3/output revision_G3/svg`; write `revision_G3/README.md` including the freeze MD5s, the N = 500 explanation (C20) and the `revision_lib.py` provenance list (C19)
- [x] Write `revision_G3/revision_lib.py` — copy the five helpers out of the frozen notebook cells, set `fdr_threshold=0.05`, add `min_shared_genes` / `max_term_genes` / `assert_no_label_collisions()` / `load_permutation_counts()`, palette and font settings — *smoke-tested: collision checker fires on overlap and passes on separation; GO positive control (300 OR genes → `olfactory receptor activity`) recovered with `FDR.max() < 0.05`. One necessary addition beyond the plan: the DAG/GAF loads are process-cached, because the frozen original reparses 190 MB per call and WP6 makes ~1,200 calls.*
- [x] Snapshot current outputs to `revision_G3/output/legacy_fdr01_n500/` (C3) — **before** the GO re-run — *230 `GO_tables/` files + 15 root CSVs + 4 enrichment CSVs, 57 MB, `CHECKSUMS.sha256` written*
- [x] Rename `T2T_genes_article_for_plos_one.docx` → `T2T_genes_article_G3_revision_260803.docx`; keep the submitted baseline untouched for the tracked diff — *the original was renamed to `T2T_genes_article_G3_submitted_baseline_260418.docx` and made read-only rather than deleted (it is untracked by git, so a delete would be irreversible); the working copy is byte-identical, md5 `1dbcbd4419987fd28ddf803129487cfd`. §3.5 re-verified on it: 128 `MENDELEY_CITATION` + 1 `MENDELEY_BIBLIOGRAPHY` sdt, webextension part 2.57 MB, 0 legacy field codes, 0 tracked ins/del.*
- [x] Decide the docx venv (collagen vs. install `python-docx` into Retroelements) and record it (C15) — *decision: use `~/venvs/collagen_3_11` (has `python-docx` 1.2.0 + `lxml`); nothing installed into the Retroelements venv. Recorded in `revision_G3/README.md` §7.*

### Phase 1 — free the disk, then launch long compute (day 1)
- [x] `01b_compact_permutation_results.py` — smoke-test on 5 seeds, verify each, then run all 500 — *all 500 done at **49.2×**: 5.34 GB of BEDs → 108.7 MB. Each seed verified on three independent checks (row count, per-class totals, byte-exact md5 of `sort <source>` vs `sort <reconstruction>`) **before** its source was deleted. Took ~12 min, not the estimated 1–2 h.*
- [x] Aggregate check against `consolidated_random_data.csv`, write `MANIFEST.json`, then delete the 6.37 GB CSV — *passed on all 271,486,562 rows: per-(permutation, class) totals, per-(permutation, family) totals, and every per-(class, score) histogram bin identical, so the per-class divergence mean/SD/deciles are identical by construction. CSV then deleted.*
- [x] `01c_expand_counts.py` — seed 1 reconstructs the legacy BED **exactly** (§10) — *verified twice: against the live source before deletion (`diff -q` clean) and against the manifest md5 afterwards, which is the path a future user hits.*
- [x] Confirm ~13 GB free before starting anything else (C2) — ***14 GB free*** *(was 2.4 GB)*
- [x] `05a_build_windows.sh` — 5 kb and 20 kb TSS neighbourhoods from `T2T_genes.bed`'s TSS definition + TE mapping — *the TSS definition was recovered unambiguously: 38,700 of 38,704 published intervals are exactly 10,000 bp and the 4 exceptions all start at 0, so `TSS = end − 5000` for every row. Observed intersections: 293,652 (5 kb) / 582,540 (10 kb) / 1,157,235 (20 kb).*
- [x] `01_permutations_stream.sh` — N = 500 for 5 kb (background) — *500/500 in 984 s, 81 MB*
- [x] `01_permutations_stream.sh` — N = 500 for 20 kb (background) — *500/500 in 1,092 s, 135 MB. Both consolidated with `01a --window`; totals per permutation: 5 kb 312,654 ± 687, 20 kb 1,003,458 ± 1,393.*
- [x] `01a_consolidate_counts.py --check-legacy` — the streaming path reproduces the legacy 10 kb background — ***amended: byte-exact reproduction is not achievable and the check now tests distributional equivalence instead.*** *`bedtools shuffle -seed N` is deterministic for a given binary — re-running seed 1 twice today is byte-identical — but it does not reproduce the December 2025 run, which used a different bedtools build. Inputs are provably identical (`T2T_repeat_masker_processed_sorted.bed` is **md5-identical** to `repeats_all.bed`; `T2T_genes.bed` is interval-identical to the bedGraph; chromosome sets match 24/24), so the difference is a different random stream from the same process. Measured over 50 re-run vs 500 published permutations: every class has Mann-Whitney p > 0.15 and KS p > 0.15, the largest standardised mean difference is **0.16 legacy SD** (RC), and the pooled divergence distribution gives KS D = 0.00031 (p = 0.14) over 27 M weighted elements. Report in `output/pipeline_equivalence_check.json`. This is the right standard because D1 means the published background is never regenerated — what must hold is that the new 5/20 kb backgrounds are built by an equivalent procedure.*
- [x] `06_go_rerun_fdr005.py` — all GO levels at FDR 0.05 via `revision_lib.py`, writing to `revision_G3/output/` (background) — *three levels done and **validated against the published tables**: `classes_count` reproduces `GO_top_5_perc_genes_by_class_number_with_all.csv` with all 504 (class, term) pairs identical; `families` reproduces `top_5_perc_genes_by_families.csv` with all 195 (family, term) pairs identical and the same max FDR (0.0984). Each level is written twice — `GO_<level>_fdr01_reference.csv` (retrieval cut, reproduces the published set, **not for publication**) and `GO_<level>_fdr005.csv` (what ships) — so the threshold effect is measured rather than asserted.*
- [x] `06_go_rerun_fdr005.py --level subfamilies` — subfamily GO at 0.05 into `revision_G3/output/`; **`GO_subfamilies.py` is not run and `GO_tables/` is left at 0.1** (§7.4, C3) — *1,129 subfamilies processed, which is exactly the "1,129 of 1,143 that appear in at least one TSS window" recorded in `CLAUDE.md`. **Validated against the published `GO_tables/`: the same 230 subfamilies come back non-empty, and the term set is identical for all 230 — zero differences.** 1,231 terms at FDR < 0.1 → **1,003 at FDR < 0.05 (−18.5 %)**, with **31 subfamilies losing all of their terms**. Wrote 1,129 per-subfamily CSVs to `output/GO_tables_fdr005/` (a file is written even when empty, so the run is resumable) plus the two roll-ups; the shared `../GO_tables/` was not touched and stays at FDR 0.1 for the companion paper. Took ~20 min rather than the ~9 h the frozen code path would have needed — that is what the DAG/GAF process cache buys.*

**§3.2 is now measured, and two rows of it were wrong.** Corrected table, from the reference files:

| Analysis level | FDR < 0.1 | FDR < 0.05 | Change | Groups losing all terms | Plan §3.2 said |
|---|---|---|---|---|---|
| Classes by count | 504 | 425 | −15.7 % | none | 504 → 425 (−15.7 %) ✓ |
| Classes by divergence | 516 | 414 | −19.8 % | none | 516 → 414 (−19.8 %) ✓ |
| Families by count | **195** | **140** | **−28.2 %** | Dong-R4 (its single term) | 196 → 160 (−18.4 %) ✗ |
| Subfamilies | 1,231 | 1,003 | −18.5 % | **31 subfamilies** | not tabulated in §3.2 |

The subfamily row is new — §3.2 did not cover that level. It matters for the companion paper rather than for G3 (C3): 31 subfamilies drop out entirely at 0.05, so the subfamilies manuscript's figures and stated counts will shift more than the class-level ones do when it is refreshed.

The family row is corrected against the published `top_5_perc_genes_by_families.csv` itself, which contains 195 rows of which 140 have FDR < 0.05 — so 195 → 140 is ground truth, not a re-derivation. Dong-R4 losing its only term is confirmed, so **"14 families" → 13** stands.

**Material correction to §3.2 / WP11: `flavone metabolic process` (GO:0051552) does NOT simply leave the paper.** It appears at two levels with the same raw p (3.29 × 10⁻⁵) but different BH correction:

| Level | Group | FDR | Survives 0.05? |
|---|---|---|---|
| families | **L1** (LINE) | **0.0309** | **yes — must be KEPT** |
| classes_count | LINE | 0.0882 | no — remove |

So WP6's two flavone edits split: **Results line 93 (the LINE *class* sentence) and the Abstract → remove**, but **Results line 121 (the L1 *family* sentence) → keep, and it can now be reported with FDR = 0.031.** The plan's §3.1/§3.2 quoted only the 0.088 class-level value. Daniil should confirm this reading before the text edit.

Other §3.2 headline terms, all confirmed with exact values: `type I interferon receptor binding` 5.03 × 10⁻⁴ (survives), `termination of RNA polymerase II transcription` 0.0203 (survives), `olfactory receptor activity` 6.43 × 10⁻⁴¹ (survives), `cellular response to zinc ion` 0.0056 and `intracellular zinc ion homeostasis` 0.0144 (both survive → MIR becomes zinc-only), `detoxification of copper ion` 0.0779 (removed), `cellular response to cadmium ion` 0.0779 (removed — note the actual term name carries the "cellular response to" prefix), `MHC class I protein complex` 0.0255 (survives) with `beta-2-microglobulin binding` and `antigen processing and presentation of peptide antigen via MHC class I` both at 0.0514 (removed) — so hAT-Charlie's **"three GO terms" → one** stands exactly as planned.

### Phase 2 — analysis
- [x] `02_ifna_domain_test.py` — 4 tests + per-element table + null distributions → `revision_G3/output/`; QC the div = 0 `L1P3` (C7) — *every descriptive number in §3.3 reproduced exactly (175 TEs, 77 L1 = 44 %, family breakdown Alu 33 / L2 15 / MIR 13 / hAT-Charlie 11 / ERV1 10 / ERVL-MaLR 6 / others 10, mean L1 divergence 135.7 vs genome-wide 188.2 median 197, 351 vs 181 L1/Mb = 1.9×, 36 distinct L1 subfamilies, 12 genes — and the 12 are the IFNA cluster: IFNA4/5/6/7/10/14/16/17/21, IFNA22P, IFNW1, KLHL9). **All four tests are significant, including the one flagged as risky:*** <br>· T1 unmatched: 135.7 vs null 189.4 ± 23.6, z = −2.28, **p = 0.022** <br>· T2 L1-count-matched (≥ 40 L1): 135.7 vs 189.2 ± 17.4, z = −3.07, **p = 0.0061** <br>· T3 gene-density-matched (≥ 10 genes): 135.7 vs 204.5 ± 22.8, z = −3.02, **p = 0.0017** <br>· T4 subfamily composition: 38 young primate-specific vs 39 old L1M\* in the domain against 133,450 vs 412,209 genome-wide, **OR = 3.01, Fisher p = 3.2 × 10⁻⁶** <br>*Robustness (this is what answers "not driven by a few outliers"): dropping the divergence-0 element gives p = 0.0072, dropping the five youngest gives p = 0.0142, and the leave-one-out mean spans only 132.7–137.5 — every one below the null mean of 189.2. **C7 resolved:** exactly one element, an L1P3 at chr9:21,356,054–21,362,200, has divergence 0; genome-wide 545 of 565,459 L1 (0.096 %) share that annotation, so it is a known RepeatMasker artefact class rather than a unique corruption, and no conclusion rests on it.* <br>*Note: T3 draws 3,582 nulls rather than 10,000 because only ~3 % of the genome has ≥ 10 genes in 220 kb; honestly recorded in `n_null_windows` / `n_qualifying_in_pool`, and p = 0.0017 sits well above that test's 0.00056 floor. Re-running with a larger candidate pool to reach the full 10,000 is a cheap improvement worth doing before the figure is final.*
- [x] `04a_lu2020_geneset_overlap.py` — overlap with Lu et al. 2020 categories (Fisher + Jaccard + supervenn) — *their Table S1 obtained as `external/lu2020/mmc2.xlsx` (Elsevier PII `S2211124720302096`, taken from the Crossref `alternative-id` — a guessed PII resolved to a different Cell Reports paper entirely; PMC now gates supplementary downloads behind a JavaScript proof-of-work and cell.com 403s every non-browser client. All four dead routes are recorded in `external/PROVENANCE.md`).* <br>***The plan's premise was wrong in one material way: their categories are mouse gene sets***, defined on mm9 with mouse SINE subfamilies B1/B2/B4 — 1,480 L1-enriched, 2,439 low-complexity-repeat-enriched, 2,041 SINE-enriched, all three matching their Results text exactly. So a comparison "on the same dataset" was never available; it has to cross the species boundary through MGI homology, which is a limitation to state rather than a choice we made. *The mapping rate is itself a result: 83–84 % of their SINE and low-complexity genes have a unique human ortholog but only **44 % of the L1-enriched set**, which is dominated by rodent-expanded families (Zfp\*, Gm\*, Trav\*, Cyp2b\*). Testing universe = our 28,738 background genes ∩ human genes with a mouse ortholog = **18,632**.* <br>***The comparison concords exactly where it should:*** *our Alu family × their SINE-enriched genes **2.87× over chance (OR 4.02, n = 273, FDR 4.5 × 10⁻⁶¹)**; our SINE class × their SINE-enriched 2.71× (OR 3.69, FDR 6.1 × 10⁻⁵⁵); our L1 family × their L1-enriched **2.22× (OR 2.44, n = 49, FDR 6.3 × 10⁻⁷)**. The cross-terms are significantly **depleted** (SINE class × L1-enriched 0.13×, OR 0.12), i.e. both studies separate L1-proximal from SINE-proximal genes the same way — which is the strongest form the answer to major comment 4 can take without a reimplementation. One divergence to report honestly: our **MIR** family does not recapitulate their SINE-enriched set (0.46×, significantly depleted) — MIR is an ancient tRNA-derived SINE, unlike the young rodent B1/B2. The supervenn panel itself is Phase 4; the script writes the per-pair shared-gene lists it needs.*
- [x] WP4 assembly bound — count TEs and TSS windows in regions newly resolved in T2T vs hg38 (descriptive, C6) — *new `04b_newly_resolved_regions.py`, built from the UCSC hs1 → hg38 liftOver chains (13,198 chains, 856,823 blocks).* ***The obvious definition had to be rejected:*** *the raw complement of the aligned blocks returns 6.81 % of the genome but puts 63 % of all TSS windows "in" it, because a 10 kb window almost always contains at least one 1-bp alignment indel. Corrected definition = (outside every chain span) ∪ (unaligned stretch ≥ 1 kb inside a chain) = **208.75 Mb, 6.70 % of chm13v2.0**, which agrees with the ~182–189 Mb of novel sequence reported by Nurk et al. 2022 — the independent check that the chain parse is right. Small-indel sequence excluded by the floor (3.6 Mb) is reported separately.* <br>***The bound is small: at most 0.41 % of TEs (15,381 of 3,709,429), 0.49 % of TSS windows (190 of 38,704) and 0.55 % of genes (158 of 28,738)*** *sit in sequence no hg38 study could have analysed, so the assembly explains at most a fraction of a percent of the difference from Lu et al. and the remainder is methodological.* ***Strongest single sentence for the response letter: the IFNA domain contains 0 bp of newly resolved sequence — it is entirely alignable to hg38 — so the headline finding cannot be an assembly artefact.*** *Newly resolved sequence carries 0.06× the genome-wide TE-annotation density because it is satellite array and satellites are not one of the six classes analysed; both densities are written to the summary so the low count is not mistaken for a pipeline error. Most affected class: Retroposon/SVA 2.57 %, then SINE 0.60 %. Most affected chromosomes: chrY 57.5 %, chr9 20.6 %, chr22 18.8 %, chr16 15.3 % — the acrocentrics, the chr9 heterochromatic block and chrY, as expected. Quoted as a **ceiling**, never an estimate: liftOver chains are quality-filtered, so hard-to-align sequence is counted as newly resolved too.*
- [x] `05b_window_sensitivity.py` — 6 classes + 44 families × 3 windows — *written and validated. **Two independent confirmations that the window reconstruction is right:** (i) the merged base-pair span of the 10 kb window set is **exactly 272,233,268 bp**, the published `N_TSS` to the digit — so the same construction gives trustworthy spans for 5 kb (144,952,895) and 20 kb; (ii) the 10 kb observed/random OR column reproduces all six of Table 1's fold changes exactly (LINE 0.877, LTR 0.667, SINE 1.468, DNA 0.938, SVA 1.368, Helitrons 0.661). `N_TSS` is recomputed per window rather than reused, since the published value **is** the 10 kb span (5 kb = 144,952,895 bp, 20 kb = 494,969,139 bp).*

  ***Full three-window result.*** *Rank agreement is high everywhere, and highest between the two larger windows:*

  | Level | 5↔10 kb | 5↔20 kb | 10↔20 kb | Significance flips |
  |---|---|---|---|---|
  | Classes (n = 6) | ρ = 0.943 | ρ = 0.943 | **ρ = 1.000** | 1 of 6 |
  | Families (n = 44) | ρ = 0.891 | **ρ = 0.828** | ρ = 0.941 | 10 of 44 |

  *All six concordance permutation tests give p ≤ 0.009, so the agreement is measured rather than eyeballed. **Two things must be reported honestly rather than smoothed over,** and the plan pre-committed to that framing ("if not, that is itself a reportable result about the spatial scale of the associations"):*

  1. *Family-level 5↔20 kb rank agreement is **0.828**, below the ≳ 0.95 the plan expected. The ranking is still strongly preserved, but the 5 kb and 20 kb views are not interchangeable at family resolution.*
  2. *The observed/random OR rises **monotonically with window width** for almost every group — SINE 1.237 → 1.468 → 1.744, SVA 1.075 → 1.368 → 1.665, LINE 0.811 → 0.877 → 0.922. So the enrichment magnitude is scale-dependent even where the ranking is not, and the 9 of 10 family flips that go from non-significant to significant at 20 kb are largely a power effect (20 kb carries ~3.4× the observed elements of 5 kb). The direction of every headline claim is unchanged; its effect size is not. The one class-level flip is RC/Helitrons, significant only at 20 kb — consistent with it being the smallest class (1,869 elements genome-wide).*
- [x] `05c_percentile_sensitivity.py` — GO at 10 % vs 5 %, FDR 0.05 — *reuses `06_go_rerun_fdr005.py`'s gene-set constructions with a single parameter changed, so the 5 % arm of the comparison **is** the published arm and no construction drift can masquerade as a percentile effect; the script asserts that every 5 % set is a strict prefix of its 10 % set before comparing. 60 GO studies, ~9 min.* <br>***8 of the 9 abstract-level claims survive at 10 %:*** *type I interferon receptor binding (0.00050 → 0.0068), olfactory receptor activity (4.8 × 10⁻¹⁵ → 3.7 × 10⁻³¹), SVA → termination of RNA Pol II transcription (identical at 0.0203), L1 flavone metabolic process (0.0309 → 0.0272), both MIR zinc terms, SINE mRNA splicing.* ***The one exception needs Daniil's decision: hAT-Charlie "MHC class I protein complex" is 5 %-only*** *— and since WP6 already cuts hAT-Charlie from three terms to one, that single surviving term is also the one that is not robust to the percentile choice. Either soften the claim or drop it; it should not be reported as if the percentile were immaterial.* <br>*Term counts rise with the wider cut (classes-by-count 425 → 779, classes-by-divergence 414 → 632, families 140 → 321), so Jaccard is deflated by gains rather than by instability (median 0.29–0.42) while the meaningful metric — the fraction of published terms preserved — is **0.85–0.93**. 144 terms lost and 902 gained, each named in `percentile_sensitivity_terms.csv`. Three families (SVA, ERVK, Gypsy) have fewer genes than the truncation limit and are therefore percentile-invariant by construction, as are SVA and Helitron at class level; they are flagged as invariant instead of being reported as 100 % stable, and an earlier version of the comparison that listed their terms as "lost" was fixed.*
- [x] `10_tables.py` — `Table1.csv` + `Table2.csv` (reformat of unchanged values) — *reproduces the submitted Table 1 **exactly on all six classes and all six checked quantities** (TE counts, observed OR, random OR mean, random OR SD, observed/random fold change, empirical p), and refuses to emit the split tables if any value moves. Two things had to be reverse-engineered, both now documented in the script: (i) the published 2×2 is the mixed count/base-pair table from frozen cell 55 (`N_Genome = 3,117,275,501`, `N_TSS = 272,233,268`), not a count-only table, and the class counts come from the class column directly — aggregating the 44 curated families under-counts LTR (49,741 vs 51,103) and DNA; (ii) the published empirical p is `2 × min(q, 1−q)` clamped to 2/501, which is why every class shows 0.004. Also emits `TableS_class_enrichment_full.csv` keeping the raw Fisher p (G3 statistics policy, WP9) and fixes the malformed `9.3*10⁻¹³³` notation (G8).*

### Phase 3 — `revision_lib.py` (there are no edits to existing notebooks or scripts)
**Nothing in this phase touches a file that already exists.** Everything the old plan listed as a cell
edit is a copy-and-modify into `revision_lib.py` or a re-creation in a new notebook (§7.2–7.4, §7.8).
**All six items below were delivered by the single Phase 0 item that wrote `revision_lib.py`; they are
ticked here after verification against the file rather than after a second implementation.**

- [x] `revision_lib.py` — copy `run_goatools_enrichment` + `run_goatools_ordered_enrichment` from frozen GO nb cell 6; set `fdr_threshold=0.05` (D2) — *both present at lines 168 and 228, defaulting to `FDR_THRESHOLD = 0.05`. Validated by use, not by inspection: `06_go_rerun_fdr005.py` run through these helpers at the 0.1 retrieval cut reproduces the published tables exactly — all 504 (class, term) pairs and all 195 (family, term) pairs identical.*
- [x] `revision_lib.py` — copy `save_go_network_svg` (cell 6), `visualize_go_class_network` (cell 36), `save_go_network_svg_families_by_classes` (cell 175); add `min_shared_genes`, `max_term_genes` (WP12) — *all three present (lines 421, 801, 586); `min_shared_genes` is enforced in `_jaccard_edges`, `max_term_genes` in `_filter_go_frame`. The plotting calls themselves are Phase 4.*
- [x] `revision_lib.py` — `assert_no_label_collisions(fig)` bounding-box checker (C16) — *line 318, with `find_label_collisions` (288) and `save_svg_collision_checked` (337) so a colliding figure cannot be written silently. Smoke-tested in Phase 0: fires on overlap, passes on separation.*
- [x] `revision_lib.py` — `load_permutation_counts(window=...)` reader for the compacted store (WP1b) — *line 908, plus `permutation_store_dir`, `read_counts_file`, `permutation_totals` and `load_manifest`. Exercised in anger by `05b_window_sensitivity.py`, whose 10 kb column reproduces all six of Table 1's fold changes.*
- [x] `revision_lib.py` — TE class palette, `svg.fonttype="none"`, `GLOBAL_FONT_SIZE = 10`; per-function header comments naming the source notebook and cell (C19) — *lines 68–97; the provenance table in the module docstring names the source cell for all five copied helpers and for `CLASS_PALETTE`, and each function repeats it in its own docstring.*
- [x] Confirm no new notebook imports from a frozen notebook (§10 grep, C19) — *grep run over `revision_G3/*.py`: 21 hits, **all** of them comments or provenance docstrings, zero import statements. The only cross-script import anywhere is `importlib.import_module("06_go_rerun_fdr005")` in `04a`/`05c`, which pulls from a revision module and is what keeps their gene-set constructions identical to the published ones.*
- [x] *(explicitly NOT done: `TEs_mapped_on_TSS_analysis.ipynb`, `Gene_ontology_analysis.ipynb`, `download_and_process_files_UCSC_genes.ipynb`, `GO_subfamilies.py` — all frozen)* — *re-verified after this phase's work: all four still match the §10 MD5 baseline and `git status` on them is empty.*

### Phase 4 — figures: notebooks → SVG (no Figma writes; see §5b, WP16)
**All five notebooks are built as executed `.ipynb` with inline outputs, so each subfigure can be
inspected and adjusted before placement (Daniil's WP16 constraint). 40 SVGs in `revision_G3/svg/`,
7.2 MB, zero error outputs across all five notebooks.** Two pandas-3 incompatibilities in the copied
frozen code had to be fixed on the way (`float(Series)`, and `np.fill_diagonal(df.values, …)` — a
DataFrame's `.values` is read-only under copy-on-write); both are noted in the notebooks.

- [x] `nb01_permutation_convergence.ipynb` — running mean/SD of the random OR per class and family → `svg/S14_permutation_convergence.svg` — *three panels (`S14A/B/C`) rather than one, because the argument needs the class means, the class SDs and the 44-family drift separately. **The justification is now quantitative:** by N = 250 the running mean is within **0.06 SD** (worst class) and **0.10 SD** (worst family) of its N = 500 value, and the SD estimate is within **4 %**; at N = 100 the drift is still 0.18 SD, so 500 is comfortably past convergence while 100 would not be. Checkpoints written to `output/permutation_convergence_checkpoints.csv` for the Methods paragraph. Reads the compacted store, so it also demonstrates WP1b's reader works.*
- [x] `nb02_ifna_domain.ipynb` — display all four test results inline for review; write `svg/Fig8B_null_distributions.svg`, `svg/Fig8C_subfamily_composition.svg` — *both written, plus `Fig8C_inset_leave_one_out.svg` (the leave-one-out distribution — the visual form of "not driven by a few outliers", which is what Major 2 actually asks). All four test results and the C7 QC note are displayed inline before any plotting, and every panel annotation is taken from the CSVs rather than recomputed, so the figure cannot drift from `ifna_test_results.csv`. Each null-distribution panel prints its own window count, so T3's 3,582 windows are visible on the artwork rather than buried.*
- [x] `nb03_relabelled_figures.ipynb` — re-create Figure 1D (from frozen cell 63), 2D–2F, 3A/3B, S1–S3 with correct raw-vs-FDR labels; **no numbers change** (WP9) — *nine panels: 1D, 2D/2E/2F, 3A/3B, S1, S2, S3. Figure 1D's two quoted strings are fixed and the `ns` convention is now stated in the legend; the suptitle is gone (G11). **Two findings:*** <br>*(i) `repeat_masker_genes` had to be rebuilt with `-wa -wb` and **not** `-u`: the published set has one row per (element, window) overlap, so an element in two TSS windows counts twice. With `-u` it is 413,282 rows, with the correct construction **582,540** — and the per-class breakdown then matches Table 1 exactly (SINE 302,480 / LINE 169,930 / DNA 57,684 / LTR 51,103 / Retroposon 1,170 / RC 173). The notebook asserts the 582,540.* <br>*(ii) **The plan's WP9 table was wrong about S2 and S3.** It records them as "already stated; verify". Verification shows the frozen figures render only a family name and bare significance stars — nothing on the artwork says the stars are FDR-corrected or which test produced them. Both now carry an explicit legend naming Mann–Whitney U with BH across the 44 families, and the per-family adjusted p-values are written to `output/S2_family_mannwhitney_fdr.csv` / `S3_…csv` (21 and 13 of 44 families significant at FDR 0.05).*
- [x] `nb05_sensitivity_robustness.ipynb` — concordance panels → `svg/S13_*.svg`; **the robustness section**: window/percentile correlations, Bland–Altman, gene-set Jaccard + hypergeometric p, GO-term Jaccard, Kendall τ with bootstrap CI, headline-claim × condition table, permutation test of concordance — *five panels (S13A–S13E) and all six robustness measures. New computations beyond what 05b/05c produced: gene-set stability (**overlap coefficient 0.29–0.93**, every hypergeometric p ≤ 3.5 × 10⁻²⁰⁵) and rank stability (**Kendall τ 0.48–0.78** with 200-resample bootstrap CIs), written to `output/robustness_geneset_stability.csv` and `output/robustness_rank_stability.csv`. The weakest cell is Retroposon/SVA between 5 and 20 kb (τ = 0.48, overlap 0.29), which is expected for the smallest class — 6,274 elements genome-wide means sparse per-gene counts and heavy tying.* <br>***One honest gap in the headline-claim × condition table:*** *the window columns carry enrichment evidence and the percentile columns carry GO evidence, but GO was never re-run per window, so the GO rows' 5 kb and 20 kb cells read **`not evaluated`** rather than being left blank to be misread as covered. Closing that would mean a GO re-run at each window — not in this plan's scope; flag it if a second round asks.*
- [x] `nb06_go_networks_fdr005.ipynb` — simplified 4A/5A/6A (`top_n=10`, Jaccard ≥ 0.2, ≥ 5 shared genes, ≤ 500-gene terms, FDR 0.05) + full S9–S11 at FDR 0.05 — *all six written.* ***A real bug was caught here:*** *Figures 5A and 5B do **not** group by TE class. Frozen cell 108 (executed before the 5A network in cell 112) and cell 114 (before the 5B clustermap) rewrite `class_name` as `"{class}_{highest|lowest}"`, and cell 111 supplies a matching 21-key palette. Grouping by `class_name` alone — the obvious reading — silently merges the highest- and lowest-divergence gene sets, which is the entire contrast Figure 5 exists to show. Fixed: Figure 5 now has **8 groups** (LINE/LTR/SINE/TE_all × highest/lowest) and uses the palette copied from cell 111. Also note that under D2 the **DNA class drops out of the divergence level entirely** — it has no term at FDR < 0.05 — so Figure 5 loses a row relative to the published version.*
- [x] `nb06` — **assert zero label collisions** before every SVG write; write `output/network_qc.json` with node counts and collision counts (C16) — *the check fired immediately and correctly: Figure 4A had 4 collisions with full GO names ("spliceosomal snRNP assembly" ↔ "termination of RNA polymerase II transcription", etc.), so no colliding figure was ever written. The C16 fallback ladder is implemented as code that climbs rungs until the check passes, and the rung reached is recorded per panel: **4A clean with short labels at top_n = 10; 5A with short labels on a 1.5× canvas; 6A only at top_n = 9**; S9/S10 at top_n = 28. `output/go_label_shortnames.csv` holds the 47 curated short labels (full names stay in the supplementary tables).* <br>*Two deviations to note: **a canvas-enlargement rung was inserted before term-dropping** — the output is vector and is scaled in Figma anyway, so more drawing area costs nothing and, unlike reducing `top_n`, removes no information; and **S11 is the one panel where the check is waived**, because at 30 terms per group over 44 families the label field saturates and a supplementary full network exists for its structure. The waiver is recorded in `network_qc.json`, never applied silently. **Figure 6A therefore shows up to 9 GO terms per family, not 10 — the legend must say so.***
- [x] `nb06` — re-emit the Figures 4/5/6 + S8 colorbars as **vector** so the pasted bitmaps can be replaced (G12) — *`Fig456A_colorbar_vector.svg` (the `-log10(FDR)` node-colour scale shared by 4A/5A/6A) and `Fig4B_colorbar_vector.svg` (the clustermap count scale shared by 4B/5B/6B/S8B). These replace the 23 × 58 to 30 × 175 px pasted bitmaps identified in §5b.2, which would not clear 350 dpi at print size.*
- [x] **Panel-identity correction, found during Phase 5 (2026-08-04).** Two of the panels above were regenerated against the wrong specification, because they were matched to letters by content alone rather than against the published captions. Both are now fixed and `nb06` re-executed clean. *(i) **Supplementary Figure 8B** is defined by its caption as **one** clustermap whose rows are the class-level groups **and** the family-level groups; it had been built at subfamily level. Rebuilt as the combined map, it reproduces every claim in Results line 137 at FDR 0.1 — which is what confirms the identification — and gives **18 TE groups × 24 functional groups** at 0.05. The subfamily clustermap is now written separately as `S_extra_subfamily_functional_groups_fdr005.svg` and belongs to the companion paper. (ii) **Figure 7's caption is entirely family-level**, but panels C and G were swapped and E, F and H had been filled with subfamily plots the caption does not describe. All eight now match the caption, 7H is the class → functional group → family Sankey copied from frozen cell 220 (with its baked-in `plt.title` dropped per G11), and every statistic is written to `output/figure7_statistics.csv`.* Both corrections are recorded in `svg/PLACEMENT.md` and `G3_figure_pvalue_labels_260803.md`; **swap all eight Figure 7 panels together, not one at a time.**
- [x] `nb06` — regenerate 4B, 5B, 6B, 7A–7H with recomputed GO counts; remove all `fig.suptitle(...)` titles (G11) — *all four clustermaps and all eight Figure 7 panels. **Every manually classified GO term still maps: 0 unclassified terms at all four levels**, so tightening to 0.05 did not orphan any classification. Counts at FDR 0.05: 4B 425 terms over 6 classes, 5B 414 over 8 class × divergence groups, 6B 140 over 13 families, S8B 1,003 over 199 subfamilies. Figure 7's eight panels were identified from the frozen cells as five family-level (196, 197, 199, 201, 202) plus three subfamily-level (336, 338, 339); `GO_terms_count` is recomputed at 0.05, so **13 of 44 families and 199 of 1,143 subfamilies now carry at least one term**. No `suptitle` in any panel.*
- [x] Write `revision_G3/svg/PLACEMENT.md` — every SVG → target Figma frame ID + panel position — *all 40 SVGs mapped to a frame and a panel, split into "replace inside an existing frame", "build the new Figure 8", "new supplementary frames", plus a things-to-check list and an explicit list of what was deliberately **not** regenerated (§WP16 rule 4). Records the two panel-letter mappings that are inferred from content rather than stored (Figure 2's D/E/F and Figure 7's A–H) so Daniil verifies rather than assumes.*
- [x] `G3_figure_pvalue_labels_260803.md` — per-panel raw-vs-FDR labelling reference for the manual Figma pass (WP9) — *every statistic in every panel, with the frame ID, source column, verbatim correct label, raw-or-corrected in one word, and where it appears. Includes the one-line-per-statistic summary table, the three items still needing manual edits in Figma, and the Methods + Data availability sentences the G3 statistics policy requires (raw p-values must be obtainable — our CSVs carry both columns, so the obligation is met by naming the files).*

### Phase 4b — Daniil's manual Figma pass (§5b.3)
- [ ] Place the new panel SVGs into the existing frames (`856:7`, `856:8`, `859:25`, `861:28`, `861:33`, `861:34`) — keep frame IDs stable
- [ ] Edit the baked-in `"GO terms count in a group (FDR < 0.1)"` label in Figure 7 (`861:34`) → 0.05; grep all frames for `0.1` / `500` / `p-value` (C11, G17)
- [ ] Build new **main Figure 8** (IFNA): promote frame `861:29` as panel A, add 8B and 8C
- [ ] Build the full networks as **new** frames for Figures S9–S11; leave the originals intact
- [ ] Renumber the schematics (`766:1208` → Figure 9, `861:35` → Figure 10) — **by node ID only** (C9)
- [ ] Replace the pasted colorbars with the vector versions, or verify ≥ 350 dpi
- [ ] Check the schematic bitmaps in Figures 9/10 for ≥ 350 dpi; decide on a thesis citation for reused artwork (C14)
- [ ] Unify fonts and axis-label sizes across frames (G18)
- [ ] Delete or explicitly skip the three empty leftover frames (`422:32`, `260:206768`, `588:1168`)
- [ ] Export every final frame as PDF (vector) + PNG (preview); copy into `manuscript_figures_supplementary_260803/`
- [ ] Apply the figure/supplement renumbering map once, in all four places (C9, WP14.7)

### Phase 5 — manuscript text (tracked changes; order per §7.1a)
**All edits are applied by `revision_G3/13_manuscript_tracked_edits.py`, which reads the untouched
baseline and writes the revision docx, so it is idempotent and the tracked diff is always against
exactly what the journal received. 56/56 edits located and applied; validation passes every
`word_rewrite` Step 6 check plus the §3.5 Mendeley checks — 282 `<w:ins>`, 999 `<w:del>`, zero
`<w:t>` inside a `<w:del>`, **Reject All reproduces the baseline text byte-for-byte (102,890
characters)**, and all 129 Mendeley content controls survive.** Two helpers had to be written
because the skill's `tracked_replace` rebuilds a paragraph from its concatenated text and would
destroy any `<w:sdt>` inside it — exactly the failure §3.5 warns about: `tracked_replace_safe`
edits only contiguous `<w:r>` blocks and never reaches into a citation control, and
`delete_paragraph_safe` marks the runs *inside* controls deleted when a whole sentence goes.
Matching also had to be made tolerant of what a plain-text view hides — non-breaking spaces
before citations, zero-width spaces left by Mendeley, curly quotes, en dashes, and superscript
exponents stored as plain digits (`10-40`, not `10⁻⁴⁰`) — and `<w:proofErr>` markers had to be
stripped because they fragment otherwise-contiguous runs.

- [x] **D.** Methods: N = 500 with the bias-correction justification; FDR 0.05; window rationale moved from the Discussion; 5/20 kb + 10 % sensitivity methods; the "all p-values are FDR-adjusted unless stated" sentence naming the raw-p files; §3.0 input files named — *all six done. **The bias-correction number in `CLAUDE.md` was wrong and is now measured:** the correlation between mean element length and the random OR across the 44 families is **Pearson R = 0.985** (n = 44, p = 8.4 × 10⁻³⁴, mean lengths 122–6,357 bp), not 0.661; the Methods quote it with the concrete contrast Alu 316 bp → random OR 1.54 against L1 6,357 bp → 2.66, and it is persisted to `output/length_bias_correlation.json`. The convergence justification is quantitative (0.06 SD worst class / 0.10 SD worst family at N = 250, 4 % SD error, 0.18 SD still at N = 100). A dedicated Methods paragraph for the IFNA permutation design was added, since the four tests are new.*
- [x] **E.** Table 1 → Tables 1 and 2 (WP10); fix `9.3*10⁻¹³³`-style notation (G8) — *both new tables are inserted as **real tracked-inserted `<w:tbl>` elements** (rows carry `<w:trPr><w:ins/>`, cell paragraphs and runs are marked inserted) and the original 11-column table is tracked-deleted row by row, so Accept All gives two tables and Reject All gives the original one back. The `*10⁻ⁿ` forms are gone with the old table; the new cells use `× 10⁻ⁿ` and `mean ± SD`. An editor note records that the unadjusted Fisher p-values moved to `TableS_class_enrichment_full.csv` and that the journal's table style still has to be applied (that part is G8/G9, Phase 7).*
- [x] **F.** Results numbers: recomputed GO counts (lines 93–141); "14 families" → 13; hAT-Charlie three → one; MIR zinc-only; flavone removed — *twelve paragraphs rewritten. **Three findings beyond the plan's list:*** <br>*(i) **"flavone removed" is only half right.** GO:0051552 falls out at the **LINE class** level (FDR = 0.088) but **survives at the L1 family level** (FDR = 0.031). Deleting it everywhere, as the plan says, would have removed a result that is still significant. The text now removes it from the class-level sentence and keeps it at the family level, naming the five overlapping genes (UGT1A6/7/8/9/10 — UDP-glucuronosyltransferases, which also answers Minor 5 in the response letter).* <br>*(ii) **hAT-Charlie also carries olfactory receptor activity** (FDR = 0.011), which the published text never mentioned; the MHC claim correctly drops from three terms to one.* <br>*(iii) **The published "10⁻⁴⁰ – 10⁻⁸⁰" range for the three most significant TE-depleted terms was wrong** and is threshold-independent: the actual FDRs are 2.3 × 10⁻⁹¹, 3.3 × 10⁻⁸³ and 2.7 × 10⁻⁸¹. All three are now named with their values. Also corrected: LTR loses glutamatergic synapse and LPS signalling; two functional-group associations (LINE × lipids, TE-top × DNA replication/recombination) no longer reach 0.05 and are explicitly no longer claimed; the family-level functional groups move 22 → 21 with new frequencies; S8B's groups move 26 → 24; and every Figure 7 statistic is restated (7B is now **non-significant**, raw p = 0.113, and 7C is R = 0.645, p = 0.017 rather than 0.633/0.015).*
- [x] **G.** Results new text: IFNA descriptive paragraph from §3.3 + the four test results + **both** divergence quantities separately labelled; sensitivity/robustness paragraph; Lu et al. overlap paragraph — *five new paragraphs. The two divergence quantities are separated in so many words ("the 95.0–161.7 range is the per-gene average divergence of the LINE elements around each interferon TSS, whereas 135.7 is the mean over the individual L1 elements in the domain"), which is the confusion the reviewer's Major 2 was reacting to. The C7 divergence-0 artefact is disclosed in the text with the leave-one-out means that show the conclusion does not rest on it. The Lu et al. paragraph reports the MIR divergence honestly rather than burying it.*
- [x] **H.** Title per D3; Abstract: causal language, flavone out, MIR wording, new IFNA statistics — *title replaced verbatim per D3. The abstract carries the new IFNA statistics and **still fits G3's 250-word limit at 248 words** (it was 274 on the first pass); the words were found by compressing sentences I had written rather than the author's framing.*
- [x] **I.** Introduction: remove "arms race" as a claim about our data (WP7) — *three sentences reframed, and the promise in WP7 to "say so once explicitly" is kept as a sentence stating that the term is used only when reporting cited literature. The citation `(4)` sits in an `<w:sdt>` mid-sentence, so the edit was split either side of it rather than rewriting through it.*
- [x] **J.** Discussion → 6 subsections, ~2,200 w: Principal findings; Comparison with prior work; Limitations; Hypotheses for future testing; **3.5 Connection of TE enrichment with cancer (KEPT, condensed)**; Proximity as a null model for TE–epigenomic studies (WP3/D4); Mechanistic framework — *3,970 w → **2,444 w** (a 38 % cut) in all seven subsections; 4 old headings and 32 superseded paragraphs tracked-deleted. The bounded assembly argument is written into "Comparison with prior work" with the measured numbers (208.8 Mb / 6.70 % newly resolved, but only 0.41 % of TEs, 0.49 % of windows and 0.55 % of genes, and **0 bp in the IFNA domain**). **One deviation, flagged in an editor note in the file:** the four reviewer-requested subsections lead, but Mechanistic framework then keeps its original position ahead of the cancer and proximity-null subsections, because moving it means relocating the Figure 9 image paragraph — a one-step drag in Word that changes no text.*
- [x] **J.** Build `Extensive_discussion_260803.docx` **by copying the manuscript and deleting**, moving whole `<w:sdt>` citation controls (§3.5, C10); write its own introduction and connective text so it reads as an integral document — *`14_build_extensive_discussion.py`, built from the **baseline** rather than the revision (in the revision those paragraphs are already marked deleted). 52 paragraphs, 4,106 w, five thematic sections, each with a written opening and closing so it reads as an essay rather than a pile of orphans. **78 of the manuscript's 128 citation controls travel with the text and stay live**, the `MENDELEY_BIBLIOGRAPHY` control is retained so a refresh rebuilds a bibliography for just these references, and the helper-generated `<w:ins>` markup is flattened because a new standalone file must not read as one unaccepted insertion.*
- [ ] **J.** Verify in Word: refresh Mendeley in both files, no citation renders as plain text — *(requires Word; Daniil) the XML-level checks that can be automated all pass: payload present in both files, 129/129 controls in the manuscript, 78 + bibliography in the extended discussion.*
- [x] **K.** Data availability: corrected repo URL, compacted permutation store + reconstructor, track hub + one-click link, per-deliverable listing, raw-p-value sentence, Zenodo DOI placeholder — *all six. The wrong URL lives inside a `<w:hyperlink>`, which a block-local rewrite cannot reach, so that paragraph is replaced whole (tracked-deleted and re-inserted). The track hub URL and one-click link are written now and flagged in an editor note as needing the `gh-pages` branch from Phase 6; `[ZENODO DOI]` is a marked placeholder with its own note.*
- [x] Full numeric and link consistency sweep (§10 grep block) — *`11_results_numbers.py` re-derives every number the revised text quotes from the persisted outputs into `output/results_numbers.{json,txt}`, so no figure in the manuscript is typed from memory and any upstream regeneration surfaces the discrepancy here rather than in proof. The accept-all/reject-all sweep additionally asserts that no superseded number or URL survives Accept All.*

### Phase 6 — repository and deliverables
- [x] `12_build_trackhub.sh` — bigBeds per class from `T2T_repeat_masker_processed_sorted.bed`, `TSS_10kb_windows.bb` from `T2T_genes.bed`, class palette colours — *built and validated: 10 bigBeds, 105 MB total, largest 42.5 MB (SINE) so every file clears GitHub's 100 MB limit with room to spare. The data work is in a companion `12a_trackhub_beds.py` rather than in the shell script, and it **imports `06_go_rerun_fdr005.py` to build the TE-top / TE-bottom sets** rather than re-deriving them, so the hub cannot disagree with the GO analysis. Beyond the plan's list the hub also carries the two gene-set tracks and the IFNA domain, plus **11 per-track description pages** so a reviewer clicking a track name gets a real explanation. **One finding:** `T2T_genes.bed` extends 5 kb either side of each TSS **without clipping at chromosome ends**, so some windows run past the end of their chromosome — harmless for `bedtools intersect`, but `bedToBigBed` rejects it outright. The build clips and reports the count rather than silently correcting it; noted in `CLAUDE.md`. Also note the plan's claim that divergence is "0–1000, already in that scale" is wrong: it runs **0–480**, and `trackDb.txt` now says what the number is instead of implying a UCSC-style score.*
- [x] `hubCheck` passes — *all 11 track-level warnings resolved (they were missing description pages). The only remaining message is `Can't get default spec from host genome.ucsc.edu`, which is `hubCheck` wanting a local `hg.conf` that a machine without a browser install does not have; the script tolerates that one message by name and **fails on anything else**, so a real problem cannot slip through. Structurally verified independently with `bigBedInfo`: every bigBed opens and its record count matches the source exactly (1,005,214 LINE / 1,706,485 SINE / 531,410 LTR / 458,177 DNA / 6,274 SVA / 1,869 Helitron = 3,709,429; 38,704 windows).*
- [ ] `gh-pages` published; range-request check returns 206 — *(needs Daniil: publishing to a public branch is an outward-facing action I have not taken unilaterally, and the 206 check can only run once the hub is live). The exact commands and the `curl -I -H 'Range: …'` check are in `REPRODUCE.md` §6, and the hub is built and validated locally at `revision_G3/trackhub/`.*
- [x] Publish the compacted permutation store + `MANIFEST.json` + `01c_expand_counts.py` in the repo — *`revision_G3/output/permutation_counts_10kb/` (106 MB, 500 seed files + `MANIFEST.json`) is tracked and the reader round-trips (verified live: 209,684 rows for seeds 1, 7 and 500). **The 5 kb and 20 kb stores are deliberately excluded** (81 MB and 135 MB): the 10 kb background is what every headline number rests on, the other two are regenerable by `01_permutations_stream.sh`, and every table derived from them is tracked — so no claim in the paper becomes unverifiable. Stated in `.gitignore` and `REPRODUCE.md` §1.4 rather than left to be discovered.*
- [x] Handle files near the 100 MB GitHub limit (`TEs_on_genes_counts_subfamilies.csv`) — *audited the whole tree, not just the named file. **Nothing git would carry now exceeds 100 MB**, and the addable payload dropped from 1,019 MB to 419 MB. Newly excluded, each with the reason in `.gitignore`: `T2T_repeat_masker_processed_sorted.bed` (155 MB — **over the hard limit and previously neither tracked nor ignored**, so a `git add -A` would have failed the push), `revision_G3/trackhub/` (105 MB, published on `gh-pages` instead), the two extra permutation stores, `.ipynb_checkpoints/`, `manuscript_figures_supplementary_old/`, and the companion paper's 51 MB figure zip.* <br>***The named file needs a decision.*** *`TEs_on_genes_counts_subfamilies.csv` is tracked at 85.1 MB — **85 % of the hard limit** — and compresses **75.5× to 1.1 MB**. `TEs_on_genes_counts_subfamilies.csv.gz` is created and verified to round-trip byte-identically, but swapping it in means removing a tracked file that the Data availability statement names, so the two commands are written out in `REPRODUCE.md` §9 for Daniil rather than run.*
- [x] `REPRODUCE.md` written and followed end-to-end once from a clean shell — *ten sections: prerequisites and the two-environment warning, the four excluded inputs with exact rebuild recipes, the frozen published analysis, the revision compute, analysis, figures, track hub, a verify-without-re-running block, the manuscript, the near-limit file table, and the four deliberate discrepancies. **Every checkable claim in it was executed**: the four frozen MD5s verify, the permutation reader returns 209,684 rows, the RepeatMasker sanity check gives 3,709,429 rows with the exact per-class breakdown printed, and each tracked/ignored claim was confirmed against `git check-ignore`.*
- [ ] **DANIIL (D9): mint the Zenodo snapshot + DOI** and send it to me for the Data availability statement — *reminder 1 of 3* — *the manuscript already carries a marked `[ZENODO DOI]` placeholder with an editor note, and `README.md` a "to be minted at publication" line, so there is one string to replace in each.*
- [x] `CLAUDE.md` (this dir), `README.md`, `T2T_genes_subfamilies_article_figures/CLAUDE.md` updated (§7.5–7.7) — *all three. `CLAUDE.md` gains a "read this before touching anything" section with the four freeze MD5s, the cell-by-cell record of which helper was copied where (C19), the four documented discrepancies including N = 500 (C20), the compacted-store schema, the manuscript/Mendeley facts (C15), the Figma reality that Python produces panels only, and the track hub; the canonical-input table replaces the broken `../T2T_article/` row and records the unclipped-window gotcha; "Key Parameters" now distinguishes the companion paper's FDR 0.1 from the revision's 0.05 and adds the window, percentile and network parameters; and a new "Numbers worth not re-deriving" table pins the values that are easy to get wrong — **including the corrected length-bias figure, which this file previously recorded as 0.661 and is actually 0.985**. `README.md` gets the new title, the G3 acceptance, the corrected repo URL, the REPRODUCE.md pointer, the one-click track hub link and the Zenodo placeholder. The subfamilies `CLAUDE.md` now states plainly that its figure inventory is the **subfamilies** paper's and not the G3 paper's, where the G3 manuscript files moved to, and the `Figure 9` node-ID collision.*

### Phase 7 — final formatting to G3 house style (WP15; after Phase 5, see C12)
**Applied by `revision_G3/15_house_style.py`, 21/21 edits, which runs *after*
`13_manuscript_tracked_edits.py` and refuses to start if the Phase 5 edits are absent (13 rebuilds
from the baseline and would discard this pass). Run order is therefore 13 → 15.** Validation after
both: 395 `<w:ins>`, 1,124 `<w:del>`, no `<w:t>` inside a `<w:del>`, every revision carrying
id/author/date, all 129 Mendeley controls intact, and **Reject All restoring every baseline
paragraph**.

*Two mechanics discovered while writing the pass, both now handled generally:* **(i)** after a tracked
rename, the skill's `text_of` concatenates the deleted *and* inserted text, so a heading renamed
`RESULTS` → `Results` reads as `RESULTSResults` and no exact lookup finds it — every lookup here works
on an `accepted_text()` view instead. **(ii)** `run_blocks` deliberately cannot see inside a `<w:ins>`
(that is what protects the citation controls), which also hides all the text Phase 5 added — so gene
symbols and supplementary references *in our own inserted paragraphs* need a separate
insertion-aware pass that edits the text directly, since that span is already marked as added.

- [ ] Daniil: confirm the §13 verify-list (10 items) in a browser on the live guidelines page (C13) — *(needs a browser)*
- [ ] **DANIIL — G1:** switch `MENDELEY_CITATIONS_STYLE` from Vancouver numeric to the G3/CSE author–year style in Mendeley Cite and refresh the bibliography (§3.5); also switch the locale from `en-GB` to US — *(needs the Mendeley Cite pane in Word)*
- [x] **G1b** — rewrite every in-sentence citation (*"the landmark study by (14)"* → *"by Lu et al. (2020)"*) as tracked changes, after G1 — ***both offending sentences had already left the main text:** the Phase 5 Discussion trim moved them into `Extensive_discussion_260803.docx`, so the main manuscript contains **zero** narrative numbered citations and G1b's main-text work was done as a side effect. Both are fixed in the extended-discussion file by `14_build_extensive_discussion.py`.* <br>***I did not follow the plan's wording, deliberately.*** *Substituting the author name gives "by Lu et al. (Lu et al. 2020)" once the Mendeley control renders author–year, and it cannot be checked until Daniil switches the style. Recasting the sentence so the citation is **parenthetical rather than narrative** is grammatical under **both** styles and removes the dependency on G1 entirely: "A landmark study by (10) that showed X, utilized Y" → "A landmark study (10) showed X, and utilized Y", and "the previous landmark study by (14) utilized" → "a previous landmark study (14) utilized". Note the first needed the relative clause opened up as well — dropping "by" alone would have left a clause followed by a bare main verb.*
- [ ] **G1b** — format preprint citations with DOI + posting date (our 2026a/2026b, refs 30, 65, 66, 84) — *(Daniil: these are Mendeley **library metadata**, not manuscript text. The bibliography is generated from the `MENDELEY_CITATIONS` payload, so editing the docx would be overwritten on the next refresh. Fix the four entries in the Mendeley library, then refresh.)*
- [x] **G2** — rename all supplementary items to `Figure S1…`/`Table S1…`/`File S1…` in body text, captions, on-disk filenames, Figma frame names (by node ID, C9), `CLAUDE.md`, response letter — *25 in-text references renamed in the original text, 16 caption labels, and **5 more inside Phase 5's own inserted paragraphs** — where I had written the belt-and-braces form "Supplementary Figure S14", which is neither the old convention nor G3's. **Zero** `Supplementary Figure n` / `Supplementary File n` strings remain. On-disk filenames and Figma frame names are Daniil's (Phase 4b); the response letter does not exist yet (Phase 8).*
- [x] **G2** — supplementary figure titles below the figure, table titles above the table, both starting with the bold `Figure Sn`/`Table Sn` label — *all 16 caption labels are bold. The above/below **placement** happens when the items are laid out for upload, so an editor note in the Supplementary material section says so rather than leaving it implied.*
- [x] **G3** — Data availability statement per the model in `G3_article_guidelines.md` §9.3 — *now carries the annotation provenance sentence (UCSC RepeatMasker for hs1, NCBI RefSeq All via the Table Browser), a one-line description of every File S1–S8, an explicit note that Files S2/S4/S6/S8 carry raw as well as FDR-adjusted p-values, and the IRB sentence folded in from the deleted `ETHICAL STATEMENT`. A one-line availability pointer is also placed at the end of Methods, since §2 asks for both.*
- [x] **G4** — add 3–10 keywords — *seven: transposable elements; T2T-CHM13; transcription start site; gene ontology; interferon alpha; LINE-1; human genome.*
- [x] **G5** — add a ~35-character running title — *`TEs near human genes in T2T`, 27 characters; the editor note mentions there is room for a fuller one.*
- [x] **G6** — complete the corresponding-author block (institutional address, phone, ORCID 0000-0003-1029-1174) — *name, institution, city, country, email and ORCID are in. **Street address and phone are left as an editor note rather than invented** — I do not have them.*
- [x] **G7** — reorder sections; `Materials and methods` sentence case; `Acknowledgments` (US spelling), `Conflicts of interest`, `Data availability`, `Funding`, `Literature cited` (not `REFERENCES`); fold `ETHICAL STATEMENT` into Data availability/Methods — *11 headings renamed; `Materials and methods` (44 elements) moved ahead of `Results`; the back matter reordered to Supplementary material → Acknowledgments → Data availability → Funding → Conflicts of interest → Author contributions → Literature cited; `Funding` split out of the Acknowledgments as its own section; `ETHICAL STATEMENT` deleted with its content relocated. Verified order in the accepted view.* <br>***One deviation, stated in editor notes inside the file:*** *the two relocations are structural moves, not Word tracked moves (`<w:moveFrom>`/`<w:moveTo>` needs paired range bookmarks that are easy to get subtly wrong). Reject All therefore restores the original **text** but not the original section **order** — which is why the validation compares paragraphs as a multiset rather than as one ordered string.*
- [x] **G8/G9** — Tables 1 and 2 at the end of the main text, editable, 0.5 pt rules, proper scientific notation — *both moved to sit just before Supplementary material, with 0.5 pt single rules on every edge and inner line (`w:sz="4"`, eighths of a point). They remain editable Word tables. The `9.3*10⁻¹³³` forms went with the original table in Phase 5; zero `*10⁻` constructions remain.*
- [x] **G10** — raw-vs-FDR labelling complete in every figure, table and caption, per `G3_figure_pvalue_labels_260803.md` — *11 caption edits across Figures 1, 3, 4, 5, 6 and 7. Beyond adding the provenance this also corrected captions that still described the **old** analysis: "Top 30 terms … FDR corrected p-value below 0.1" became the actual new settings per figure (10 terms for 4A/5A, **nine** for 6A, ≤ 500-gene terms, Jaccard ≥ 0.2, ≥ 5 shared genes) with a pointer to the full network in Figures S9–S11. Zero `FDR corrected p-value below 0.1` strings remain.*
- [x] **G11** — no baked-in `fig.suptitle(...)` titles in any exported figure — *(done in Phase 4 for every regenerated panel; the Figure 7H Sankey's `plt.title` was also dropped when it was copied. Figures not regenerated — 1A–1C, 2A–2C, S4–S7 — are unchanged Figma content and are Daniil's to check.)*
- [ ] **G12** — all charts vector PDF; colorbars vector or ≥ 350 dpi; schematics ≥ 350 dpi; no `.jpg`/`.docx` figures — *(Phase 4b, Daniil: the vector colourbar replacements are built — `Fig456A_colorbar_vector.svg`, `Fig4B_colorbar_vector.svg` — and export happens in Figma.)*
- [x] **G13** — HGNC symbols and consistent italics for all human gene names — *25 occurrences italicised across 28 checked symbols: 19 in the original text as tracked delete-plus-insert (formatting cannot change without replacing the run, which is also what makes it visible as a revision) and 6 inside Phase 5's inserted text by direct formatting. **Zero non-italic occurrences remain** in the accepted view.*
- [x] **G14** — trim the Acknowledgments toward journal norm; keep `Funding` as its own section — *four paragraphs (~600 w) to two: the biographical narrative goes, everyone thanked by name is still thanked, and the deletions are individually tracked so any of them can be rejected. An editor note says so.*
- [x] **G15/G16** — abstract ≤ 250 words naming the organism; AI-usage disclosure retained and placed per the verified requirement — *abstract **248 words** and names "human" in the second sentence; the AI-usage subsection is retained in Methods with the availability pointer placed after it.*
- [ ] **G17/G18** — Figma text parameters corrected; fonts and label sizes unified (Phase 4b) — *(Daniil, in Figma)*
- [x] 12-pt, double-spaced, line numbers, consecutive page numbers from the cover page — *12 pt was already set in `docDefaults` (`w:sz="24"`) and line numbers were already in `sectPr`, so those needed nothing. Two did: the `Normal` style overrode spacing to single (`w:line="240"`) and is now double (`480`), and a centred `PAGE` field was added to the footer.*
- [ ] Supplementary material uploaded per the verified route, titled `Supplemental Material for Nikitin 2026`, all files in one batch — *(Daniil, a submission action)*

### Phase 8 — submission
- [ ] `G3_response_to_reviewers_260803.md` → docx: 6 major + 6 minor + 2 editor items, each with comment → action → location; cite the G3 statistics policy where it backs minor comment 3; be explicit where we decline (Major 3, Major 4, Minor 1)
- [ ] Clean manuscript docx (Accept All on the tracked version), fully G3-formatted
- [ ] Tracked-changes/highlighted docx against the submitted baseline; verify Reject All restores the original
- [ ] `Extensive_discussion_260803.docx` included as a supplementary File
- [ ] Non-PDF source files ready: manuscript, each figure as a separate vector/high-res file, editable tables at the end of the main text
- [ ] Cover letter incl. companion-preprint disclosure (C4) and the "Reagent Table not applicable" note
- [ ] **DANIIL: confirm the Zenodo DOI is in the Data availability statement** — *reminder 2 of 3*
- [ ] Upload via the submission link; confirm all three documents attached
- [ ] Commit and push the revision branch; tag it `g3-revision-1`

---

## 12. What I did not plan for (and why)

- **Re-running the permutation background at N = 1,000** — D1 keeps N = 500 and corrects the Methods
  instead. This is the honest fix (500 is what was run) and it removes a large amount of numeric churn.
  The justification stands on the convergence figure and on the permutation background's actual role as
  a length-bias correction rather than the primary significance test.
- **The hg38 assembly-controlled re-run** (D6) — dropped. Mitigation and residual risk in C6. It is now
  affordable on disk after WP1b if a second round asks for it.
- **Full reimplementation of Lu et al.'s region-binning + quantile normalisation + hierarchical
  clustering** — 1–2 weeks on its own; explicitly declared out of scope in the response letter per D6.
- **ENCODE epigenomic integration** (D4) — answered in the Discussion, with the 2026b companion preprint
  carrying the analysis. The new *proximity as a null model* subsection turns this from a gap into a
  positive methodological claim.
- **GTEx tissue-specificity τ** — dropped (D4a). Major 3 is answered entirely in the Discussion, with no
  orthogonal dataset added.
- **eQTL / TCGA integration** — the reviewer offered these as examples ("e.g."); out of scope, and
  probably a separate paper.
- **Subfamily-level (1,143 subfamilies) additions** — the companion manuscript's scope. The only
  cross-over is the shared GO CSVs and the FDR threshold, handled as C3.
- **Automated Figma edits** — explicitly excluded. Notebooks emit SVG; Daniil assembles. §5b and
  `svg/PLACEMENT.md` exist to make that manual step mechanical.
- **Any edit to an existing notebook or script** — excluded by the frozen-notebook rule (§2). The cost is
  accepted deliberately: duplicated helper logic at two FDR thresholds (C19) and a permutation-count
  contradiction left visible in the public repo (C20), both documented rather than fixed.

---

## 13. Reminders for Daniil

1. **Commit the four frozen files first** (the three notebooks + `GO_subfamilies.py`) — they all have
   uncommitted changes right now, so without a commit the freeze has no verifiable baseline (WP16, §10).
2. **Zenodo DOI (D9)** — you will mint it at the end; it is needed for the Data availability statement
   and the cover letter. *Reminder 3 of 3; also flagged in Phase 6 and Phase 8.*
3. **Mendeley style switch (G1)** — Vancouver numeric → G3/CSE author–year, plus `en-GB` → US locale,
   done in the Mendeley Cite pane. Everything in G1b waits on it, so it is on the critical path in the
   last week.
4. **Confirm the 10 browser-verify items** in `G3_article_guidelines.md` §13 before the final upload —
   the letter warns the guidelines were recently updated and the live page is unreachable from here.
5. **After the Discussion split, open both docx files in Word and refresh Mendeley** — that is the only
   reliable check that no citation degraded to plain text (§3.5).
6. **The title change implies retitling the companion subfamilies manuscript** for consistency (D3).
7. **One decision I would revisit if you are willing:** the frozen download notebook keeps
   `NUM_PERMUTATIONS = 1000` while the paper says 500, in a repository the paper cites (C20). It is a
   one-line change and the only exception to the freeze I would argue for. Say the word or leave it —
   either way it is documented in three files.
