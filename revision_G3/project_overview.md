# `revision_G3/` — project overview

Background and detail for the G3 revision package. `CLAUDE.md` beside this file is the operating
manual (rules, run order, API, the numbers not to re-derive); this file explains *why* the package is
shaped the way it is, what every work package answers, what every script found, and which of the
plan's assumptions turned out to be wrong.

Sources: `../G3_revision_implementation_plan_260803.md` (the plan, with a per-item completion log),
`../G3_reviewer_report_260802.md` (the decision letter), `../G3_article_guidelines.md` (house style +
gaps G1–G18), `../REPRODUCE.md` (clean-checkout run order).

---

## 1. The manuscript and its status

**Title (revised, per decision D3):** *"Telomere-to-telomere co-mapping of transposable elements and
human genes identifies a cluster of young L1 elements in the interferon-alpha domain"* — replacing
*"Evolutionary arms race between transposable elements and human genes…"*.

**Journal:** *G3: Genes|Genomes|Genetics*, manuscript **G3-2026-406828**, **conditionally accepted**
with minor revisions. Letter dated 2026-07-29, 30-day window → target submission on or before
**2026-08-28**. Preprint of the same work: DOI [10.32942/X2FM2M](https://doi.org/10.32942/X2FM2M).

Nothing in the reviewer report challenges the core result. The requests are (a) reframing causal
language, (b) statistical rigour around the headline IFNA–L1 finding, (c) orthogonal validation,
(d) sensitivity analyses, (e) a shorter Discussion, and (f) cosmetic/consistency fixes.

The companion **subfamilies** manuscript (in preparation, figures in
`../T2T_genes_subfamilies_article_figures/`) shares data files and the Figma canvas with this one.
That coupling is the source of several of the caveats below.

---

## 2. The four decisions that shape everything

1. **The permutation background stays at N = 500** (D1) and the Methods text is corrected *down*
   from its erroneous "1,000". 500 is what was actually run; Results, Discussion and Figure 1D
   already said 500; Table 1's empirical p = 0.004 = 2/501 confirms it independently. Consequence:
   **no enrichment CSV, table or figure is regenerated for permutation reasons** — which removes a
   whole class of numeric churn. Justification is *demonstrated* (a convergence figure) rather than
   asserted, and rests on the correct argument: the permutation background exists to **correct the
   odds ratio for element-length bias**, not to serve as the primary significance test (that is the
   FDR-adjusted Fisher exact p). What N must deliver is a stable mean and SD of the random OR, not a
   small p-value floor.
2. **GO FDR tightens to 0.05 everywhere** (D2), with no "suggestive 0.05–0.1 band" in main text,
   supplementary figures or supplementary tables. Every GO count in every figure, table and caption
   is recomputed.
3. **No new epigenomic or genome-wide orthogonal analysis** (D4, D4a). Major 3 is answered in the
   Discussion with a positive methodological argument (below), and the requested epigenomic
   integration already exists as the companion 2026b preprint. GTEx tissue-specificity τ was dropped
   too — a single τ analysis would have invited "why expression but not chromatin or eQTL?" without
   settling the comment.
4. **Figures are produced only up to SVG.** All revision plotting happens in notebooks so each
   subfigure can be inspected; the notebooks write SVGs to `svg/`, and Daniil places them into Figma
   frames by hand. No script writes to Figma.

Also decided: neutral title (D3); full sensitivity at 5/10/20 kb × 5/10 % with a dedicated robustness
comparison (D5); Lu et al. comparison by gene-set overlap only, with the hg38 re-run and the
clustering reimplementation declared out of scope (D6); Discussion cut to ~2,200 words but **keeping**
subsection *3.5 Connection of TE enrichment with cancer*, with the excised material rebuilt as a
standalone "Extensive discussion" docx (D7); GitHub Pages track hub (D8); Zenodo DOI minted by Daniil
at the end (D9).

### The frozen-notebook rule

Daniil's constraint was explicit: new notebooks only, in a new subdirectory, covering only what the
reviewer asked for; do not touch the existing ones. Three structural consequences hold throughout:

1. **New notebooks only, minimal scope.** Five notebooks produce every figure this revision needs and
   nothing else. A figure that does not change is not regenerated (Figures 1A–1C, 2A–2C, S4–S7 keep
   their current SVGs).
2. **Helpers are copied, not edited in place** → `revision_lib.py`, where the FDR threshold becomes
   0.05 and the new arguments are added. The originals keep working unchanged at 0.1.
3. **Nothing regenerates data the frozen notebooks own.** `enrichment_*.csv` and the class/family
   enrichment tables stand as published; the new notebooks *read* them.

The accepted cost: duplicated helper logic at two FDR thresholds (caveat C19) and a
permutation-count contradiction left visible in the public repository (C20). Both are documented
rather than fixed, in three files each.

---

## 3. Reviewer comment → work package map

| Reviewer item | WP | Type | Delivered as |
|---|---|---|---|
| Minor 1 — 500 vs 1,000 permutations | WP1 | text + justification figure | `nb01` → S14A–C; Methods correction |
| (D1 companion) 11 GB permutation store | WP1b | new script | `01b`, `01c`, `output/permutation_counts_10kb/` |
| Major 2 — IFNA statistical test + counts/subfamilies | WP2 | new analysis + new figure | `02_ifna_domain_test.py`, `nb02` → new main Figure 8 |
| Major 3 — orthogonal data (ENCODE/GTEx) | WP3 | **text only** (D4) | new Discussion subsection *Proximity as a null model* |
| Major 4 — direct comparison with Lu et al. | WP4 | gene-set overlap only (D6) | `04a`, `04b` |
| Major 5 — window + percentile sensitivity | WP5 | new analysis + tables + notebook | `05a`–`05c`, `nb05` → S13A–E |
| Minor 2 — FDR 0.1 → 0.05 | WP6 | re-run GO + text | `06_go_rerun_fdr005.py`, `nb06` |
| Major 1 — reframe "arms race"/causality | WP7 | text only | Title, Abstract, Intro, Results, Discussion |
| Major 6 — shorten Discussion + extras docx | WP8 | text + new deliverable | `13_…`, `14_build_extensive_discussion.py` |
| Minor 3 — p-value provenance in every label | WP9 | figure/caption edits + new doc | `nb03`, `../G3_figure_pvalue_labels_260803.md` |
| Minor 4 — Table 1 too wide | WP10 | table restructure | `10_tables.py` → Table1/Table2 |
| Minor 5 — "flavone metabolism" | WP11 | text only | GO:0051552 confirmed; see §6.4 |
| Minor 6 — simplify Figures 4A/5A/6A | WP12 | re-plot + supplementary | `nb06` → simplified + S9–S11 |
| Editor — data/code links | WP13a | repo + text | corrected URL, `REPRODUCE.md`, published store |
| Editor — browser/track instance | WP13b | new deliverable | `12_build_trackhub.sh` |
| Response letter, tracked changes, package | WP14 | manuscript production | Phase 8, open |
| Final formatting to G3 house style | WP15 | manuscript production | `15_house_style.py` (G1, G12, G17/G18 are Daniil's) |
| Notebook figure surface + `revision_lib.py` | WP16 | infrastructure | this directory |

≈ 22 person-days plus ~6–10 h of compute. D1 and D6 removed the two big re-runs, so compute was never
on the critical path.

### The Major 3 answer, in full

The new Discussion subsection makes four points in order: (1) the design is correlative and static —
it measures where elements are, not what they do; (2) **the positive methodological argument** — a
large literature reports that some TE family carries some epigenomic mark near genes of some
function and reads that as regulatory recruitment, but any such claim needs a **positional
baseline**, because a share of the signal follows from nothing more than certain TE families sitting
physically closer to certain gene classes than chance allows. This paper supplies exactly that
baseline, genome-wide, at T2T resolution, length-bias-corrected. Framed this way the absence of
epigenomic data here is the point rather than a gap; (3) epigenomics is the declared next step and it
exists — the 2026b companion preprint (DOI 10.64898/2026.03.19.712972; 7 marks × 12 cell lines ×
1,122 subfamilies on T2T); (4) named limitations and concrete follow-ups (CRISPRi of IFNA-domain L1s
under IFN stimulation, SVA_B deletion at the SSU72L cluster, conditioning published mark-based
enrichments on this proximity baseline). The companion preprint is disclosed in the cover letter so
the editor reads it as a declared companion rather than undisclosed overlap (C4).

---

## 4. The compacted permutation store (WP1b)

`../../epigenomic_files/permutation_results/` held 11 GB — 500 per-seed BEDs (5.34 GB) plus a
6.37 GB `consolidated_random_data.csv` — on a filesystem with **2.4 GB free**. Compaction is what made
the rest of the plan fit.

Each per-seed file is an *unordered multiset* of 4-tuples `(score, subfamily, family, class)`, one row
per intersecting TE in arbitrary order. Row order carries no information, so run-length encoding by
tuple is **lossless for every downstream statistic**. Measured on seed 1: 545,099 rows → 70,040 unique
tuples; 10.7 MB raw → 217 KB at `zstd -19`.

| | |
|---|---|
| Result | 5.34 GB → **108.7 MB, 49.2×**, in ~12 min (not the estimated 1–2 h). Free disk 2.4 GB → **14 GB** |
| Schema | `seed, score, subfamily_name, family_name, class_name, n` — one zstd-19 file per seed. **`n` is a weight**, not a row count. `family_name` is empty for a few thousand rows per seed, faithfully to the source, so read with explicit tab separation and NA-keeping off |
| Verification | every seed passed three independent checks — row count, per-class totals, and byte-exact md5 of `sort <source>` vs `sort <reconstruction>` — **before** its source was deleted |
| Aggregate check | all 271,486,562 rows of the consolidated CSV: per-(permutation, class) totals, per-(permutation, family) totals and every per-(class, score) histogram bin identical, so per-class divergence mean/SD/deciles are identical by construction. Only then was the CSV deleted |
| Provenance | `MANIFEST.json` per store: script, version, creation time, source directory, window, seed list, byte counts, ratio, aggregate-check result |
| Reconstruct | `01c_expand_counts.py --seed N` rebuilds an exact legacy BED, so the old format is never truly lost |
| Other windows | `permutation_counts_5kb/` 81 MB, `permutation_counts_20kb/` 135 MB — gitignored, regenerable in ~16–18 min each |

Fallback had any check failed: stop, keep the source, `zstd -19` the BEDs in place (~6×). Never delete
a source whose verification did not pass (C2).

**Publication choice:** the 10 kb store is tracked (106 MB) because every headline number rests on it;
the 5 kb and 20 kb stores are excluded because they are regenerable and every table derived from them
*is* tracked, so no claim becomes unverifiable. Stated in `.gitignore` and `REPRODUCE.md` §1.4 rather
than left to be discovered.

---

## 5. What each analysis script found

### `02_ifna_domain_test.py` — Major 2

Domain chr9:21,150,692–21,370,055 (220 kb). Every descriptive number in the plan reproduced exactly:
**175 TEs, 77 L1 (44 %)** from **36 distinct L1 subfamilies**, family breakdown Alu 33 / L2 15 /
MIR 13 / hAT-Charlie 11 / ERV1 10 / ERVL-MaLR 6 / others 10; mean L1 divergence **135.7** vs
genome-wide 188.2 (median 197); **351 vs 181 L1/Mb = 1.9×**; 12 genes, and the 12 are the IFNA cluster
(IFNA4/5/6/7/10/14/16/17/21, IFNA22P, IFNW1, KLHL9).

All four tests are significant, including the one the plan flagged as risky:

| Test | Result |
|---|---|
| T1 divergence vs unmatched genome-wide background | 135.7 vs null 189.4 ± 23.6, z = −2.28, **p = 0.022** |
| T2 vs L1-count-matched windows (≥ 40 L1) | 135.7 vs 189.2 ± 17.4, z = −3.07, **p = 0.0061** |
| T3 vs gene-density-matched windows (≥ 10 genes) | 135.7 vs 204.5 ± 22.8, z = −3.02, **p = 0.0017** |
| T4 subfamily composition (2×2 Fisher) | 38 young primate-specific vs 39 old L1M\* in-domain against 133,450 vs 412,209 genome-wide, **OR = 3.01, p = 3.2 × 10⁻⁶** |

Robustness — this is what answers "not driven by a few outliers": dropping the divergence-0 element
gives p = 0.0072, dropping the five youngest gives p = 0.0142, and the leave-one-out mean spans only
132.7–137.5, every value below the null mean of 189.2.

**C7 resolved.** Exactly one element, an L1P3 at chr9:21,356,054–21,362,200, has divergence 0.
Genome-wide, 545 of 565,459 L1 (0.096 %) share that annotation, so it is a known RepeatMasker artefact
class rather than a unique corruption, and no conclusion rests on it. Disclosed in the text with the
leave-one-out means.

**Honest limit:** T3 draws 3,582 nulls rather than 10,000 because only ~3 % of the genome has ≥ 10
genes in 220 kb. Recorded in `n_null_windows` / `n_qualifying_in_pool`, printed on the figure panel
itself, and p = 0.0017 sits well above that test's 0.00056 floor. Enlarging the candidate pool to
reach 10,000 is a cheap improvement worth doing before the figure is final.

**Framing fix.** The manuscript conflated two quantities; both now appear, separately labelled — the
per-gene mean LINE divergence across the 8 IFNA TSS windows (95.0–161.7), and the mean divergence of
the 77 L1 elements in the 220 kb domain (135.7). That conflation is what Major 2 was reacting to.

### `04a_lu2020_geneset_overlap.py` — Major 4

**The plan's premise was wrong in one material way: Lu et al.'s categories are *mouse* gene sets**,
defined on mm9 with mouse SINE subfamilies B1/B2/B4 — 1,480 L1-enriched, 2,439
low-complexity-repeat-enriched, 2,041 SINE-enriched, matching their Results text exactly. A comparison
"on the same dataset" was therefore never available; it has to cross the species boundary through MGI
homology, which is a limitation to state rather than a choice. The mapping rate is itself a result:
83–84 % of their SINE and low-complexity genes have a unique human ortholog but only **44 % of the
L1-enriched set**, which is dominated by rodent-expanded families (Zfp\*, Gm\*, Trav\*, Cyp2b\*).
Testing universe = our 28,738 background genes ∩ human genes with a mouse ortholog = **18,632**.

The comparison concords exactly where it should: our Alu family × their SINE-enriched genes **2.87×
over chance** (OR 4.02, n = 273, FDR 4.5 × 10⁻⁶¹); SINE class × SINE-enriched 2.71× (FDR 6.1 × 10⁻⁵⁵);
L1 family × L1-enriched **2.22×** (OR 2.44, n = 49, FDR 6.3 × 10⁻⁷). The cross-terms are significantly
**depleted** (SINE class × L1-enriched 0.13×, OR 0.12) — both studies separate L1-proximal from
SINE-proximal genes the same way, which is the strongest form the answer to Major 4 can take without
a reimplementation. One divergence reported honestly: our **MIR** family does not recapitulate their
SINE-enriched set (0.46×, significantly depleted) — MIR is an ancient tRNA-derived SINE, unlike the
young rodent B1/B2.

Obtaining their Table S1 was itself work: `mmc2.xlsx` came via the Elsevier PII from the Crossref
`alternative-id` (a guessed PII resolved to a different Cell Reports paper entirely). PMC now gates
supplementary downloads behind a JavaScript proof-of-work and cell.com 403s every non-browser client.
All four dead routes are in `external/PROVENANCE.md`.

### `04b_newly_resolved_regions.py` — the assembly ceiling (C6)

D6 drops the assembly-controlled experiment, so the "methodology or assembly?" question is answered by
bounding the assembly contribution descriptively, from the UCSC hs1 → hg38 liftOver chains (13,198
chains, 856,823 blocks).

**The obvious definition had to be rejected:** the raw complement of aligned blocks returns 6.81 % of
the genome but puts 63 % of all TSS windows "in" it, because a 10 kb window almost always contains at
least one 1-bp alignment indel. Corrected definition = (outside every chain span) ∪ (unaligned stretch
≥ 1 kb inside a chain) = **208.75 Mb, 6.70 % of chm13v2.0**, which agrees with the ~182–189 Mb of novel
sequence Nurk et al. 2022 report — the independent check that the chain parse is right. Small-indel
sequence excluded by the floor (3.6 Mb) is reported separately.

The bound is small: at most **0.41 % of TEs** (15,381 of 3,709,429), **0.49 % of TSS windows** (190 of
38,704) and **0.55 % of genes** (158 of 28,738) sit in sequence no hg38 study could have analysed.
**Strongest single sentence for the response letter: the IFNA domain contains 0 bp of newly resolved
sequence — it is entirely alignable to hg38 — so the headline finding cannot be an assembly artefact.**

Newly resolved sequence carries 0.06× the genome-wide TE-annotation density because it is satellite
array and satellites are not one of the six classes analysed; both densities are written out so the
low count is not mistaken for a pipeline error. Most affected class: Retroposon/SVA 2.57 %, then SINE
0.60 %. Most affected chromosomes: chrY 57.5 %, chr9 20.6 %, chr22 18.8 %, chr16 15.3 % — the
acrocentrics, the chr9 heterochromatic block and chrY, as expected. Quoted as a **ceiling**, never an
estimate: liftOver chains are quality-filtered, so hard-to-align sequence counts as newly resolved too.

### `05b_window_sensitivity.py` — Major 5, windows

Two independent confirmations that the window reconstruction is right: the merged base-pair span of the
10 kb set is **exactly 272,233,268 bp**, the published `N_TSS` to the digit; and the 10 kb
observed/random OR column reproduces all six of Table 1's fold changes exactly (LINE 0.877, LTR 0.667,
SINE 1.468, DNA 0.938, SVA 1.368, Helitrons 0.661). The TSS definition was recovered unambiguously:
38,700 of 38,704 published intervals are exactly 10,000 bp and the 4 exceptions all start at 0, so
`TSS = end − 5000` for every row. `N_TSS` is recomputed per window (5 kb 144,952,895 bp, 20 kb
494,969,139 bp) rather than reused.

| Level | 5↔10 kb | 5↔20 kb | 10↔20 kb | Significance flips |
|---|---|---|---|---|
| Classes (n = 6) | ρ = 0.943 | ρ = 0.943 | **ρ = 1.000** | 1 of 6 |
| Families (n = 44) | ρ = 0.891 | **ρ = 0.828** | ρ = 0.941 | 10 of 44 |

All six concordance permutation tests give p ≤ 0.009, so the agreement is measured rather than
eyeballed. Two things are reported honestly rather than smoothed over, and the plan pre-committed to
that framing:

1. Family-level 5↔20 kb agreement is **0.828**, below the ≳ 0.95 expected. The ranking is still
   strongly preserved, but the 5 kb and 20 kb views are not interchangeable at family resolution.
2. The observed/random OR rises **monotonically with window width** for almost every group — SINE
   1.237 → 1.468 → 1.744, SVA 1.075 → 1.368 → 1.665, LINE 0.811 → 0.877 → 0.922. Enrichment magnitude
   is scale-dependent even where ranking is not, and 9 of the 10 family flips go non-significant →
   significant at 20 kb, largely a power effect (20 kb carries ~3.4× the elements of 5 kb). Every
   headline claim's direction is unchanged; its effect size is not. The one class-level flip is
   RC/Helitrons, significant only at 20 kb — consistent with it being the smallest class (1,869
   elements).

### `05c_percentile_sensitivity.py` — Major 5, percentiles

Reuses `06_go_rerun_fdr005.py`'s gene-set constructions with a single parameter changed, so the 5 % arm
of the comparison **is** the published arm and no construction drift can masquerade as a percentile
effect; the script asserts every 5 % set is a strict prefix of its 10 % set before comparing. 60 GO
studies, ~9 min.

**8 of the 9 abstract-level claims survive at 10 %:** type I interferon receptor binding
(0.00050 → 0.0068), olfactory receptor activity (4.8 × 10⁻¹⁵ → 3.7 × 10⁻³¹), SVA → termination of RNA
Pol II transcription (identical at 0.0203), L1 flavone metabolic process (0.0309 → 0.0272), both MIR
zinc terms, SINE mRNA splicing. **The one exception needs Daniil's decision: hAT-Charlie "MHC class I
protein complex" is 5 %-only** — and since WP6 already cuts hAT-Charlie from three terms to one, that
single surviving term is also the one not robust to the percentile choice. Either soften or drop it;
it should not be reported as if the percentile were immaterial.

Term counts rise with the wider cut (classes-by-count 425 → 779, classes-by-divergence 414 → 632,
families 140 → 321), so Jaccard is deflated by gains rather than instability (median 0.29–0.42) while
the meaningful metric — fraction of published terms preserved — is **0.85–0.93**. 144 terms lost and
902 gained, each named in `percentile_sensitivity_terms.csv`. Three families (SVA, ERVK, Gypsy) have
fewer genes than the truncation limit and are **percentile-invariant by construction**, as are SVA and
Helitron at class level; they are flagged as invariant rather than reported as 100 % stable (an earlier
version listing their terms as "lost" was fixed).

`nb05` adds gene-set stability (overlap coefficient 0.29–0.93, every hypergeometric p ≤ 3.5 × 10⁻²⁰⁵)
and rank stability (Kendall τ 0.48–0.78, 200-resample bootstrap CIs). The weakest cell is
Retroposon/SVA between 5 and 20 kb (τ = 0.48, overlap 0.29), expected for the smallest class.

**The honest gap recorded here — GO never re-run per window, so the claim table's 5 kb and 20 kb
cells read `not evaluated` — is now closed.** `07a_build_gene_tables.py` builds the per-window
equivalent of `TEs_on_genes.csv` (its `--verify-10kb` gate proves the 10 kb rebuild reproduces the
published table exactly, down to the per-window multiset of family names), and `07b_go_grid.py` runs
the four missing cells of the 3 × 2 grid by **calling the existing gene-set builders** with a
different `df`, so no third copy of the construction exists and a window effect cannot be confused
with a construction artefact. The two 10 kb cells are copied, not recomputed, and
`07b --check-reuse` proves it two ways: the six copied files are row-identical to their sources, and
a full re-run of `classes_count` at 10 kb from the rebuilt gene table reproduces the published table's
425 (group, term) pairs with zero differences.

What the grid measures (`nb07`, `output/go_grid_summary.json`), and it is not the tidy result:

* **Widening the percentile always finds more terms** (9 of 9 level × window combinations): more
  genes, more power. **Widening the window does not**, and the direction is not even the same at
  every level — classes-by-count *falls* 510 → 425 → 299, classes-by-divergence *peaks at 10 kb*
  (209 / 414 / 303), families *rises* 137 → 140 → 201. A wider window adds elements to every gene,
  so it does not simply add power: it changes which genes are in the top 5 % and dilutes the
  promoter-proximal signal with more distal sequence.
* **Preservation of the published terms is markedly weaker across windows than across
  percentiles**: at worst **0.440**, median **0.677**, against 0.85–0.93 for the percentile arm. So
  the TSS window matters more to GO than the gene-set percentile does.
* **3 of the 9 headline claims survive all six conditions**; all 9 hold in the published cell. The
  two weakest are `SVA / termination of RNA polymerase II transcription` (2 of 6) and
  `hAT-Charlie / MHC class I protein complex` (1 of 6) — the latter reinforcing the decision already
  flagged in §8.
* Concordance of per-group term counts against the published cell: lowest Spearman ρ **0.614**,
  every label-shuffling permutation p ≤ **0.022**.
* **The grid is GO only.** No permutations were re-run, so a difference between cells is a gene-set
  effect, not a background effect, and the enrichment odds ratios of Table 1 and Figure 2 are
  unchanged. The Results text must say so (edit M4) or the two analyses will be conflated.

### `06_go_rerun_fdr005.py` — Minor 2

Validated against the published tables before anything shipped: `classes_count` reproduces
`GO_top_5_perc_genes_by_class_number_with_all.csv` with all 504 (class, term) pairs identical;
`families` reproduces `top_5_perc_genes_by_families.csv` with all 195 pairs identical and the same max
FDR (0.0984); subfamilies return the same 230 non-empty subfamilies with identical term sets — zero
differences. Each level is written twice, `GO_<level>_fdr01_reference.csv` (retrieval cut, reproduces
the published set, **not for publication**) and `GO_<level>_fdr005.csv` (what ships), so the threshold
effect is measured rather than asserted. 1,129 subfamilies processed — exactly the "1,129 of 1,143 that
appear in at least one TSS window" figure. ~20 min instead of the ~9 h the frozen code path needed;
that is what the DAG/GAF process cache buys.

**Plan §3.2 was wrong in one row and silent in another:**

| Level | FDR < 0.1 | FDR < 0.05 | Change | Loses all terms | Plan said |
|---|---|---|---|---|---|
| Classes by count | 504 | 425 | −15.7 % | none | ✓ |
| Classes by divergence | 516 | 414 | −19.8 % | none | ✓ |
| Families by count | **195** | **140** | **−28.2 %** | Dong-R4 | 196 → 160 ✗ |
| Subfamilies | 1,231 | 1,003 | −18.5 % | **31 subfamilies** | not tabulated |

The family row is ground truth from the published `top_5_perc_genes_by_families.csv` itself (195 rows,
140 with FDR < 0.05), so "14 families" → **13** stands. The subfamily row matters for the companion
paper (C3), not for G3: 31 subfamilies drop out entirely at 0.05.

### `10_tables.py` — Minor 4

Reproduces the submitted Table 1 exactly on all six classes and all six checked quantities, and
**refuses to emit the split tables if any value moves**. Two things had to be reverse-engineered:
(i) the published 2×2 is the mixed count/base-pair table from frozen cell 55
(`N_Genome = 3,117,275,501`, `N_TSS = 272,233,268`), not a count-only table, and the class counts come
from the class column directly — aggregating the 44 curated families under-counts LTR (49,741 vs
51,103) and DNA; (ii) the published empirical p is `2 × min(q, 1−q)` clamped to 2/501, which is why
every class shows 0.004. Also emits `TableS_class_enrichment_full.csv` keeping the raw Fisher p (G3
statistics policy) and fixes the malformed `9.3*10⁻¹³³` notation (G8).

### `12_build_trackhub.sh` + `12a_trackhub_beds.py` — Editor request (D8)

10 bigBeds, 105 MB total, largest 42.5 MB (SINE), so every file clears GitHub's 100 MB limit. Tracks:
one bigBed per TE class coloured with the project palette, the 38,704 TSS windows, the TE-top and
TE-bottom gene sets, and the IFNA domain, plus **11 per-track description pages** so a reviewer
clicking a track name gets a real explanation. `12a` **imports `06_go_rerun_fdr005.py`** to build the
gene sets rather than re-deriving them, so the hub cannot disagree with the GO analysis.

`hubCheck` passes; the only remaining message is `Can't get default spec from host genome.ucsc.edu`,
which is `hubCheck` wanting a local `hg.conf`. The script tolerates that one message **by name** and
fails on anything else. Verified independently with `bigBedInfo`: every file opens and its record
count matches the source (1,005,214 LINE / 1,706,485 SINE / 531,410 LTR / 458,177 DNA / 6,274 SVA /
1,869 Helitron = 3,709,429; 38,704 windows).

Two plan corrections: `T2T_genes.bed` windows are **not clipped at chromosome ends**, which
`bedToBigBed` rejects outright — the build clips and reports the count rather than silently
correcting; and divergence runs **0–480**, not the "0–1000, already in that scale" the plan claimed,
so `trackDb.txt` says what the number is instead of implying a UCSC-style score.

One-click URL:
`https://genome.ucsc.edu/cgi-bin/hgTracks?db=hs1&hubUrl=https://nikit357.github.io/T2T_transposons_genes/trackhub/hub.txt&position=chr9:21150692-21370055`

---

## 6. The manuscript work

### 6.1 Citations — read before any text edit

| Fact | Value |
|---|---|
| Mechanism | **Mendeley Cite** (the Word *web add-in*), not the legacy COM plugin |
| In-text citations | **128** `<w:sdt>` content controls tagged `MENDELEY_CITATION…` |
| Bibliography | **1** `<w:sdt>` tagged `MENDELEY_BIBLIOGRAPHY`, flagged `IS_DIRTY` |
| Payload | one `MENDELEY_CITATIONS` property in `word/webextensions/webextension1.xml` (2.57 MB) |
| Style / locale | NLM/Vancouver numeric, `en-GB` — both must change for G3 (G1) |
| Legacy field codes | **none** — zero `ADDIN CSL_CITATION`, `instrText`, `fldSimple`, `fldChar` |

Three consequences: moving or deleting a citation means moving the entire `<w:sdt>` subtree; the second
document must be built by **copying** the manuscript and deleting, never by pasting into a fresh file
(the payload lives in a part a fresh document does not have); and G1 is a Mendeley *style switch*, not
a text edit — what remains ours is the in-sentence grammar that numbered citations allowed and
author–year does not.

### 6.2 Tracked changes

`13_manuscript_tracked_edits.py` reads the untouched baseline and writes the revision, so it is
idempotent and the tracked diff is always against exactly what the journal received. **56/56 edits
located and applied**; validation passes every `word_rewrite` Step 6 check plus the Mendeley checks —
282 `<w:ins>`, 999 `<w:del>`, zero `<w:t>` inside a `<w:del>`, **Reject All reproduces the baseline
text byte-for-byte (102,890 characters)**, all 129 content controls survive.

Two helpers had to be written because the skill's `tracked_replace` rebuilds a paragraph from its
concatenated text and would destroy any `<w:sdt>` inside it — exactly the failure §6.1 warns about:
`tracked_replace_safe` edits only contiguous `<w:r>` blocks and never reaches into a citation control;
`delete_paragraph_safe` marks the runs *inside* controls deleted when a whole sentence goes. Matching
also had to tolerate what a plain-text view hides: non-breaking spaces before citations, zero-width
spaces left by Mendeley, curly quotes, en dashes, and superscript exponents stored as plain digits
(`10-40`, not `10⁻⁴⁰`); `<w:proofErr>` markers had to be stripped because they fragment otherwise
contiguous runs.

`15_house_style.py` adds Phase 7 on top: **21/21 edits**, taking the totals to 395 `<w:ins>` /
1,124 `<w:del>`. Two mechanics it discovered, both now handled generally: after a tracked rename the
skill's `text_of` concatenates deleted *and* inserted text (`RESULTS` → `Results` reads as
`RESULTSResults`), so every lookup works on an `accepted_text()` view; and `run_blocks` deliberately
cannot see inside a `<w:ins>` — which is what protects the citation controls but also hides all the
text Phase 5 added, so edits to our own inserted paragraphs need a separate insertion-aware pass.

**One structural limit, stated in editor notes inside the file:** the two section relocations are
structural moves, not Word tracked moves (`<w:moveFrom>`/`<w:moveTo>` needs paired range bookmarks
that are easy to get subtly wrong). Reject All therefore restores the original *text* but not the
original section *order*, which is why validation compares paragraphs as a multiset.

#### The working file changed on 2026-08-04, and the contract with it

Daniil edited the revision by hand and **accepted 594 of the 1,124 tracked deletions**, resolving one
editor note. `T2T_genes_article_G3_revision_260804_manual.docx` is therefore the current state of the
paper — measured at 129 `<w:sdt>` / 380 `<w:ins>` / 530 `<w:del>` / 9 editor notes, so all 129
Mendeley controls survive and nothing is broken. `T2T_genes_article_G3_revision_260803.docx` is kept
as the record of what the scripts produced before that pass.

`13_…`'s idempotency came *entirely* from `shutil.copyfile(BASELINE, TARGET)`, so running it unchanged
against the working file would have discarded every acceptance. It now **edits in place**, and its
idempotency is per-edit instead: `applied` / `already present` / **`not found`**, with only the last
fatal. Three mechanics that turned out to be needed, none of them obvious from the plan:

1. **Insertions are not self-limiting.** A text replacement cannot apply twice — its search string is
   gone. An insertion can, and a re-run would have added a second copy of the 16 new Discussion
   paragraphs, a third Table 1, and a duplicate repository-URL paragraph. Each insertion now names a
   marker sentence it is checked against.
2. **13's own search strings had to learn 15's rename.** `15_house_style.py` rewrites
   `Supplementary File 4` → `File S4` (G2). Against a document 15 has already touched, several of
   13's targets therefore no longer matched even though the edit was still wanted — which is exactly
   what produced the three `NOT FOUND` results on the first in-place run. Locating and matching now
   try both forms. This is *not* folded into `soften()`: that function returns an index map used to
   splice runs, and a length-changing rewrite would corrupt it.
3. **A skipped insertion must return the paragraphs that are already there**, not an empty list — the
   Lu et al. paragraph is anchored to the return value of the robustness insertion, and an empty list
   made it report a false failure.

Verified state after the conversion: **49 edits, 2 applied, 47 already present, 0 not found.** The two
applied are a locate diagnostic and the Lu et al. overlap paragraph, which is genuinely absent from
the working file.

**`Reject All` no longer reproduces the baseline**, and this is by design rather than a regression:
accepted deletions are permanent. The clean-plus-tracked submission pair is produced from the working
file forward, and the diff a reviewer sees covers the edits made *after* the acceptance pass, not the
whole revision. If the journal wants a full tracked diff against the submitted version, produce it as
a Word **Compare** of the baseline against the final clean file.

### 6.3 The Discussion restructure (WP8/D7)

3,970 w → **2,444 w** (a 38 % cut) in seven subsections: Principal findings; Comparison with prior
work; Limitations; Hypotheses for future testing; Mechanistic framework; **3.5 Connection of TE
enrichment with cancer (kept per D7, condensed)**; Proximity as a null model for TE–epigenomic
studies. Four old headings and 32 superseded paragraphs are tracked-deleted. The bounded assembly
argument is written into "Comparison with prior work" with the measured numbers. One deviation, flagged
in an editor note: the four reviewer-requested subsections lead, but Mechanistic framework keeps its
original position ahead of the cancer and proximity-null subsections, because moving it means
relocating the Figure 9 image paragraph — a one-step drag in Word that changes no text.

`Extensive_discussion_260803.docx` is built from the **baseline** rather than the revision (in the
revision those paragraphs are already marked deleted): 52 paragraphs, 4,106 w, five thematic sections,
each with a written opening and closing so it reads as an essay rather than a pile of orphans. **78 of
the 128 citation controls travel with the text and stay live**, the `MENDELEY_BIBLIOGRAPHY` control is
retained so a refresh rebuilds a bibliography for just these references, and the helper-generated
`<w:ins>` markup is flattened because a standalone file must not read as one unaccepted insertion.

### 6.4a Three findings from the 2026-08-05 review pass

1. **Figures 2D–2F were at the wrong analysis level and had two letters swapped.** They plotted the
   subfamily table (n = 1,129); the published caption is family-level, and it assigns
   D = observed/random OR, E = observed OR, F = random OR — the reverse of what was drawn. This is the
   **third** instance of the same defect class, after Figure 7 and Supplementary 8B: a panel matched
   to a letter by content rather than against the caption the journal already holds. At family level
   the panels reproduce the three published significance sentences exactly (2D SINE–DNA only; 2E
   LINE–DNA and SINE–DNA; 2F LINE–SINE, LINE–DNA, LTR–SINE), so **no manuscript text changes** — and
   `nb03` cell 17 now asserts those three sentences, which turns the check from a habit into a test.
2. **Figure 7H did not apply the ribbon filter its own published caption promised.** The caption says
   *"Connecting ribbons were filtered by at least 5 GO terms. This filtering was applied to the
   visualization only."* The code applied the threshold to whether the count *label* was printed;
   every ribbon was still drawn. So this was a conformance bug, not a new request — and the caption
   settles the design question it raises: bar heights stay unfiltered, so retained ribbons do not fill
   their bars. Filtering hides 36 class → group and 50 group → family ribbons, 146 GO terms. The
   unfiltered version is Supplementary **S8C**, whose published caption already describes it with no
   filter sentence.
3. **Supplementary Figures S9–S14 have no captions, and S12 did not exist.** Only S1–S8 are captioned
   in the working file. S9, S10, S11, S13 and S14 are cited — inside the Figure 4/5/6 legends and the
   sensitivity paragraphs — but were never captioned. S12 was neither cited nor captioned: the WP14.7
   numbering map had reserved it for a Lu et al. overlap figure that was never produced. So the GO
   grid takes **S12**, the inventory becomes a contiguous S1–S14, and edit M2 is **six** captions
   rather than one. *(Superseded 2026-08-07: layout promoted old S6 to Figure 8A, so the delivered
   inventory is a contiguous **S1–S13** and the GO grid is **S11**. `svg/PLACEMENT.md` §0 has the
   full mapping; `16_figure_alignment_edits.py` applied it to the manuscript.)*

### 6.4 Text findings beyond the plan's list

- **Flavone is only half-removed.** GO:0051552 falls out at the **LINE class** level (FDR 0.088) but
  **survives at the L1 family level** (FDR 0.031) — same raw p (3.29 × 10⁻⁵), different BH correction.
  Deleting it everywhere, as the plan said, would have removed a still-significant result. The text
  removes it from the class-level sentence and the Abstract, keeps it at family level, and names the
  five overlapping genes (UGT1A6/7/8/9/10 — UDP-glucuronosyltransferases), which also answers Minor 5.
- **hAT-Charlie also carries olfactory receptor activity** (FDR 0.011), which the published text never
  mentioned; the MHC claim correctly drops from three terms to one (MHC class I protein complex, FDR
  0.0255; β2-microglobulin binding and antigen processing MHC-I both 0.0514, removed).
- **The published "10⁻⁴⁰ – 10⁻⁸⁰" range for the three most significant TE-depleted terms was wrong**
  and is threshold-independent: the actual FDRs are 2.3 × 10⁻⁹¹, 3.3 × 10⁻⁸³ and 2.7 × 10⁻⁸¹.
- Also corrected: LTR loses glutamatergic synapse and LPS signalling; two functional-group
  associations (LINE × lipids, TE-top × DNA replication/recombination) no longer reach 0.05 and are
  explicitly no longer claimed; family-level functional groups move 22 → 21; S8B's groups 26 → 24;
  Figure 7B is now **non-significant** (raw p = 0.113) and 7C is R = 0.645, p = 0.017 (not
  0.633/0.015).
- **G1b turned out to be already done in the main text.** Both offending narrative citations left with
  the Discussion trim, so the main manuscript contains zero narrative numbered citations; both are
  fixed in the extended-discussion file. The fix does **not** follow the plan's wording, deliberately:
  substituting the author name gives "by Lu et al. (Lu et al. 2020)" once Mendeley renders author–year.
  Recasting the citation as **parenthetical rather than narrative** is grammatical under both styles and
  removes the dependency on G1 entirely.
- Abstract is **248 words** against G3's 250 cap and names "human" in the second sentence; keywords,
  running title (`TEs near human genes in T2T`, 27 chars) and the corresponding-author block are in,
  with street address and phone left as an editor note rather than invented.

---

## 7. Verification, with expected values

Plan §10 is the canonical list. The checks that matter most:

- **Freeze:** `md5sum -c` on the four frozen files → 4× OK; `git status --short` on them → empty.
- **No frozen imports:** grep for `nbimporter|import_ipynb|Gene_ontology_analysis|TEs_mapped_on_TSS_analysis`
  across `revision_G3/` → hits only in comments and provenance docstrings (21 such hits).
- **Canonical inputs:** `wc -l` → 3,709,429 and 38,704; the TSS interval sets `diff`-clean against the
  legacy bedGraph.
- **Compaction:** 500 `.tsv.zst` files, ~106 MB, `01c_expand_counts.py --seed 1` reconstructs exactly.
- **Empirical p floor:** minimum of the empirical columns = 0.003992… (2/501), **not** 0.002.
- **Windows:** 500 files each in `permutation_counts_5kb/` and `_20kb/`;
  `05b_window_sensitivity.py --summary` → Spearman ρ > 0.9.
- **FDR:** max FDR < 0.05 in every `GO_*fdr005*.csv`; no class-level `flavone`, `copper ion`,
  `response to cadmium` or `beta-2-microglobulin` in the shipped tables.
- **Networks:** `output/network_qc.json` → `label_collisions == 0` for every panel except the recorded
  S11 waiver (7 overlapping pairs at the rung the waiver replays, down from 11, with all 30 terms per
  family kept), with the fallback rung reached per panel.
- **Network page geometry (S9–S11):** `svg_pt` / `aspect` in the same file, measured from the written
  file → S9 and S10 at 815.04 × 1608.48 pt (aspect 0.507, half the published width at the same
  height), S11 at 804.24 × 1267.92 pt (aspect 0.634, the aspect of reference frame `861:33`).
  `rl.assert_svg_geometry` refuses anything more than 3 % off target and `rl.assert_labels_on_page`
  refuses a label crossing the page edge; 0 off-page labels in all three.
- **Network panels are not byte-reproducible:** `adjust_text` is stochastic and matplotlib stamps a
  `<dc:date>`, so compare panels with `svg/compare_panels.py` (page size, mark count, exact label
  set), never with `md5sum`. The three main-text panels are the exception, and only because they are
  restored from the approved run rather than re-derived: md5 `dee24f6b…` (4A), `b9c18242…` (5A),
  `2382b91c…` (6A).
- **Mendeley:** webextension part present in **both** docx files; 129 controls in the manuscript;
  78 + bibliography in the extended discussion; tracked ins/del > 0.
- **Manuscript greps:** no `T2T_genes_evolution`, no `1000 random`, no `FDR < 0.1`, no
  `REFERENCES`/`ETHICAL STATEMENT`/`ACKNOWLEDGEMENTS` headings; `0.004` still present; `arms race` only
  where citing others' work.
- **Track hub:** `hubCheck -checkSettings` clean apart from the tolerated `hg.conf` message; once
  published, `curl -sI` → 200 on `hub.txt` and 206 on a ranged `.bb` request.

---

## 8. What remains, and who owns it

Everything computational is complete. Open items, all requiring a browser, Word, Figma or an
outward-facing action:

| Item | Owner | Note |
|---|---|---|
| ~~Phase 4b Figma pass~~ **DONE 2026-08-07** — exported to `current_figures_260807/`. Layout renumbered the supplementary set to **S1–S13** (old S6 became Figure 8A) and gave Figure 8 a **fourth panel**; the manuscript was realigned by `16_figure_alignment_edits.py` | Daniil | `svg/PLACEMENT.md` §0 is the as-built map; `figures_text_alignment_plan_260807.md` the reconciliation |
| Figma remainder: panel-letter **A** missing on Figure 4 · `Fisher p = 3.2⁻⁶` malformed in Figure 8C · `-log10(FDR)` colourbar missing from Figures 5A, 6A and Supplementary S9 · S13 panel C legend says `Retroposon`/`RC` not `SVA`/`Helitron` · place `S14C…svg` as **S13D** · rename `Supplementary Fgiure 11.pdf` · re-export afterwards | Daniil | plan §4.3, items F-1 … F-10 |
| ~~Edit the baked-in `"GO terms count in a group (FDR < 0.1)"` in frame `861:34`~~ **DONE** — Figure 7 now reads `FDR < 0.05` and no exported PDF contains `FDR < 0.1`, verified against all 23 files | Daniil | G17/C11 |
| Publish the track hub, in three steps: `12b_publish_trackhub.sh --push` → enable Pages on `gh-pages` (*Settings → Pages → Deploy from a branch → `gh-pages`, `/ (root)`*) → `12c_verify_trackhub_live.sh` | Daniil | not taken unilaterally: publishing to a public branch is outward-facing, so the push sits behind `--push`. **Diagnosed 2026-08-14: the URL 404s because the hub was never published** — `has_pages: false` on the repos API, remote branches `[main]`, no Pages API record. The artefact itself is intact (10/10 bigBeds valid, 105 MB), so nothing needs rebuilding. Full plan in `../trackhub_ghpages_plan_260814.md`; commands in `REPRODUCE.md` §6 |
| Mint the Zenodo DOI (D9) and replace `[ZENODO DOI]` | Daniil | one string in the docx, one in `README.md` |
| G1: switch `MENDELEY_CITATIONS_STYLE` to G3/CSE author–year and the locale to US; then fix the four preprint entries (refs 30, 65, 66, 84) in the Mendeley **library**, not the docx | Daniil | editing the docx would be overwritten on the next refresh |
| Open both docx files in Word, refresh Mendeley, confirm no citation renders as plain text | Daniil | the only reliable check; all XML-level checks already pass |
| Confirm the 10 browser-verify items in `../G3_article_guidelines.md` §13 | Daniil | the live OUP page 403s every automated route (C13) |
| Phase 8: response letter (6 major + 6 minor + 2 editor items), clean + tracked docx, cover letter with the companion-preprint disclosure and the "Reagent Table not applicable" note, upload, tag `g3-revision-1` | — | be explicit where we decline: Major 3, Major 4, Minor 1 |

**One decision left before the text is final** (the hAT-Charlie one is closed — see below):

1. ~~**hAT-Charlie MHC class I**~~ **CLOSED 2026-08-07.** Softened in the Results by
   `16_figure_alignment_edits.py` edit T-8: the sentence now says it is significant only in the
   5 % set *and survives only one of the six window × cut-off conditions*, and reports it as the
   weakest of the family-level associations rather than an established one.
2. **`TEs_on_genes_counts_subfamilies.csv`** — tracked at 85.1 MB, 85 % of GitHub's hard limit,
   compresses 75.5× to 1.1 MB. The verified `.gz` exists; swapping it in means removing a tracked file
   that Data availability names, so the two commands are written out in `REPRODUCE.md` §9 rather than
   run.

**One cheap improvement, noted and not done:** enlarge IFNA test 3's candidate pool to reach the full
10,000 null windows.

**Closed since:** re-running GO per window. The second round did ask, and the full 3 × 2 grid is now
built (§5) — `robustness_headline_claims_by_condition.csv` no longer contains a single
`not evaluated` cell, and `nb05` asserts that.

---

## 9. Caveat index

The plan's §9 carries the full text of C1–C20. The ones that change how work is done here:

| | |
|---|---|
| C2 | The compaction delete is the one irreversible step. Verify before every delete; keep `MANIFEST.json` and `01c` |
| C3 | `GO_tables/`, `top_5_perc_genes_by_*.csv` and the GO summaries are **shared with the subfamilies manuscript**. Tightening to 0.05 silently invalidates its figures — hence the `legacy_fdr01_n500/` snapshot and the deferred refresh |
| C4 | Disclose the 2026b companion preprint in the cover letter — it is load-bearing in the Major 3 answer |
| C5 | Do not re-download `go-basic.obo` / `goa_human.gaf`. Open item: the Methods cite a 2025-12-31 snapshot but the files self-report `rel(2025-10-10)`; **correct the Methods date to match the files**, do not refresh the files |
| C6 | Major 4 is answered without the controlled experiment. A second round could ask for the hg38 re-run; it is ~6–10 h and now affordable on disk |
| C8 | The SSU72L1–L5 + POLR2A cluster spans 116 kb, so the SVA → termination association was expected to weaken at ±2.5 kb. It did not — the term is identical at 0.0203 across percentiles |
| C9 | The Figma canvas is shared with the subfamilies paper and `Figure 9` exists twice. Scope by node ID; a name-based sweep corrupts the companion paper's figures |
| C10 | Mendeley content controls are fragile — see §6.1 |
| C11 | Analysis parameters are baked into figure text and do not update when the CSVs do |
| C12 | House style cannot be parallelised with the text work; it must come last |
| C15 | `python-docx` lives in `~/venvs/collagen_3_11`, not the Retroelements venv |
| C16 | Label-collision checking is **enforced, not attempted**. If the assertion fires, climb the documented ladder (short labels → larger canvas → fewer terms); do not loosen the check |
| C18 | The reviewer never asked for new biology. Every new analysis answers a specific numbered comment; anything else belongs in the subfamilies paper |
| C19 | Copied helpers will drift from their originals. `revision_lib.py` is authoritative for anything G3-related; retire the notebook copies when the subfamilies paper is refreshed |
| C20 | The repository keeps one visible contradiction deliberately: `NUM_PERMUTATIONS = 1000` in the frozen notebook against 500 in the paper. Documented in three files; wording must be unambiguous |
