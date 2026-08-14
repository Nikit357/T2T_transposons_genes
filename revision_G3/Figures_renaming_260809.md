# Figure renaming and panel moves — 2026-08-09

**For manual execution in Figma.** Nothing in `revision_G3/` writes to Figma; this document is the
instruction set. File key `WRNeTzKZObdmAQ8QG1EZlq`, single page `0:1`.
Deep-link any frame: `https://www.figma.com/design/WRNeTzKZObdmAQ8QG1EZlq/Figures-T2T-genes?node-id=<id with ':' replaced by '-'>`

**Scope any Figma work by node ID, never by name** — the same canvas carries the subfamilies
paper's labels and `Figure 9` exists twice (caveat C9).

## Why this is needed

Principle P3: figures, supplementary figures, tables and the panels inside them are numbered in
the order the text first cites them. The audit of
`T2T_genes_article_G3_revision_260807_manual.docx` (legends excluded, elided lists expanded) found
the supplementary series cited in this order:

```
S1  S2  S3  S4  S5  S8  S9  S6  S7  S10  S13  S12  S11
```

Main figures 1–10 and Tables 1–2 are already in citation order and **do not move**.

Two further inputs fix the plan:

* Your 260807 edit cites **Figure S1A/B/C** for the permutation-convergence curves in Methods and
  **Figure S1D/E** for the TE length distributions in Results. That merges the old S13B–D into S1.
  The legends in the docx were never updated to match — the text pointed at panels whose legend
  described something else. This plan closes that gap.
* The rewritten sensitivity subsection cites the three robustness figures in the order
  window → gene sets → GO grid, which is old **S12 → S13A → S11**.

Everything in `revision_G3/Revised_manuscript/T2T_genes_article_G3_revision_260809.docx` — body
citations, panel letters and all thirteen supplementary legends — already uses the **new** numbers
below. The Figma canvas is the only thing left to change.

---

## 1. The renaming map

Read it as: *the frame currently exporting as `Old` becomes `New`.*

| New | Old | Content | Action |
|---|---|---|---|
| **S1 A, B, C** | **S13 B, C, D** | permutation convergence: running mean, running SD, family drift | **move 3 panels into S1, re-letter to A/B/C** |
| **S1 D, E** | **S1 A, B** | TE length distributions, all vs TSS-mapped | **re-letter A→D, B→E** |
| S2 | S2 | divergence ridge plots, all vs mapped | rename only |
| S3 | S3 | length ridge plots by class, full | rename only |
| S4 | S4 | supervenn gene-set intersections | rename only |
| S5 | S5 | TSS distributions by count (A) and divergence (B) | rename only |
| **S6** | **S8** | full GO network, classes by element count | rename |
| **S7** | **S9** | full GO network, divergence-stratified | rename |
| **S8** | **S6** | log-scaled TSS distributions by TE family | rename |
| **S9 A, B, C** | **S7 A, B, C** | family gene-set network / clustermap / Sankey | rename, panels unchanged |
| S10 | S10 | full GO network, families | rename only (number unchanged) |
| **S11 A, B, C, D** | **S12 A, B, C, D** | window and percentile concordance | rename, panels unchanged |
| **S12** | **S13 A** | gene-set and rank stability across windows | **becomes a single-panel figure**; drop the panel letter |
| **S13 A, B, C** | **S11 A, B, C** | GO grid across 3 windows × 2 percentiles | rename, panels unchanged |

Three numbers are **collisions** — do the renames through temporary names, or Figma will end up
with two frames called the same thing:

* S6 ⇄ S8 swap (old S8 → new S6, old S6 → new S8)
* S7 → S9 while old S9 → S7
* S11 → S13 while old S13's remaining panel → S12, and old S12 → S11

Suggested order: rename every affected frame to `TMP_<old number>` first, then to its new number.

## 2. The one figure that changes structure: new S1

New **Supplementary Figure 1** is a five-panel technical-controls figure. Its legend in the
manuscript now reads:

> **Figure S1.** Technical controls on the permutation background and on element length.
> (A) Running mean of the random odds ratio for each TE class as permutations accumulate, as a
> fraction of its value at N = 500; the shaded band is ± 1 %. (B) Running standard deviation of the
> same quantity, on the same scale. (C) Drift of the running mean for all 44 families, expressed in
> units of the final standard deviation. (D, E) Ridge plots for length distribution comparison
> between all (blue) and TSS-neighborhood-mapped TEs (red), (D) for all classes together and
> (E) for individual classes; distributions were compared by two-sample Kolmogorov-Smirnov test.

Actions:

1. Move the three convergence panels out of old S13 into S1, **above** the two existing panels.
2. Re-letter: old S13B → **A**, old S13C → **B**, old S13D → **C**; old S1A → **D**, old S1B → **E**.
3. Old S13's remaining panel A becomes the whole of new S12; delete its `A` label.

Source SVGs, for reference — **the filenames are deliberately not renamed**, because the notebooks
write them and renaming breaks re-execution:

| Panel | SVG |
|---|---|
| new S1A | `svg/S14A_permutation_convergence_mean.svg` |
| new S1B | `svg/S14B_permutation_convergence_sd.svg` |
| new S1C | `svg/S14C_permutation_convergence_families.svg` |
| new S1D, S1E | the two existing S1 panels |
| new S12 | `svg/S13E_geneset_and_rank_stability.svg` |

Note that `S14C_permutation_convergence_families.svg` was still an open placement item as of
2026-08-07 — it needs to be placed as new **S1C**.

## 3. Main figures — no renumbering, but five open corrections

Figures 1–10 keep their numbers. These items were open on 2026-08-07 and are still open; four of
them are things the manuscript text asserts, so they are traceability defects, not cosmetics.

| # | Figure | Item |
|---|---|---|
| 1 | **Figure 4** | panel letter **A** is missing from the frame |
| 2 | **Figure 8C** | `Fisher p = 3.2⁻⁶` is malformed — should read `Fisher p = 3.2 × 10-6` (plain text, not a Unicode superscript) |
| 3 | **Figures 5A, 6A and new S7** | the `-log10(FDR)` colourbar their captions promise is absent. Vector colourbars exist: `svg/Fig456A_colorbar_vector.svg`, `svg/Fig4B_colorbar_vector.svg` |
| 4 | **new S1C** (old S13D) | its legend says `Retroposon` / `RC` where the panel above says `SVA` / `Helitron`. The manuscript uses **SVA** and **Helitrons**; make the panel match |
| 5 | export filename | `current_figures_260807/Supplementary Fgiure 11.pdf` is misspelled, and under this plan that file becomes **Supplementary Figure 13** |

## 4. Export checklist

After the renames, re-export to a new folder `current_figures_260809/` and confirm:

- [ ] 10 main figures, `Figure 1.pdf` … `Figure 10.pdf`
- [ ] 13 supplementary figures, `Supplementary Figure 1.pdf` … `Supplementary Figure 13.pdf`, all
      spelled correctly
- [ ] new S1 has five panels lettered A–E, convergence first
- [ ] new S12 has no panel letter
- [ ] no frame is called by an old number, and no two frames share a name
- [ ] the four corrections in §3 are applied

## 5. What this plan deliberately does **not** do

**It does not regroup panels across the three robustness figures**, and there is a case for doing
so that I am flagging rather than acting on, because it is your call and it costs real Figma work.

As numbered above, new S11 mixes two subjects — window sensitivity (A, B) and gene-set cut-off
sensitivity (C, D) — and new S12 is a single panel. A regrouping by subtheme would give:

* **window sensitivity**: old S12A, old S12B, old S13A → a 3-panel figure
* **GO robustness across window and cut-off**: old S11A/B/C, old S12C, old S12D → a 5-panel figure

That is 12 supplementary figures instead of 13, each internally coherent (P4d, one subtheme per
object), and it removes the single-panel figure. It also means re-lettering seven panels and
rewriting two legends.

**I did not do this**, because it changes what the figures *are* rather than what they are called,
and the manuscript legends have to match whichever you choose. If you want it, say so and I will
rewrite the two legends and the citing sentences to match; the numbering in the 260809 docx would
then be S11 and S12 with S13 removed.
