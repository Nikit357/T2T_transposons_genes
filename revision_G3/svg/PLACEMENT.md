# `revision_G3/svg/` → Figma: what goes where

> **Superseded in part, 2026-08-07.** The figures have been placed and exported to
> `revision_G3/current_figures_260807/`, and layout made two changes this file did not
> anticipate. **The exports are authoritative; §2 and §3 below are the plan, not the
> outcome.** See §0 for what was actually built and
> `../figures_text_alignment_plan_260807.md` for the full reconciliation.

Every SVG in this directory, its destination frame, and what it replaces. Frame IDs are from
plan §5b.1 and are **authoritative** — place into the existing frame, do not re-create it, or
the IDs referenced by the plan, `G3_article_guidelines.md`,
`G3_figure_pvalue_labels_260803.md` and the response letter go stale (caveat C9).

---

## 0. What was actually built (2026-08-07) — read this before §2 and §3

**The supplementary set is S1–S13, not S1–S14.** Old Supplementary Figure 6 (the UCSC view of
the interferon-alpha domain) was consumed by the new main Figure 8, so everything above it
shifted down one. Established by matching every SVG's label set against every exported PDF,
all-vs-all:

| **Final** | Content | **§3 called it** |
|---|---|---|
| S1–S5 | unchanged | S1–S5 |
| **S6** | log-scaled TSS by TE family | S7 |
| **S7** | family gene-set network / clustermap / Sankey (A, B, C) | S8 |
| **S8** | full GO network, classes by element count | S9 |
| **S9** | full GO network, divergence-stratified | S10 |
| **S10** | full GO network, families | S11 |
| **S11** | GO grid (A, B, C) | S12 |
| **S12** | window / percentile concordance (A, B, C, D) | S13 |
| **S13** | gene-set & rank stability (A) + permutation convergence (B, C, D) | S13E + S14 |

**Figure 8 has four panels, not three.** `Fig8C_inset_leave_one_out.svg` was promoted from an
optional inset to a full **panel D**, so Figure 8 is A (browser view) · B (null distributions) ·
C (subfamily composition) · D (leave-one-out).

**The SVG filenames are deliberately NOT renamed.** The notebooks write them, and renaming
breaks re-execution; `S13E_geneset_and_rank_stability.svg` is panel **S13A**, and
`S14A/B/C_*.svg` are panels **S13B/C/D**.

**The colourbar situation is the reverse of §4.4.** In the exports, Supplementary S10 (the
family network) carries a colourbar and a size legend, while **Figure 6A carries neither** —
§4.4 predicted the opposite. Figures 5A and 6A and Supplementary S9 are missing the
`-log10(FDR)` colourbar their captions promise; adding it is an open Figma item.

**Still open in Figma** (see the plan's §4.3): the missing panel-letter **A** on Figure 4; the
malformed `Fisher p = 3.2⁻⁶` in Figure 8C; the three missing colourbars; S13 panel C's legend
saying `Retroposon`/`RC` where panel B says `SVA`/`Helitron`; placing
`S14C_permutation_convergence_families.svg` as **S13D**; and the misspelled
`Supplementary Fgiure 11.pdf`.

**Already done and verified from the exports** — do not redo: the baked-in `FDR < 0.1` in
Figure 7 is gone, Figures 9 and 10 are renumbered, all eight Figure 7 panels are placed, and
Supplementary S10 is a 100 % label match to `S11_full.svg` with all 30 terms per family.

Deep-link any frame: `https://www.figma.com/design/WRNeTzKZObdmAQ8QG1EZlq/Figures-T2T-genes?node-id=<id with ':' replaced by '-'>`

**Nothing in `revision_G3/` writes to Figma.** Placement, panel lettering and export are the
manual step (§5b.3).

---

## 1. Replace a panel inside an existing frame

| SVG | Figure / panel | Frame | Replaces | Why it changed |
|---|---|---|---|---|
| `Fig1D_relabelled.svg` | 1D (all five sub-panels) | `856:7` | the current 1D block | Minor 3: axis said `-log10(p)` for an FDR-corrected column; `ns` convention now stated; suptitle removed (G11) |
| `Fig2D_relabelled.svg` | 2D | `856:8` | current 2D | **now family-level and D/E swapped to match the caption — see §"Figure 2D–2F"**; colourbar and star legend name Mann–Whitney U + BH across the testable class pairs |
| `Fig2E_relabelled.svg` | 2E | `856:8` | current 2E | same |
| `Fig2F_relabelled.svg` | 2F | `856:8` | current 2F | same |
| `Fig3A_relabelled.svg` | 3A | `856:9` | current 3A | K-S test now named and declared **raw** |
| `Fig3B_relabelled.svg` | 3B | `856:9` | current 3B | same |
| `Fig4A_simplified.svg` | 4A | `859:25` | the dense 4A network | Minor 6: 10 terms/group, Jaccard ≥ 0.2, ≥ 5 shared genes, terms ≤ 500 genes, FDR 0.05 |
| `Fig4B_fdr005.svg` | 4B | `859:25` | current 4B | GO counts recomputed at FDR 0.05 (D2) |
| `Fig5A_simplified.svg` | 5A | `861:28` | the dense 5A network | as 4A |
| `Fig5B_fdr005.svg` | 5B | `861:28` | current 5B | as 4B |
| `Fig6A_simplified.svg` | 6A | `861:33` | the dense 6A network | as 4A, at **5 terms/group** after the 1.2× font increase (was 9 — see §4 and §"Packing") |
| `Fig6B_fdr005.svg` | 6B | `861:33` | current 6B | as 4B |
| `Fig7A_fdr005.svg` … `Fig7H_fdr005.svg` | 7A–7H | `861:34` | all eight panels | every panel plots a GO term count, so all move at FDR 0.05 |
| `S1_relabelled.svg` | Supplementary Figure S1 | `856:10` | current S1 | per-class K-S declared raw, 6 tests stated |
| `S2_relabelled.svg` | S2 | `856:11` | current S2 | **stars had no legend at all**; now names Mann–Whitney U + BH across 44 families |
| `S3_relabelled.svg` | S3 | `856:12` | current S3 | same |
| `S8B_fdr005.svg` | Supplementary Figure S8, panel B | `861:32` | current S8B | GO counts recomputed at FDR 0.05; **panel identity corrected — see below**; 10 % smaller per side (review item 7) |
| `S8C_sankey_full_fdr005.svg` | Supplementary Figure S8, panel **C** | `861:32` | current S8C | the class → functional group → family Sankey with **every** ribbon drawn. Figure 7H is the same plot filtered at ≥ 5 GO terms per ribbon, as its published caption says. Note `Figure S8C` is **not** `Figure 8C` (the IFNA subfamily composition panel) |

### Figure 2D–2F — level and panel letters corrected 2026-08-05

The first version of these three panels was matched to letters by content and plotted the
**subfamily** table (n = 1,129, one point per subfamily). The published Figure 2 caption is
family-level and assigns the letters differently — *"Box plot of **observed to random OR** by
TE **families** between TE classes"* is (D), observed OR is (E), random OR is (F). So the level
was wrong **and** D/E were swapped. Both are fixed:

| File | Quantity plotted | Level |
|---|---|---|
| `Fig2D_relabelled.svg` | observed / random odds ratio | 44 families |
| `Fig2E_relabelled.svg` | observed odds ratio | 44 families |
| `Fig2F_relabelled.svg` | random odds ratio (mean of 500 permutations) | 44 families |

**No manuscript text changes.** At family level the three panels reproduce the published
Results sentences exactly, and `nb03` cell 17 asserts it rather than trusting it:

| Panel | Significant class pairs at FDR 0.05 |
|---|---|
| 2D | SINE–DNA only |
| 2E | LINE–DNA, SINE–DNA |
| 2F | LINE–SINE, LINE–DNA, LTR–SINE |

Two presentation consequences of the level change: the `n=` on the heatmap x-axis now reads
1–22 instead of 2–600, and Retroposon and RC have a **single family each**, so they are drawn
as points without a box and the 9 class pairs involving them are annotated `n/a` rather than
`ns`. None of the published claims depends on those pairs.

### Figure 7 and Supplementary 8B — two panel identities corrected 2026-08-04

The first version of these panels was matched to letters by content alone. Checking them
against the **published captions** (which are unambiguous) showed two mismatches, both now
fixed and re-executed in `nb06`:

**Figure 7** — its caption is entirely family-level. Three panels were wrong: C and G were
swapped, and E, F and H had been filled with subfamily-level plots that the caption does not
describe. The panels now match the caption exactly:

| Panel | Caption says | Plot | Statistic |
|---|---|---|---|
| 7A | GO count by obs/random OR | scatter | Pearson R = 0.167, raw p = 0.585 |
| 7B | families with/without significant enrichment, by GO count | box | Mann–Whitney U, raw p = 0.113 |
| 7C | GO count by total TE number in a family | scatter | Pearson R = 0.645, raw p = 0.017 |
| 7D | GO count by average family divergence | scatter | Pearson R = −0.283, raw p = 0.348 |
| 7E | families with/without GO terms, by average divergence | box | Mann–Whitney U, raw p = 0.208 |
| 7F | families with/without GO terms, by obs/random OR | box | Mann–Whitney U, raw p = 0.029 |
| 7G | families with/without GO terms, by TE count | box | Mann–Whitney U, raw p = 1.2 × 10⁻⁶ |
| 7H | Sankey: GO term groups in TE classes (left) and families (right) | Sankey | Fisher per family × group vs its class, BH-corrected |

Values are in `output/figure7_statistics.csv`.

**Panels 7B, 7E, 7F and 7G are now half the size per side** (review item 4, decision D-a) with
the strip dot size unchanged. At that size a title and two-line tick labels do not fit, so the
group sizes and the Mann–Whitney result **move into the figure legend** (manuscript edit M9).
`nb06` prints the four sentences to paste:

```
(B) Mann–Whitney U, raw p = 0.113 (single test); n = 11 (Significant) vs 2 (Non-Significant).
(E) Mann–Whitney U, raw p = 0.208 (single test); n = 13 (Yes) vs 31 (No).
(F) Mann–Whitney U, raw p = 0.029 (single test); n = 13 (Yes) vs 31 (No).
(G) Mann–Whitney U, p < 0.001 (single test); n = 13 (Yes) vs 31 (No).
```

If M9 is skipped, those four panels report a test with no p-value anywhere — the legend edit is
load-bearing, not cosmetic.

**Panel 7H now applies the ribbon filter its caption already promised.** The first version
applied the ≥ 5 GO terms threshold only to whether a ribbon's count *label* was printed; every
ribbon was still drawn. Filtering the ribbons themselves hides **36 class → group and 50 group
→ family ribbons, 146 GO terms** (`output/sankey_ribbon_filter.json`). Because the filter is
visual only, the bars keep their full heights, so retained ribbons do not fill their bars —
that is correct and the published caption already says so ("applied to the visualization
only"). The unfiltered version is the new **S8C**.

**Supplementary Figure 8B** — the caption defines it as *one* clustermap whose rows are the
class-level groups **and** the family-level groups ("TE top, TE bottom, the four classes with
significant GO terms … and TE families with significant GO terms"). It had been regenerated
at subfamily level. Rebuilding it as the combined map at FDR 0.1 reproduces every claim in
the Results, which confirms the identification; at FDR 0.05 it is **18 TE groups × 24
functional groups**. The subfamily-level clustermap belongs to the companion subfamilies
paper, which reproduces it from its own pipeline; it and the three other
`S_extra_subfamily*` panels were **removed on 2026-08-05** (review item 6) because this
manuscript never used them.

---

## 2. Build the new main Figure 8 (WP2) — *the plan; see §0 for what was built*

The old **Supplementary Figure 6** is promoted. Its frame `861:29` is already vector with no
raster fills, so panel A needs no regeneration. **As built the inset became panel D**, so the
figure has four panels rather than three.

| SVG | Panel | Source |
|---|---|---|
| — | 8A | existing UCSC browser view, frame `861:29`, reuse as-is |
| `Fig8B_null_distributions.svg` | 8B | new — null distributions for tests 1–3 with the observed value marked |
| `Fig8C_subfamily_composition.svg` | 8C | new — L1 subfamily composition + per-element divergence |
| `Fig8C_inset_leave_one_out.svg` | 8C inset (optional) | leave-one-out means; the visual form of "not driven by a few outliers" |

There is **no panel D** — decision D4 removed the planned chromatin panel.

Consequent renumbering, **by node ID only** (caveat C9 — `Figure 9` exists twice on the canvas,
once as this frame and once as a loose label belonging to the subfamilies paper):

| Frame | Was | Becomes |
|---|---|---|
| `766:1208` | Figure 8 (mechanisms schematic) | **Figure 9** |
| `861:35` | Figure 9 (integrated map) | **Figure 10** |

---

## 3. New supplementary frames to create — *the plan; the numbers below all shifted down one, see §0*

| SVG | New figure | Note |
|---|---|---|
| `S9_full.svg` | Supplementary Figure S9 | full version of 4A; keep the original 4A frame intact. Portrait **815 × 1608 pt, aspect 0.507** — create the frame at that aspect, not square (§4c) |
| `S10_full.svg` | S10 | full version of 5A; same **815 × 1608 pt, aspect 0.507** |
| `S11_full.svg` | S11 | full version of 6A; **804 × 1268 pt, aspect 0.634**, the aspect of reference frame `861:33`. Frame **`935:11876`** ("S11_full 1", 1072 × 1127, aspect 0.95) already holds the *old* square S11 — re-import and resize it to 0.634 |
| `S13A_or_by_window.svg` | S13A | obs/random OR at 5 / 10 / 20 kb |
| `S13B_bland_altman_families.svg` | S13B | Bland–Altman between window pairs |
| `S13C_go_term_stability_percentile.svg` | S13C | GO term stability, 5 % vs 10 % gene sets |
| `S13D_headline_claims_by_condition.svg` | S13D | each headline claim under **all six** conditions (3 windows × 2 percentiles). Renamed from `..._by_percentile.svg`: it was two columns, it is now six (decision D-b) |
| `S13E_geneset_and_rank_stability.svg` | S13E | gene-set overlap + Kendall τ with bootstrap CI |
| `S14A_permutation_convergence_classes.svg` | S14A | running mean of the random OR per class |
| `S14B_permutation_convergence_sd_classes.svg` | S14B | running SD per class |
| `S14C_permutation_convergence_families.svg` | S14C | drift in SD units, all 44 families |
| `S12A_go_grid_term_counts.svg` | **S12A** | GO term counts across the 3 windows × 2 percentiles grid, one heatmap per GO level |
| `S12B_go_grid_preservation.svg` | S12B | fraction of published terms preserved, and Jaccard, per cell |
| `S12C_group_survival.svg` | S12C | which TE groups keep ≥ 1 GO term in each cell |

S12 answers the GO half of Major 5, S13 the enrichment half, S14 answers Minor 1.

**S12 takes the number 12, not 15.** The WP14.7 numbering map had reserved S12 for a Lu et al.
overlap figure that was never produced — that paragraph ends "provided as a supplementary
table", and no such SVG exists. S12 was therefore neither cited nor captioned anywhere, so the
GO grid fills the gap and the final inventory is a contiguous **S1–S14 with no gap and no
invented fifteenth item**. *(Superseded: layout removed old S6, so the delivered inventory is a
contiguous **S1–S13** — §0. The GO grid is S11, not S12.)*

Note also that **only Figures S1–S8 have captions in the manuscript.** S9, S10, S11, S13 and
S14 are cited but were never captioned, so manuscript edit M2 writes **six** captions
(S9–S14), not one. *(All fourteen captions exist as of edit M2; `16_figure_alignment_edits.py`
then renumbered them to S1–S13 and deleted the old S6.)*

---

## 4. Things to check while the file is open

1. **The baked-in threshold string.** Frame `861:34` (Figure 7) contains the literal
   `"GO terms count in a group (FDR < 0.1)"`. The new panels say `FDR < 0.05`; the frame text
   must be edited to match or the figure contradicts the revised Methods (caveat C11, G17).
   Grep every frame for `0.1`, `500` and `p-value` while you are there.
2. **Vector colourbars (G12).** `Fig456A_colorbar_vector.svg` replaces the pasted
   `-log10(FDR)` node-colour bitmaps in Figures 4, 5 and 6; `Fig4B_colorbar_vector.svg`
   replaces the clustermap count bars in 4B/5B/6B/S8B. The pasted bitmaps are 23 × 58 to
   30 × 175 px and will not clear 350 dpi at print size.
3. **Figure 6A shows 5 terms per group, not 9 or 10.** Read the number off
   `output/network_qc.json`, never from memory. The 1.2× font increase costs this panel four
   terms: at 12 pt the collision checker cannot place 9 labels for 44 family groups at any
   canvas size the ladder tries, and the overlaps it rejects are real (93–348 px², not float
   noise). Say "up to five GO terms per family" in the legend (manuscript edit M10). The other
   panels keep their full term counts: 4A and 5A at 10, **S9 and S10 both at 30** — S10 gained
   one over its earlier 29 when the panels were re-laid out on the pinned canvas (§4c).
4. **S11 is the one panel with the collision check waived**, now at **7 overlapping label pairs
   instead of 11**, with all 30 terms per family kept. At 30 terms across 44 families the label
   field saturates at any page size; the supplementary full network is there for structure.
   `network_qc.json` records the waiver per panel, the collision count at the rung the waiver
   replays (`collisions_at_waived_rung`) and the off-page label count — it was never applied
   silently, and its cost is now a measured number.
   **S11 carries no colourbar and no legends**: reclaiming the legend strip is what gives the
   network the full pinned width, and both legends then landed on top of node labels (measured
   on the first run). Figure 6A is the family-level twin and carries the same node-size scale
   and the same Jaccard line widths; the caption points at it. S9 and S10 keep their node-size
   legend and, like S11, take the `-log10(FDR)` colourbar from
   `Fig456A_colorbar_vector.svg` (G12).
5. **Panel letters.** G3 wants a letter in the upper-left of each panel; none of these SVGs
   carries one, because letters belong to the composed frame.
6. **Fonts (G18) — every SVG was re-exported at 1.2× for this.** Matplotlib sizes are in
   points and SVG consumers render them at 96 dpi, so 1 pt = 96/72 px: the old
   `GLOBAL_FONT_SIZE = 10` pt arrived in Figma as **13.33 px**. `revision_lib.FONT_SCALE = 1.2`
   makes the base 12 pt = **exactly 16 px**, and every hard-coded size in the notebooks goes
   through `rl.fs()` so the whole figure scales together rather than leaving 8 pt axis labels at
   10.67 px against 16 px titles. All five notebooks were re-executed for this; no SVG here is
   at the old scale. Inter everywhere except Figure 1 and old Supplementary 6, which inherit
   Helvetica from the UCSC exports.
7. **Three empty leftover frames** (`422:32`, `260:206768`, `588:1168`) would export as blank
   artboards in a bulk export — delete or skip explicitly.

---

## 4b. Packing: what the review asked for, and what is achievable

The review asked for the network panels to be packed 30–50 % denser. **They cannot be, and the
reason is geometric rather than a tuning failure.** A fixed number of labels at a fixed size
needs a minimum area to avoid overlap; raising every font by 1.2× inflates each label box by
1.2× per side, i.e. **1.44× in area**. Asking for 30 % *less* area at the same time is asking
for the same labels in 0.49× the area they already needed. With the published term counts and a
zero-overlap requirement (caveat C16, which forbids loosening the check) the two requests are
incompatible, and the font increase wins — 16 px legible text is the reviewable requirement,
and canvas area costs nothing in a vector figure Figma rescales anyway.

Measured, from `output/network_qc.json`:

| Panel | canvas vs nominal | terms/group | collisions | best any compaction rung reached |
|---|---|---|---|---|
| 4A | × 1.56 | 10 | 0 | 3 collisions at × 0.90 |
| 5A | × 2.25 | 10 | 0 | 8 collisions at × 0.90 |
| 6A | × 2.25 | **5** | 0 | 9 collisions at × 0.90 |
| S9 | pinned (§4c) | 30 | 0 | — |
| S10 | pinned (§4c) | **30** | 0 | — |
| S11 | pinned (§4c) | 30 | **waived, 7** | — |

The three supplementary panels are no longer laid out by the canvas ladder at all, so
`canvas_area_vs_baseline` is 1.0 for them by construction; their geometry is pinned and asserted
instead (§4c). S11 remains the one panel whose collision check is waived, and that waiver is
recorded per panel in `network_qc.json` rather than applied silently.

Two things were measured and rejected rather than assumed: raising the `adjust_text` forces
makes packing **worse** (at 2× and 5× the defaults every panel gained collisions, because the
solver overshoots and drives labels into each other), and more `spring_layout` iterations
changes almost nothing. What the reworked ladder does buy over the first version is the
*smallest clean canvas* and the *most terms at it* — trying tight label padding alongside the
loose padding at every rung is worth four extra terms on S9 (30 instead of 26).

**The "cannot be packed denser" conclusion above holds for 4A, 5A and 6A, but it was reached
without ever trying two levers, and both of them work.** They are what makes §4c possible, so
this section should not be read as "compaction was exhausted":

- **label wrapping** — nothing in this repository wrapped label text before 2026-08-07, while
  every reference frame does. A 107-character label on one line is 642 pt; wrapped at 18
  characters it is a fraction of that, and the width it frees is width the layout can use.
- **`spring_layout`'s k** — the first ladder tied k to the canvas factor, so it was never
  raised above the published 0.6 at a fixed canvas size. Raising it spreads nodes at constant
  page area, which is precisely "denser packing of the same information".

Neither was applied to 4A/5A/6A, because those three panels are approved and placed and
regenerating them perturbs them (§4c, last paragraph). If the packing of the main-text panels
is revisited, these are the two levers to try first.

---

## 4c. The pinned canvas for S9, S10 and S11

The three supplementary networks were re-laid out on **2026-08-07** to the page geometry the
figures need, rather than to whatever the collision ladder arrived at:

| Panel | measured | aspect | reference | terms/group | collisions |
|---|---|---|---|---|---|
| S9 | 815.04 × 1608.48 pt | 0.507 | half the published width, same height | 30 | 0 |
| S10 | 815.04 × 1608.48 pt | 0.507 | half the published width, same height | 30 | 0 |
| S11 | 804.24 × 1267.92 pt | 0.634 | frame `861:33` is 0.633 | 30 | waived, 7 |

Read off `output/network_qc.json` (`svg_pt`, `aspect`), which is written from a measurement of
the file, not from `figsize`: the panels are saved with `bbox_inches` **off** so `figsize × 72`
is exactly the page, and `rl.assert_svg_geometry` refuses a panel more than 3 % from its target.
A pinned page crops what a tight bounding box would have grown to include, so
`rl.assert_labels_on_page` also refuses a label that crosses the page edge — **0 off-page labels
in all three panels**.

How each panel got clean, in the order the ladder tries things (both levers cost no
information, which is why they come before any reduction of `top_n`):

- **S9** — short labels, no wrapping needed, k = 0.9, tight label padding.
- **S10** — the same at k = 1.2, and it now carries **30** terms per group rather than 29.
- **S11** — labels wrapped at 18 characters, k = 1.2, `same_class_weight` 10.0 → 2.0. The last
  of those is why its GO circles used to clump: it is a layout-only edge between every pair of
  same-class families at ten times the maximum Jaccard weight, so each class collapsed onto a
  point and dragged its 30 terms with it. `min_top_n` is pinned to 30 for this panel so the
  waiver cannot buy legibility with terms — measured, dropping to 20 terms per family would
  have bought only two fewer overlapping pairs.

Two things to know before regenerating any of these:

1. **Label placement is stochastic.** `adjust_text` gives different coordinates on every run —
   measured: two renders with the same `PYTHONHASHSEED` differ — and matplotlib stamps a
   `<dc:date>` into every SVG. So **no md5 gate on these panels can pass**, and the reproducible
   comparison is `svg/compare_panels.py`, which compares page size, mark count and the exact set
   of label strings.
2. **A re-run perturbs the panels that are already approved.** Re-running `nb06` re-lays out
   4A, 5A and 6A too, and on both 2026-08-07 runs Figure 6A came out at **7** terms per family
   instead of the approved 5. The three main-text SVGs and their `network_qc.json` records were
   therefore restored from the approved run (md5 `dee24f6b…`, `b9c18242…`, `2382b91c…`), so M10
   still reads 5. If you *want* the extra two terms, re-place 6A from a fresh run and update
   M10 — it is a real improvement, not noise, but it is a Figma job.

---

## 5. Figures deliberately NOT regenerated

Per plan §WP16 rule 4, a figure that does not change is not regenerated. Untouched: **1A–1C**,
**2A–2C**, **3 (nothing beyond A/B)**, **Supplementary Figures S4–S7**, and both schematics
(which only get renumbered). Their current SVGs and frames stand.
