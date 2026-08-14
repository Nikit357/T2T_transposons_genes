# Review corrections implementation plan — G3 revision, 2026-08-04

**Input:** `../claude_results_review_260804.txt` (Daniil's review of every file, folder and figure in
`revision_G3/`), plus his inline annotations on the first draft of this plan.
**Scope:** 13 figure changes + a global font change, the supplementary package as **five thematic
Excel workbooks**, and a new GO robustness analysis across window sizes × percentiles.
**Manuscript to edit:** `Revised_manuscript/T2T_genes_article_G3_revision_260804_manual.docx`
(Daniil's manually edited file — see §2.5, this changes how `13_…` works).
**Status:** plan only. Nothing below is implemented yet.

---

## 1. Overview

The review splits into three blocks that need different kinds of work:

1. **Figures (13 items + one global change).** Mostly geometry and colour: panel sizes, packing
   density, colour-by-class, and a global **1.2× font increase** so imported SVG text lands at
   16 px in Figma instead of 13.33 px. Two items are not cosmetic at all — Figures 2D–2F are at the
   wrong analysis level, and Figure 7H is missing a filter its own published caption already
   promises. Both are corrections that make the artwork agree with the manuscript.
2. **Tables.** There is no assembled supplementary package anywhere in the repository. The only
   supplementary files on disk are the **April baseline** in `manuscript_figures_supplementary_old/`,
   whose GO-derived members are still at FDR 0.1 and are therefore superseded by decision D2. The
   revision's own new tables exist as loose CSVs in `output/` with no `File Sn` number and no
   directory. The deliverable does not exist yet, and it is now built as **five thematic workbooks**
   rather than fourteen separate files, so a reader has five things to open instead of fourteen.
3. **Analysis.** GO enrichment currently exists at **10 kb / 5 %** (published) and **10 kb / 10 %**
   (`05c`). The four remaining cells of the 3 × 2 grid — 5 kb and 20 kb at both percentiles — have
   never been run, so the robustness statement in the Results is narrower than it claims and the
   headline-claim table honestly reads `not evaluated` in its window columns. This block builds the
   grid, compares all six cells against the published one for all three GO levels, emits supplementary
   panels, and updates the manuscript.

The **key design decision** for block 3: do not write a third copy of the gene-set construction.
`06_go_rerun_fdr005.py`'s three `run_level_*(df)` builders and `05c_percentile_sensitivity.py`'s three
`build_*(df)` builders already take the per-window gene table as an argument, so the new work only has
to supply a *different `df`* — one per window — and call the existing functions. That keeps all six
grid cells built by exactly the code that produced the published arm, which is what makes a
window effect distinguishable from a construction artefact.

---

## 2. Background — five findings that change the work

Each of these was verified against the manuscript and the data **before** writing this plan, because
each one changes what the correct fix is.

### 2.1 Figures 2D–2F: wrong analysis level *and* wrong panel letters

`nb03_relabelled_figures.ipynb` cell 7 plots `enrichment_df_subfamilies` — one box per class, one point
per **subfamily** (n = 1,129). The published caption and the published Results text are both
**family-level**, and they also assign the letters differently:

| Panel | Caption (Figure 2 legend) | `nb03` currently plots |
|---|---|---|
| 2D | Box plot of **observed to random OR** by TE families between TE classes | `Observed_Odds_Ratio` |
| 2E | Box plot of **observed OR** by TE families between TE classes | `OR_Observed_to_Random` |
| 2F | Box plot of **random OR** by TE families between TE classes | `Random_Odds_Ratio_Mean` |

So **D and E are swapped** as well as being at the wrong level. `svg/PLACEMENT.md`'s "Panel-letter
mapping to confirm" table repeats the same wrong mapping and must be corrected with them.

I re-ran the panels at family level (44 families, Mann–Whitney U per class pair, BH across the pairs)
and the result reproduces the published Results text **exactly**:

| Quantity | Significant pairs at FDR 0.05 | Results text |
|---|---|---|
| obs/random OR → **2D** | SINE–DNA (0.024) only | "DNA elements were significantly depleted … compared to SINEs … all the rest pairwise comparisons were non-significant after the FDR correction" ✔ |
| observed OR → **2E** | LINE–DNA (0.032), SINE–DNA (0.040) | "significantly lower for DNA elements compared to both LINEs and SINEs" ✔ |
| random OR → **2F** | LINE–SINE (0.009), LINE–DNA (0.009), LTR–SINE (0.022) | "the highest number of significant differences was found by the random background OR: between DNA elements and LINEs, between LINEs and SINEs and between LTRs and SINEs" ✔ |

**Consequence: no manuscript text changes for this item.** The fix restores agreement instead of
creating a new claim, and the three sentences above become its verification test. This is the third
instance of the same defect class as the Figure 7 and Supplementary 8B panel-identity corrections of
2026-08-04 — panels matched to letters by content rather than against the published caption — so the
implementation must add the caption check as an assertion rather than a habit.

Family counts per class are **DNA 22, LINE 9, SINE 6, LTR 5, Retroposon 1, RC 1**. Retroposon and RC
have a single family each, so they have no distribution: they must be drawn as single points without a
box, and the 9 class pairs that involve them must be annotated `n/a` rather than `ns`. None of the
published claims depends on those pairs.

### 2.2 Figure 7H: the caption already promises the ribbon filter

The published Figure 7 caption reads: *"(H) Sankey plot visualization of GO term groups found in TE
classes (left) and families (right). **Connecting ribbons were filtered by at least 5 GO terms. This
filtering was applied to the visualization only.**"*

`nb06` cell 35 defines `RIBBON_LABEL_THRESHOLD = 5` but applies it only to whether the **count label**
is printed; every ribbon is still drawn. So this is a **conformance bug**, not a new request — and the
caption also settles the design question it raises: *"applied to the visualization only"* means the
bars keep their full heights and the underlying counts are unchanged. Retained ribbons will therefore
not fill their bars, and that is correct and already documented.

Separately, Supplementary Figure S8's caption already describes a panel C: *"(C) Sankey plot of GO
terms count comparison by groups between the large-scale TE groups … and TE families."* — with no
filter sentence. So the current unfiltered output is exactly S8C.

**Naming hazard:** the p-value labelling document already uses the bare string `8C` for **Figure 8C**
(the IFNA subfamily composition panel). The new panel is **Figure S8C**. Every reference must carry
the `S`.

### 2.3 The 10 kb / 10 % GO tables exist — they are just undiscoverable

`05c_percentile_sensitivity.py` already wrote `GO_classes_count_p10_fdr005.csv`,
`GO_classes_divergence_p10_fdr005.csv` and `GO_families_p10_fdr005.csv` (plus `_fdr01_reference`
twins) into `output/`, where they sit among 60-odd other files with no index. What is genuinely
missing is the four window cells. The fix is therefore *partly* a naming and packaging fix: move the
whole grid into `output/GO_grid/` under one systematic name, with an index CSV, and keep a copy of
the existing `*_p10_*` files' content under the grid name so nothing is recomputed twice.

### 2.4 The font factor is a unit conversion, and 1.2 is exact

Matplotlib font sizes are in **points**; SVG consumers render them in **pixels** at 96 dpi, so
1 pt = 96/72 = 1.3333 px. `GLOBAL_FONT_SIZE = 10` pt therefore arrives in Figma as **13.33 px**,
which is exactly what Daniil sees. To land on **16 px** the point size must be 16 × 72/96 = **12 pt**,
i.e. exactly **1.2×**. So the requested factor is not a fudge — it is the right number, and it applies
to *every* font size literal in the notebooks, not only to `GLOBAL_FONT_SIZE`.

There are **101 font-size literals** across the five notebooks and `revision_lib.py`
(`fontsize=`, `title_fontsize=`, `labelsize=`, `annot_kws={"size": …}`, `markersize=` where it scales
with text). Scaling only the global would leave 8 pt axis labels at 10.67 px against 16 px titles.

### 2.5 The working manuscript changed, and it breaks `13_…`'s idempotency contract

The file to edit from now on is
`Revised_manuscript/T2T_genes_article_G3_revision_260804_manual.docx`, which carries Daniil's manual
edits and his acceptances of part of the tracked diff. Measured state of all four documents:

| File | `<w:sdt>` | `<w:ins>` | `<w:del>` | Editor notes |
|---|---|---|---|---|
| `…_submitted_baseline_260418.docx` (read-only) | 129 | 0 | 0 | 0 |
| `…_revision_260803.docx` (script output) | 129 | 395 | 1,124 | 10 |
| **`…_revision_260804_manual.docx`** (**the working file**) | **129** | **380** | **530** | **9** |

All 129 Mendeley content controls survive in the manual file, so nothing is broken — but **594 of the
1,124 tracked deletions have been accepted** and one editor note has been resolved.

`13_manuscript_tracked_edits.py` line 301 does `shutil.copyfile(BASELINE, TARGET)` before editing. Its
idempotency comes entirely from that line. **Running it unchanged would overwrite the manual file's
predecessor and throw away every acceptance and manual edit.** Three consequences:

1. `13_…` gains an **in-place mode** that edits the working file directly (as `15_house_style.py`
   already does) and, per edit, **skips** rather than fails when the target string is already present.
2. `15_house_style.py`'s `TARGET` must point at the same working file.
3. The documented claim *"Reject All reproduces the baseline text byte-for-byte"* is **no longer true**
   for the accepted subset — accepted deletions are permanent by design. The submission pair (clean +
   tracked) is produced from this file forward, and `README.md`, `CLAUDE.md` and `project_overview.md`
   must say so instead of repeating the old claim. `14_build_extensive_discussion.py` still builds
   from the pristine baseline and is unaffected; leave it alone.

### 2.6 Supplementary figures S9–S14 have no captions, and S12 does not exist

Checked in the working file: **only Figures S1–S8 have captions.** S9, S10, S11, S13 and S14 are cited
(inside the Figure 4/5/6 legends and the Results/Methods sensitivity paragraphs) but **no caption was
ever written for any of them**. And **S12 is neither cited nor captioned**: the WP14.7 numbering map
reserved it for a Lu et al. overlap figure that was never produced — that paragraph ends "provided as a
supplementary table", and no such SVG exists in `svg/`.

Two consequences:

- The GO-grid figure takes **S12**, not S15. The inventory then runs S1–S14 with no gap and no
  invented fifteenth item.
- Manuscript edit M2 becomes **six captions** (S9–S14), not one.

---

## 3. Item-by-item disposition

| # | Review item | Disposition | Where |
|---|---|---|---|
| 1 | Fig2D–2F at subfamily level, need families | **Correction** — switch to `enrichment_families_with_random.csv`, **and swap D/E** to match the caption (§2.1) | `nb03` cell 7 |
| 2 | Fig6A packing too sparse | tighter layout: lower `spring_layout k`, lower `adjust_text expand`, shrink canvas; ladder reordered to compact-first | `revision_lib.py`, `nb06` cells 5, 10 |
| 3 | Fig6B 30 % smaller both directions | `figsize (14, 9)` → `(9.8, 6.3)` | `nb06` cell 19 |
| 4 | Fig7B/E/F/G two times smaller, same dot size | `figsize (4.2, 4)` → `(2.1, 2.0)`, `stripplot size` unchanged, statistics move to the caption | `nb06` cell 27 (`go_count_box`) |
| 5 | Fig7H filter ribbons by ≥ 5 terms; current version is S8C | add `min_ribbon_count`; emit `Fig7H` (=5) and `S8C` (=1) (§2.2) | `nb06` cell 35 |
| 6 | `S_extra_subfamily*` not needed — remove everywhere | delete 4 SVGs, delete the 2 producing cells, drop from `PLACEMENT.md` | `nb06` cells 22, 36, 37 |
| 7 | S8B 10 % smaller both directions | `figsize (14, 9)` → `(12.6, 8.1)` | `nb06` cell 21 |
| 8 | S9–S11 denser packing by 30–50 % | same compaction as item 2, applied at `FULL` settings | `nb06` cells 12–14 |
| 9 | S13A: classes subpanel smaller, dots larger; colour by class then window shade | `width_ratios [1, 2.4]` → `[1, 3.6]`, class `s` 30 → 110; new `shade()` helper | `nb05` cell 3 |
| 10 | S13B: colour dots by TE class | join the family → class map | `nb05` cell 5 |
| 11 | S13C: left-side vertical colour annotation of the comparison groups | new annotation strip, colours from Fig 4B/5B/6B row palettes | `nb05` cell 9 |
| 12 | S13D: bigger circles, colour by adjusted p with the GO-network colourbar, tighter packing, **full 3 windows × 2 percentiles grid** | `s` 260 → 700, `coolwarm` on −log10(FDR), row pitch 0.42 → 0.30, 6 columns | `nb05` cell 11 |
| 13 | All figures: font 1.2× larger | `FONT_SCALE = 1.2`, `GLOBAL_FONT_SIZE = 12`, `fs()` helper applied to all 101 literals; re-execute all five notebooks | `revision_lib.py` + all notebooks |
| 14 | No folder with the new supplementary tables | new `supplementary/` package: **five thematic Excel workbooks, Files S1–S5** | new `08_build_supplementary.py` |
| 15 | GO for 5/10 % × 5/10/20 kb + robustness notebook + panels + manuscript | new `07a`/`07b` scripts, `nb07`, **Figure S12**, manuscript edits | §5.6–5.8 |

---

## 4. The two conflicts this creates, and how they are resolved

**4a. Larger fonts fight tighter packing, and both fight the collision assert.**
`assert_no_label_collisions` refuses to write an SVG with overlapping labels (caveat C16), and the
current fallback ladder resolves collisions by *enlarging the canvas* — which is precisely what made
Figure 6A and S9/S10 sparse (they landed at 19.5 × 16.5 in and 25.5 × 22.5 in). Raising every font by
1.2× makes every label box 20 % wider, so re-executing unchanged code would push those panels further
down the ladder — sparser, or fewer terms.

Resolution: **invert the ladder.** Add layout parameters to the two network functions and make the
first rungs *compaction* rungs, so canvas growth becomes the last resort rather than the second:

```
rung 0  short labels, k=0.35, expand=(1.02,1.05), canvas x0.70   <- target: "denser by 30 %"
rung 1  short labels, k=0.35, expand=(1.05,1.10), canvas x0.80
rung 2  short labels, k=0.45, expand=(1.10,1.15), canvas x0.90
rung 3  short labels, k=0.60, expand=(1.15,1.25), canvas x1.00   <- today's rung 1
rung 4+ …existing enlargement and top_n reduction rungs, unchanged
```

`network_qc.json` gains `canvas_area_vs_baseline` per panel so "30–50 % denser" is a **measured**
number in the QC record, not an assertion. If a panel cannot reach ≥ 30 % compaction cleanly, the
achieved value is reported and Daniil decides.

**Sequencing matters:** do the font change **first**, then tune packing against the final font size.
Tuning packing at 10 pt and then raising the font would invalidate every rung.

**4b. `top_n` may move again, and the captions quote it.**
Figure 6A currently shows **nine** terms per family because that is where the ladder stopped, and its
caption says so ("nine"). If the font increase or the compaction changes the achieved `top_n`, the
caption must change with it (M10) — and the six captions being written for the first time (M2, §2.6)
must quote the final values too, including S9–S11's term counts. Therefore: **read every quoted term
count off the final `network_qc.json`** after the last re-execution, and treat every caption edit as
the last step of the figure block, not the first.

---

## 5. Files to change

### 5.1 `revision_G3/revision_lib.py`

**5.1a — font scale (new, near line 68).**

```python
# before
GLOBAL_FONT_SIZE = 10

plt.rcParams["svg.fonttype"] = "none"
plt.rcParams.update({"font.size": GLOBAL_FONT_SIZE})
```

```python
# after
# Matplotlib sizes are in points; SVG consumers render at 96 dpi, so 1 pt = 96/72 px.
# A 10 pt label therefore arrives in Figma as 13.33 px, against the 16 px the frames use.
# FONT_SCALE = 1.2 makes the base 12 pt = exactly 16 px. Every hard-coded size in the
# notebooks goes through fs() so the whole figure scales together, not just the base.
FONT_SCALE = 1.2
BASE_FONT_SIZE = 10
GLOBAL_FONT_SIZE = BASE_FONT_SIZE * FONT_SCALE  # 12.0 pt -> 16.0 px in Figma


def fs(points):
    """Scale a font size expressed at the original 10 pt baseline."""
    return round(points * FONT_SCALE, 2)


plt.rcParams["svg.fonttype"] = "none"
plt.rcParams.update({"font.size": GLOBAL_FONT_SIZE})
```

**5.1b — network layout parameters (`save_go_network_svg` line 421, `save_go_network_svg_families_by_classes` line 586).**
Add four keyword arguments to both, defaulting to today's values so nothing changes until a caller
opts in:

```python
    layout_k=0.6,              # spring_layout k; lower = tighter
    layout_iterations=None,    # None keeps the current 50 / 80
    label_expand=(1.15, 1.25),
    label_forces=((0.6, 0.8), (0.4, 0.6)),   # (force_text, force_static)
```

and thread them through, replacing the two hard-coded call sites:

```python
# before (line 497 / 680)
pos = nx.spring_layout(G, k=0.6, iterations=50, seed=42)
...
adjust_text(texts, ax=ax, arrowprops=..., force_text=(0.6, 0.8),
            force_static=(0.4, 0.6), expand=(1.15, 1.25))
```

```python
# after
pos = nx.spring_layout(G, k=layout_k, iterations=layout_iterations or 50, seed=42)
...
adjust_text(texts, ax=ax, arrowprops=..., force_text=label_forces[0],
            force_static=label_forces[1], expand=label_expand)
```

**5.1c — colour shading helper (new, after `CLASS_PALETTE_SHORT`).** Needed by S13A.

```python
def shade(colour, amount):
    """Lighten (amount < 0) or darken (amount > 0) a hex colour by blending.

    S13A needs three tints of each TE class colour, one per TSS window: 5 kb faint,
    10 kb the palette colour itself, 20 kb dark.
    """
    r, g, b = mcolors.to_rgb(colour)
    if amount < 0:                      # blend toward white
        f = -amount
        return (r + (1 - r) * f, g + (1 - g) * f, b + (1 - b) * f)
    return (r * (1 - amount), g * (1 - amount), b * (1 - amount))


WINDOW_SHADES = {"5kb": -0.55, "10kb": 0.0, "20kb": 0.35}
```

**5.1d — the GO-network colour scale, exported (new).** S13D must use "the same colourbar as the GO
networks". The networks use `coolwarm` on −log10(FDR) but auto-scale per panel, so the shared quantity
must be named explicitly rather than copied by eye:

```python
GO_FDR_CMAP = "coolwarm"


def go_fdr_norm(fdr_values, floor=FDR_THRESHOLD):
    """Normalizer for -log10(FDR) with the publication threshold as vmin."""
    values = -np.log10(np.asarray(fdr_values, dtype=float))
    return mcolors.Normalize(vmin=-np.log10(floor), vmax=float(np.nanmax(values)))
```

Use it in both network functions in place of the inline `Normalize(vmin=min, vmax=max)` so the
networks and S13D really do share one scale, and record the range in `network_qc.json`.

**5.1e — move `CLASS_PALETTE_DIVERGENCE` in from `nb06` cell 3**, keeping its provenance comment
(verbatim copy of frozen GO cell 111). `nb05` needs it for the S13C row annotation and `nb06` keeps
using it; one copy, two consumers.

### 5.2 `revision_G3/nb03_relabelled_figures.ipynb`

**5.2a — cell 7, Figures 2D–2F to family level with the caption's letter order.**

```python
# before
def boxplot_with_significance_matrix(y_name, panel_name, y_label, log_y=True):
    """One Figure-2 panel: subfamily distributions per class + BH-corrected pair matrix."""
    x_name = "class_name"
    frame = enrichment_df_subfamilies[enrichment_df_subfamilies[y_name] > 0].copy()
...
for panel_name, y_name, y_label, log_y in [
    ("Fig2D", "Observed_Odds_Ratio", "Observed odds ratio", True),
    ("Fig2E", "OR_Observed_to_Random", "Observed / random odds ratio", True),
    ("Fig2F", "Random_Odds_Ratio_Mean", "Random odds ratio (mean of 500)", True),
]:
```

```python
# after
def boxplot_with_significance_matrix(y_name, panel_name, y_label, log_y=True):
    """One Figure-2 panel: FAMILY distributions per class + BH-corrected pair matrix.

    The published caption and Results text are family-level ("Box plot of ... by TE
    families between TE classes"), and the letters are D = obs/random OR, E = observed
    OR, F = random OR. The first version of this notebook used the subfamily table and
    the D/E order swapped; both are corrected here. The three assertions at the end of
    the notebook check the panels against the published sentences.
    """
    x_name = "class_name"
    frame = enrichment_df_families[enrichment_df_families[y_name] > 0].copy()
...
# Panel letters follow the published caption, not the order the quantities happen to
# appear in the enrichment table.
for panel_name, y_name, y_label, log_y in [
    ("Fig2D", "OR_Observed_to_Random", "Observed / random odds ratio", True),
    ("Fig2E", "Observed_Odds_Ratio", "Observed odds ratio", True),
    ("Fig2F", "Random_Odds_Ratio_Mean", "Random odds ratio (mean of 500)", True),
]:
```

Three further edits inside the function:

- `Status` does not exist in the families table. Derive it exactly as `nb06` cell 26 does, in cell 3:
  `enrichment_df_families["Status"] = np.where(enrichment_df_families["p_adjusted_empirical_bh"] < 0.05, "Significant", "Non-Significant")`.
- Legend text `subfamily significant (FDR-corrected empirical p < 0.05)` → `family significant …`.
- Singleton classes: draw Retroposon and RC as points only (a filtered `sns.boxplot` call on classes
  with n ≥ 2), and write `n/a` instead of `ns` in the matrix for any pair whose smaller group has
  n < 2. The heatmap x-tick labels already carry `n=`, which now reads 1–22 instead of 2–600 and makes
  the change self-evident.

**5.2b — cell 17, add the caption-conformance assertions.**

```python
# after (new, replacing the existing single invariant check)
EXPECTED_SIGNIFICANT = {                      # published Results, Figure 2 paragraph
    "Fig2D": {("DNA", "SINE")},
    "Fig2E": {("DNA", "LINE"), ("DNA", "SINE")},
    "Fig2F": {("DNA", "LINE"), ("LINE", "SINE"), ("LTR", "SINE")},
}
for panel, expected in EXPECTED_SIGNIFICANT.items():
    got = {tuple(sorted(pair)) for pair in significant_pairs[panel]}
    assert got == {tuple(sorted(p)) for p in expected}, (
        f"{panel} no longer reproduces the published sentence: {got} vs {expected}")
print("Figures 2D-2F reproduce the three published significance statements")
```

**5.2c — fonts.** 38 literals in this notebook go through `rl.fs(...)`.

### 5.3 `revision_G3/nb05_sensitivity_robustness.ipynb`

**5.3a — cell 1, load the family → class map (it is not loaded today).**

```python
# after (append)
family_to_class = (
    pd.read_csv(os.path.join(rl.REPO_DIR, "families_by_classes_TE.csv"))
    .dropna().set_index("family_name")["class_name"].to_dict()
)
```

**5.3b — cell 3, panel S13A.**

```python
# before
fig, (ax_class, ax_family) = plt.subplots(
    1, 2, figsize=(11.0, 5.2), gridspec_kw={"width_ratios": [1, 2.4]}
)
...
    for window in wide.columns:
        ax.scatter(wide[window], positions, s=30, zorder=3,
                   color=WINDOW_COLOURS[window], edgecolor="black", linewidth=0.3,
                   label=window)
```

```python
# after
# The classes panel is deliberately narrow with large dots and the families panel wide
# with small ones, so the six classes read as the headline and the 44 families as the
# supporting detail.
fig, (ax_class, ax_family) = plt.subplots(
    1, 2, figsize=(11.0, 5.2), gridspec_kw={"width_ratios": [1, 3.6]}
)
DOT_SIZE = {"class_name": 110, "family_name": 30}
...
    for window in wide.columns:
        colours = [rl.shade(class_of(name, level), rl.WINDOW_SHADES[window])
                   for name in wide.index]
        ax.scatter(wide[window], positions, s=DOT_SIZE[level], zorder=3,
                   color=colours, edgecolor="black", linewidth=0.3)
```

with `class_of(name, level)` returning `rl.CLASS_PALETTE[name]` at class level and
`rl.CLASS_PALETTE[family_to_class[name]]` at family level. The legend becomes two blocks — TE class
swatches (from `CLASS_PALETTE_SHORT`) and three grey tints labelled `5 kb (faint) · 10 kb · 20 kb
(dark)` — with one sentence in the caption saying the tint encodes the window and the hue the class.
`WINDOW_COLOURS` is kept only for S13C bar fills.

**5.3c — cell 5, panel S13B.** Replace the single `color="#70453c"` with per-point class colours:

```python
# after
    colours = [rl.CLASS_PALETTE.get(family_to_class.get(f, ""), "#cccccc")
               for f in pair.index]
    ax.scatter(mean_of_pair, difference, s=26, color=colours,
               edgecolor="black", linewidth=0.3, zorder=3)
```

plus a shared class legend on the right-hand axes.

**5.3d — cell 9, panel S13C row annotation.** Add a 0.12-inch axes to the left of `ax_counts` and draw
one coloured cell per row, keyed by the row's level so the colours are the same ones the reader sees in
Figures 4B, 5B and 6B:

```python
# after (new helper above the figure)
def row_annotation_colour(level, group):
    """Row colour matching the clustermap side annotations of Figures 4B / 5B / 6B."""
    if level == "classes_count":
        return rl.CLASS_PALETTE.get(group, "#cccccc")          # as Figure 4B
    if level == "classes_divergence":                            # "LINE / highest"
        return rl.CLASS_PALETTE_DIVERGENCE.get(group.replace(" / ", "_"), "#cccccc")  # 5B
    return rl.CLASS_PALETTE.get(family_to_class.get(group, ""), "#cccccc")           # 6B
```

**5.3e — cell 11, panel S13D: the full 3 × 2 grid.** Six columns (5/10/20 kb × 5/10 %) instead of two,
larger circles, colour by adjusted p on the GO-network scale, tighter row pitch. The data source moves
from `percentile_sensitivity_headline.csv` to the grid's
`output/go_grid_headline_by_condition.csv` (§5.8), so this edit **must follow** Phase B.

```python
# before
fig, ax = plt.subplots(figsize=(8.4, 0.42 * len(headline) + 1.4))
...
    for j, arm in enumerate(["5pct", "10pct"]):
        ...
        ax.scatter(j, i, s=260,
                   color="#195f90" if survives else "white",
                   edgecolor="black", linewidth=0.8, zorder=3)
```

```python
# after
# Six conditions, not two: 3 windows x 2 percentiles. Larger circles, tighter rows, and
# the GO-network colour scale so a reader moving between this panel and Figures 4A/5A/6A
# reads one colourbar rather than two conventions.
CONDITIONS = [(w, p) for w in ("5kb", "10kb", "20kb") for p in (5, 10)]
fig, ax = plt.subplots(figsize=(10.6, 0.30 * len(headline) + 1.2))
norm = rl.go_fdr_norm(all_fdr_values_shown)
cmap = plt.get_cmap(rl.GO_FDR_CMAP)
...
    for j, (window, pct) in enumerate(CONDITIONS):
        ...
        ax.scatter(j, i, s=700,
                   color=cmap(norm(-np.log10(fdr))) if survives else "white",
                   edgecolor="black", linewidth=0.8, zorder=3)
```

plus: x tick labels `5 kb / 5 %` … `20 kb / 10 %` with the published condition (10 kb / 5 %) marked;
a horizontal colourbar labelled `-log10(FDR-corrected GO enrichment p-value)`; in-circle text at
`rl.fs(6)`; and the legend reduced to one entry (`open circle = not significant / absent`) since colour
now carries significance.

**5.3f — fonts.** 16 literals through `rl.fs(...)`.

### 5.4 `revision_G3/nb06_go_networks_fdr005.ipynb`

**5.4a — cell 5, compaction ladder.** Replace the `attempts` list with the rungs in §4a, add
`layout_k` / `label_expand` / `label_forces` to each rung's `options`, and record density:

```python
# after (inside write_network, replacing the attempts list)
COMPACT_RUNGS = [
    ("compact k=0.35 x0.70", 0.35, (1.02, 1.05), 0.70),
    ("compact k=0.35 x0.80", 0.35, (1.05, 1.10), 0.80),
    ("compact k=0.45 x0.90", 0.45, (1.10, 1.15), 0.90),
    ("baseline k=0.60 x1.00", 0.60, (1.15, 1.25), 1.00),
    ("relaxed k=0.60 x1.25", 0.60, (1.15, 1.25), 1.25),
    ("relaxed k=0.60 x1.50", 0.60, (1.15, 1.25), 1.50),
]
...
        record["canvas_area_vs_baseline"] = round(factor ** 2, 3)
        record["layout_k"] = k
```

`min_top_n` reduction rungs stay last and unchanged.

**5.4b — cells 10, 12, 13, 14.** No signature change needed — the ladder does the work. S11 keeps its
waived collision check (§8.1) but now runs at the compact rungs, so it gets denser; the waiver stays
recorded in `network_qc.json`.

**5.4c — cell 19, Figure 6B 30 % smaller.**

```python
# before
counts_6b, _ = clustermap_panel(go_families, family_classification, "family_name",
                                "Fig6B", figsize=(14, 9))
# after  (0.7 x 0.7 = 30 % smaller in both directions)
counts_6b, _ = clustermap_panel(go_families, family_classification, "family_name",
                                "Fig6B", figsize=(9.8, 6.3))
```

**5.4d — cell 21, Supplementary 8B 10 % smaller.** `figsize=(14, 9)` → `(12.6, 8.1)`.

**5.4e — cell 27, `go_count_box` half size.** Decision (a): half per side, dot size unchanged,
statistics into the caption.

```python
# before
    fig, ax = plt.subplots(figsize=(4.2, 4))
    ...
    if len(order) == 2:
        ...
        ax.set_title(f"Mann–Whitney U, {p_text} (raw, single test)", fontsize=8.5)
    ax.set_xticklabels([f"{g}\nn={int((data[x_column] == g).sum())}" for g in order],
                       fontsize=8)
```

```python
# after — half the width and height; dot size deliberately unchanged (review item 4).
# At this size a 12 pt title and two-line tick labels do not fit, so the group sizes and
# the test result move into the manuscript legend. Both are already recorded in
# output/figure7_statistics.csv, so nothing leaves the record.
    fig, ax = plt.subplots(figsize=(2.1, 2.0))
    ...
    ax.set_xticklabels(order, fontsize=rl.fs(8))
    # no in-panel title; the caption carries "Mann-Whitney U, raw p = ..., n = a vs b"
```

`sns.stripplot(..., size=3)` stays at 3. The four panels' caption sentences are generated from
`figure7_statistics.csv` and printed by the notebook so they can be pasted into the legend edit.

**5.4f — cell 35, Figure 7H ribbon filter and the new S8C.**

```python
# before
RIBBON_LABEL_THRESHOLD = 5  # the caption's ">= 5 GO terms per ribbon" filter
...
def plot_sankey_classes_families(panel="Fig7H"):
```

```python
# after
# The published Figure 7 caption says "Connecting ribbons were filtered by at least 5 GO
# terms. This filtering was applied to the visualization only." The first version applied
# the threshold to the ribbon LABELS only, so every ribbon was still drawn. Filtering the
# ribbons themselves is what the caption describes; because the filter is visual, bar
# heights and the underlying counts stay at their full values, so retained ribbons do not
# fill their bars. Supplementary Figure S8C is the same plot with no filter.
RIBBON_LABEL_THRESHOLD = 5


def plot_sankey_classes_families(panel="Fig7H", min_ribbon_count=5):
    ...
    dropped = {"class_group": 0, "group_family": 0, "terms": 0}
    ...
    for _, row in counts_by_group_classes.sort_values(...).iterrows():
        ...
        if row["count"] < min_ribbon_count:      # advance the offsets, draw nothing
            dropped["class_group"] += 1
            dropped["terms"] += int(row["count"])
            class_offsets[class_name] -= height
            group_offsets_left[group] -= height
            continue
```

The same guard goes into the second (group → family) loop. The footer strip gains
`Ribbons with fewer than N GO terms are not drawn (visualisation only; bar heights are unfiltered).`
and the function prints the dropped counts. Two calls at the end:

```python
plot_sankey_classes_families(panel="Fig7H", min_ribbon_count=5)
plot_sankey_classes_families(panel="S8C_sankey_full", min_ribbon_count=1)
```

Both dropped counts go into a new `output/sankey_ribbon_filter.json` so the caption can state exactly
what is hidden.

**5.4g — cells 22, 36, 37 deleted** (the four `S_extra_subfamily*` panels). Cell 20's markdown keeps
its explanation of why S8B is the combined class + family map, but the sentence introducing the
subfamily companion panel goes. Delete the four SVG files from `svg/`. The `subfamily_classification`
load in cell 16 becomes unused — remove it too, and keep `go_subfamilies` only where cell 26 still
uses it for `enrichment_df_subfamilies["GO_terms_count"]`.

**5.4h — cell 41** SVG listing loses the `S_extra` prefix and gains `S8C`.

**5.4i — fonts.** 18 literals through `rl.fs(...)`, including `annot_kws={"size": rl.fs(9)}` in
`clustermap_panel` and `font = rl.GLOBAL_FONT_SIZE` in the Sankey (already correct once the constant
changes).

### 5.5 `nb01_permutation_convergence.ipynb`, `nb02_ifna_domain.ipynb`

Font-only: 7 literals each through `rl.fs(...)`, then re-execute. No content change. Their SVGs
(`S14A–C`, `Fig8B`, `Fig8C`, `Fig8C_inset`) are re-emitted purely for the 16 px text.

### 5.6 New — `revision_G3/07a_build_gene_tables.py`

Builds the per-window equivalent of `TEs_on_genes.csv` for all three windows, which is what the GO
grid needs and what does not exist today.

```
Inputs   ../T2T_repeat_masker_processed_sorted.bed        (3,709,429 elements)
         output/windows_{5,10,20}kb.bed                    (38,704 windows each, from 05a)
Output   output/TEs_on_genes_{5,10,20}kb.csv               (38,704 rows each)
```

Construction, deliberately mirroring the published table's semantics:

1. `bedtools intersect -a <repeats> -b windows_W.bed -wa -wb` keeping the **window** chrom/start/end
   as the key, so two TSS windows of the same gene stay distinct rows (the published table is
   per-window, not per-gene — this is the multiple-TSS property flagged in the Limitations).
2. Aggregate per window: `{LINE,LTR,SINE,DNA,Retroposon,RC}_number`, `TE_number`,
   `Divergence_Avg_{class}`, `Average_Divergence_Score`, and the comma-joined
   `TE_subfamilies` / `TE_families` / `TE_classes` / `Divergence_scores` strings in the same format
   `add_family_counts()` parses.
3. **Left join onto the full window list**, so windows with no intersecting TE are retained with
   zeros. This is essential: the GO background is `df["Gene_name"].unique()`, and dropping empty
   windows would shrink the background from 28,738 genes and make FDRs incomparable across windows.

**The regression gate that makes the 5 kb / 20 kb tables trustworthy** — rebuilding 10 kb must
reproduce the published table exactly. Assert all of:

| Quantity | Expected |
|---|---|
| rows | 38,704 |
| unique `Gene_name` | 28,738 |
| windows with `TE_number == 0` | 343 |
| `TE_number` total | 582,540 |
| per class | LINE 169,930 · LTR 51,103 · SINE 302,480 · DNA 57,684 · Retroposon 1,170 · RC 173 |
| `Average_Divergence_Score` | equal to `TEs_on_genes.csv` within 1e-9 on every window |

The script exits non-zero if any of these fails, and prints the 5 kb / 20 kb totals (expected 293,652
and 1,157,235, from `05a`'s log) as the cross-check that the other two windows used the same code path.

### 5.7 New — `revision_G3/07b_go_grid.py`

Runs the 3 × 2 grid for the three GO levels **by calling the existing builders**, so no third copy of
the gene-set construction is created:

```python
go_rerun   = importlib.import_module("06_go_rerun_fdr005")          # 5 % arm
percentile = importlib.import_module("05c_percentile_sensitivity")   # 10 % arm

BUILDERS = {
    (5,  "classes_count"):      go_rerun.run_level_classes_count,
    (5,  "classes_divergence"): go_rerun.run_level_classes_divergence,
    (5,  "families"):           go_rerun.run_level_families,
    (10, "classes_count"):      percentile.build_classes_count,
    (10, "classes_divergence"): percentile.build_classes_divergence,
    (10, "families"):           percentile.build_families,
}
```

```
Outputs  output/GO_grid/GO_<level>_<window>_p<pct>_fdr005.csv          (6 x 3 = 18 files)
         output/GO_grid/GO_<level>_<window>_p<pct>_fdr01_reference.csv (retrieval cut)
         output/GO_grid/INDEX.csv        one row per cell: level, window, percentile,
                                         n_foreground, n_terms_005, n_terms_01, path
         output/GO_grid/MANIFEST.json    script version, builders used, run time, gene-table md5s
```

Rules:

- **The 10 kb cells are not recomputed.** `--reuse-10kb` (default on) copies the existing
  `GO_<level>_fdr005.csv` and `GO_<level>_p10_fdr005.csv` into the grid naming and asserts the copied
  content is row-identical to what a re-run produces on one spot-checked level. This keeps the
  published arm byte-stable and answers §2.3's discoverability half.
- `CLASS_TOP_N = 1436` (5 %) and `2872` (10 %) are counts of genes, and the gene count is 28,738 at
  every window, so the cut sizes do not drift between windows. Assert it.
- `SVA` and `Helitron` remain percentile-invariant by construction at every window (they are
  "any gene carrying such an element", not a ranked cut) — carried through and flagged, not re-run.
- Runtime estimate: 4 new cells × ~62 studies ≈ 250 GO studies ≈ **35–45 min** with the DAG/GAF
  process cache. Backgroundable.

### 5.8 New — `revision_G3/nb07_go_grid_robustness.ipynb`

The comparison Daniil asked for: **10 kb / 5 % against all five alternatives**, for all three GO
levels, with the same measures the existing robustness section uses so the two read as one analysis.

| Section | Content | Output |
|---|---|---|
| 1 | Term counts per cell, per level — a 3 × 2 grid heatmap per level | `S12A_go_grid_term_counts.svg` |
| 2 | Jaccard and **fraction of published terms preserved** for every cell vs 10 kb / 5 %, per level | `S12B_go_grid_preservation.svg` |
| 3 | Every term gained or lost, named, per cell | `output/go_grid_terms.csv` |
| 4 | Headline claims × 6 conditions — the table whose window cells currently read `not evaluated`. **Feeds S13D**, which becomes the 6-condition panel (§5.3e) rather than a second figure of the same thing | `output/go_grid_headline_by_condition.csv` |
| 5 | Group-level survival: which classes/families keep ≥ 1 term in each cell | `S12C_group_survival.svg` |
| 6 | A concordance statement: Spearman of per-group term counts between each cell and the published one, with a label-shuffling permutation test, matching the method already used for the OR concordance | `output/go_grid_concordance.csv` |

Three SVG panels (`S12A`–`S12C`), because the headline-claim grid lives in S13D. Panels reuse
`rl.fs()`, `rl.CLASS_PALETTE`, `rl.CLASS_PALETTE_DIVERGENCE` and `rl.go_fdr_norm`, and every SVG goes
through `save_svg_collision_checked`.

**This closes the honest gap recorded in `project_overview.md` §5**: once section 4 exists, the
`not evaluated` cells in `robustness_headline_claims_by_condition.csv` can be filled.

### 5.9 New — `revision_G3/08_build_supplementary.py` and `revision_G3/supplementary/`

The missing deliverable, assembled by one script so it is reproducible rather than hand-collected.
**Five thematic Excel workbooks**, so a reader opens five files instead of fourteen, and each workbook
is one subject with one sheet per table:

```
revision_G3/supplementary/
  File_S1_TE_TSS_map_and_enrichment.xlsx
      README                    what each sheet is, and which script wrote it
      TSS_TE_intersections      38,704 TSS windows x 25 columns   (was File S1)
      enrichment_classes        6 classes, raw AND adjusted Fisher p (was TableS_class_enrichment_full)
      enrichment_families       44 families, Fisher + permutation columns   (was File S2, corrected)
      enrichment_subfamilies    1,143 subfamilies, same columns             (resolves the S2 mismatch)
  File_S2_gene_sets.xlsx
      README
      by_TE_group               long format: group, gene, rank              (was File S3, 8 sheets)
      by_divergence             long format: group, divergence_tail, gene   (was File S5, 10 sheets)
      by_family                 long format: family, gene, rank             (was File S7, 44 sheets)
  File_S3_gene_ontology.xlsx
      README
      GO_TE_groups              FDR 0.05 + functional-group classification  (was File S4)
      GO_by_divergence          FDR 0.05 + classification                   (was File S6)
      GO_by_family              FDR 0.05 + classification                   (was File S8)
  File_S4_IFNA_domain_and_prior_work.xlsx
      README
      IFNA_elements             all 175 TEs in the 220 kb domain
      IFNA_tests                the four tests: observed, null mean/SD, empirical p
      IFNA_subfamily_composition
      assembly_bound            newly resolved sequence, by class and chromosome
      prior_work_overlap_matrix category x group, Fisher p, Jaccard
      prior_work_categories     mouse->human ortholog mapping
      prior_work_shared_genes
  File_S5_sensitivity_and_robustness.xlsx
      README
      window_classes / window_families / window_concordance / window_flips
      percentile_summary / percentile_terms
      GO_grid_index / GO_grid_preservation / GO_grid_terms / GO_grid_concordance
      headline_by_condition     9 claims x 6 conditions
      geneset_stability / rank_stability
  CHECKSUMS.sha256
  README.md                     the five workbooks, their sheets, and provenance
  figures/PLACEHOLDERS.md       the 14 expected Figure_S1.pdf ... Figure_S14.pdf exports
```

Five things to get right:

1. **Every GO sheet is re-emitted at FDR 0.05.** The baseline copies in
   `manuscript_figures_supplementary_old/` are at 0.1 and are superseded by D2. Assert
   `FDR.max() < 0.05` on each GO sheet before writing.
2. **The old File S2 caption/content mismatch is resolved by carrying both.** The manuscript says
   *"enrichment statistics of TE subfamilies"*, but `Supplementary File 2.csv` contains the **44
   families** and lacks the permutation columns, while Data availability claims S2 carries raw and
   adjusted p-values. Workbook S1 now has `enrichment_families` **and** `enrichment_subfamilies`, both
   with the full column set, so both the caption and the availability claim become true without
   dropping anything.
3. **`by_family` must be long format, not 44 sheets.** The point of the regrouping is fewer things to
   open; 44 sheets inside one workbook would defeat it. Same for `by_TE_group` and `by_divergence`.
4. **Sheet-name limits.** Excel caps sheet names at 31 characters and forbids `[]:*?/\`. Assert both
   before writing, and assert every sheet is under Excel's 1,048,576-row limit — `TSS_TE_intersections`
   is the only one close to being large (38,704 rows, fine) but its comma-joined columns are wide, so
   check the written file opens and report its size.
5. **Figures stay placeholders.** The current figures are exported from Figma by Daniil (Phase 4b), so
   the script writes `figures/PLACEHOLDERS.md` listing the 14 expected filenames rather than inventing
   content.

### 5.10 `revision_G3/13_manuscript_tracked_edits.py`

**5.10a — the working file changes, and with it the script's contract (§2.5).**

```python
# before
BASELINE = os.path.join(DOCX_DIR, "T2T_genes_article_G3_submitted_baseline_260418.docx")
TARGET = os.path.join(DOCX_DIR, "T2T_genes_article_G3_revision_260803.docx")
...
        shutil.copyfile(BASELINE, TARGET)          # line 301 — the idempotency mechanism
    document = docx.Document(TARGET if not DRY_RUN else BASELINE)
```

```python
# after
BASELINE = os.path.join(DOCX_DIR, "T2T_genes_article_G3_submitted_baseline_260418.docx")
# Daniil edited the revision by hand and accepted part of the tracked diff (594 of 1,124
# deletions), so the file below — not the baseline — is the current state of the paper.
# Rebuilding from the baseline would silently discard that work, so IN_PLACE is the default
# and the copyfile path is only reachable with --from-baseline on an explicit new filename.
WORKING = os.path.join(DOCX_DIR, "T2T_genes_article_G3_revision_260804_manual.docx")
TARGET = WORKING
IN_PLACE = True
...
    if not IN_PLACE:
        shutil.copyfile(BASELINE, TARGET)
    document = docx.Document(TARGET)
```

Because the script no longer starts from a pristine copy, its per-edit failure mode must change:
today a missing target string is a hard error (exit 1) that proves no edit was silently skipped. In
place, an edit that is **already present** is the normal case on a re-run. So each edit gains a
three-way outcome — `applied` / `already present, skipped` / **`not found, ERROR`** — and the run
report in `output/manuscript_edit_report.json` records which, with the exit code still non-zero for
the third. That keeps the "a silent skip is impossible" guarantee without making re-runs destructive.

**5.10b — the new tracked edits**, all as `<w:sdt>`-safe replacements through `tracked_replace_safe`:

| Edit | Location | Change |
|---|---|---|
| M1 | Supplementary material section | replace the `File S1`–`File S8` block with five workbook descriptions, each listing its sheets |
| M2 | Supplementary material section | **write the six missing captions, `Figure S9`–`Figure S14`** (§2.6). S9–S11 full networks, **S12 the GO grid**, S13 window/percentile concordance, S14 permutation convergence |
| M3 | Figure S8 caption | panel (C) sentence stays; add that S8C shows **all** ribbons while Figure 7H is filtered at ≥ 5, naming the dropped-ribbon count from `sankey_ribbon_filter.json` |
| M4 | Results, sensitivity paragraph | replace "All comparisons are given in Figure S13 and the accompanying tables." with named items (Figures S12, S13; File S5) and add two sentences reporting the GO grid result |
| M5 | Methods, sensitivity paragraph | state that GO was repeated at both percentiles **at each window size**, i.e. six conditions, replacing the current wording which implies percentiles only |
| M6 | Editor note on the Table 1/2 split | replace the bare filename `TableS_class_enrichment_full.csv` with `File S1, sheet enrichment_classes` |
| M7 | Data availability | rewrite the "supplementary items are as follows" sentence for the five workbooks; keep the raw-p sentence, now naming `File S1` and `File S3` |
| M8 | **27 `File Sn` citations** throughout Results, Methods and the figure legends | renumber to the new workbooks and name the sheet, e.g. `File S3` → `File S3, sheet GO_by_family`. Do this from a written old→new map, in one pass |
| M9 | Figure 7 caption, panels B/E/F/G | add the group sizes and the raw Mann–Whitney p to each sentence, since decision (a) moves them out of the artwork (§5.4e) |
| M10 | Figure 6A caption | if `network_qc.json` reports a `top_n` other than 9 after the font + compaction pass, update the quoted term count. Read it off the file; do not guess |
| M11 | Figure 2 caption | **no change** — §2.1 shows the caption is already right and the panels are being corrected to match it |
| M12 | Figure 7 caption, panel H | **no change** to the ribbon sentence — it already describes the corrected behaviour |

The old→new File S map for M8:

| Old | New |
|---|---|
| File S1 (TSS/TE coordinates) | File S1, sheet `TSS_TE_intersections` |
| File S2 (enrichment) | File S1, sheets `enrichment_families` / `enrichment_subfamilies` |
| File S3 (gene sets by TE group) | File S2, sheet `by_TE_group` |
| File S4 (GO by TE group) | File S3, sheet `GO_TE_groups` |
| File S5 (gene sets by divergence) | File S2, sheet `by_divergence` |
| File S6 (GO by divergence) | File S3, sheet `GO_by_divergence` |
| File S7 (gene sets by family) | File S2, sheet `by_family` |
| File S8 (GO by family) | File S3, sheet `GO_by_family` |
| — (was unnamed: Lu overlap matrix) | File S4, sheet `prior_work_overlap_matrix` |
| — (was unnamed: "the accompanying tables") | File S5 |

M4's new sentences must quote numbers from `output/go_grid_concordance.csv`, via
`11_results_numbers.py` (§5.12), never typed by hand.

### 5.11 `revision_G3/15_house_style.py`

- `TARGET` moves to `T2T_genes_article_G3_revision_260804_manual.docx` so it edits the same file
  `13_…` now edits in place.
- The G2 renaming sweep covers `File S1`–`File S5` and `Figure S9`–`Figure S14`, so the assertion that
  zero `Supplementary Figure n` / `Supplementary File n` strings remain still holds over the new text.
- The existing "already applied" sentinel (the running-title check) is unchanged and still correct.

### 5.12 `revision_G3/11_results_numbers.py`

Add a `go_grid` block re-deriving: term counts per cell per level, preservation fraction per cell,
the six-condition headline table, and the concordance statistics — so every number M4 introduces is
re-derivable and any regeneration surfaces here rather than in proof. Add a `supplementary` block
listing each workbook's sheets and row counts, so the Supplementary material section can be checked
against what the package actually contains.

### 5.13 Documentation

| File | Change |
|---|---|
| `svg/PLACEMENT.md` | correct the Figure 2 D/E/F mapping (§2.1); drop the four `S_extra` rows; add `S8C_sankey_full.svg` → frame `861:32` panel C; add the three `S12*` SVGs as a new frame; state the 16 px font target and that every SVG was re-exported for it; record the achieved compaction per network panel |
| `../G3_figure_pvalue_labels_260803.md` | the 2D/2E/2F row: "per-subfamily empirical p" → "per-family empirical p", label string → `family significant (FDR-corrected empirical p < 0.05)`, and the D/E order corrected; add rows for `S8C` and `S12A–C`; move the Figure 7B/E/F/G statistics from "in-panel annotation" to "caption"; disambiguate `8C` (IFNA) from `S8C` (Sankey) everywhere |
| `CLAUDE.md` (this dir) | the caption-conformance rule as a hard rule ("match panels to the published caption, and assert it"); the font contract (`FONT_SCALE`, pt→px); the inverted collision ladder; `07a`/`07b`/`nb07`/`08` in the run order and layout; `supplementary/` and `output/GO_grid/` in the layout; the corrected Figure 2 D/E/F facts; drop the `S_extra` mention; **replace the manuscript section's idempotency and Reject-All claims per §2.5** and name the working file |
| `project_overview.md` | §5: add the GO-grid result and close the `not evaluated` gap note; §6.2: the `13_…` contract change and the weakened Reject-All guarantee; §6.4: the Figure 2 and Figure 7H caption-conformance findings and the missing-captions finding; §8: move the two closed items out of "cheap improvements" |
| `README.md` (this dir) | run order gains `07a`, `07b`, `nb07`, `08`; §4's idempotency paragraph rewritten per §2.5; §6 layout gains `supplementary/` and `output/GO_grid/` |
| `../REPRODUCE.md` | §4 gains `07a`/`07b`, §5 gains `nb07`, a new §8b for the supplementary package; §8 names the working docx; the near-limit table gains `File_S1_TE_TSS_map_and_enrichment.xlsx` |
| `../.gitignore` | `revision_G3/output/TEs_on_genes_*.csv` (3 × ~25 MB, regenerable by `07a`). The five workbooks are tracked — the largest is well under the 100 MB limit |

---

## 6. Files that do NOT change

| File | Why |
|---|---|
| The four frozen files (`TEs_mapped_on_TSS_analysis.ipynb`, `Gene_ontology_analysis.ipynb`, `download_and_process_files_UCSC_genes.ipynb`, `GO_subfamilies.py`) | frozen; verify the MD5s again at the end |
| `06_go_rerun_fdr005.py`, `05c_percentile_sensitivity.py` | **called, not edited.** `07b` passes them a different `df`; their own outputs stay byte-stable so the published arm is untouched |
| `02_ifna_domain_test.py`, `04a`, `04b`, `05a`, `05b`, `10_tables.py` | no result changes; only their figures get the font pass |
| `01*`, `12*` (permutations, track hub) | unaffected |
| `T2T_genes_article_G3_submitted_baseline_260418.docx` | still read-only, still the reference for the tracked diff. It is no longer the *input* to `13_…` (§2.5), but it is never written |
| `T2T_genes_article_G3_revision_260803.docx` | superseded by the manual file; keep it as the record of what the scripts produced before Daniil's pass. Do not edit or delete |
| `14_build_extensive_discussion.py` and `Extensive_discussion_260803.docx` | built from the pristine baseline, already complete and verified; unaffected by the working-file change |
| `enrichment_*.csv`, `Table1.csv`, `Table2.csv` | D1 keeps N = 500; Figure 2's move to family level **reads** the existing family table, it does not recompute it |
| `../GO_tables/` (FDR 0.1) | still the companion paper's, still not touched (caveat C3) |
| `output/legacy_fdr01_n500/` | the C3 snapshot stays as it is |

---

## 7. Side effects and caveats

- **S1 — the supplementary figure inventory is completed, not extended.** S12 was reserved for a Lu
  et al. figure that was never produced and is cited nowhere, so the GO grid takes S12 and the final
  inventory is a contiguous S1–S14. Six captions (S9–S14) are written for the first time. The
  numbering sweep therefore touches manuscript text, captions, `PLACEMENT.md`, the p-value labels doc
  and the Figma frame names — **by node ID only** (caveat C9).
- **S2 — `Figure 8C` vs `Figure S8C`.** Two different panels one letter apart. Every mention must carry
  the `S`, including in the response letter.
- **S3 — re-executing all five notebooks re-runs the collision ladder.** Achieved `top_n` and canvas
  size may move for 4A, 5A, 6A, S9, S10. Nothing is wrong if they do, but two captions quote those
  numbers (§4b), so read the final values from `network_qc.json` before editing captions.
- **S4 — Figure 6B at 70 % with 1.2× fonts may not be legible.** It carries 13 family rows × 21
  functional-group columns plus stars. If tick labels collide, the fallbacks in order are: shorten the
  functional-group labels, drop the row dendrogram (`row_cluster=False` loses information — avoid),
  or accept 80 % instead of 70 %. Report which was used; do not silently shrink the font back.
- **S5 — Figure 7B/E/F/G lose their in-panel statistics.** Decision (a) buys the half-size panels at
  the cost of moving `n=` and the Mann–Whitney p into the legend (M9). Nothing leaves the record —
  both are in `output/figure7_statistics.csv` — but the legend edit is now load-bearing: if M9 is
  skipped, those four panels report a test with no p-value anywhere.
- **S6 — the supplementary package changes what reviewers saw, in two ways.** The GO sheets move from
  FDR 0.1 to 0.05 (required by D2), and fourteen candidate files become five workbooks, which
  renumbers every `File Sn` citation. Both are improvements, but the response letter must say so
  explicitly, otherwise a reviewer comparing old and new supplementary material finds silent
  differences and a numbering that no longer matches their notes.
- **S7 — repository size.** The five workbooks add roughly 25–30 MB; `TEs_on_genes_counts_subfamilies.csv`
  remains at 85 MB. Nothing crosses GitHub's 100 MB per-file limit, and the `.gz` decision in
  `REPRODUCE.md` §9 is still open.
- **S8 — the GO grid does not re-run permutations.** Enrichment ORs at 5/20 kb come from the existing
  `05b` outputs; the grid is GO only. So a "window effect" in the grid is a gene-set effect, not a
  background effect — state that in M4 so the two analyses are not conflated.
- **S9 — cross-manuscript coupling is unchanged.** Removing the subfamily panels (item 6) removes them
  from *this* paper only; they are reproducible from the companion paper's own pipeline, and
  `GO_subfamilies_fdr005.csv` stays on disk. Note that workbook S1 still ships
  `enrichment_subfamilies` — that is the *enrichment* table the old File S2 caption promised, not the
  subfamily GO analysis, so it does not reintroduce companion-paper content.
- **S10 — `Reject All` no longer restores the baseline** (§2.5). Accepted deletions are permanent. The
  clean-plus-tracked submission pair is produced from the working file, and the diff a reviewer sees
  covers the edits made after Daniil's acceptance pass, not the whole revision. If the journal wants a
  full tracked diff against the submitted version, produce it as a Word *Compare* of the baseline
  against the final clean file rather than from the revision marks.
- **S11 — compute and wall clock.** `07a` ~5 min; `07b` ~35–45 min (backgroundable); `nb07` ~5 min;
  five notebook re-executions ~20–30 min total, dominated by `nb06`'s ladder; `08` ~3 min. No new
  permutations, no new downloads, ~100 MB of new intermediates.

---

## 8. Decisions — resolved

| # | Question | Daniil's answer | Effect on the plan |
|---|---|---|---|
| D-a | Figure 7B/E/F/G: how "two times smaller" resolves against 1.2× fonts | **(a) half per side** | `figsize=(2.1, 2.0)`, dot size unchanged, `n=` and the Mann–Whitney p move into the caption (M9). §5.4e, caveat S5 |
| D-b | Does S13D stay 2 columns or become 6? | **the full 3 × 2 grid** | S13D becomes the 6-condition panel fed by `go_grid_headline_by_condition.csv`; `nb07` emits three SVGs (S12A–C), not four. §5.3e, §5.8 |
| D-c | File S numbering | **regroup thematically** | five workbooks, all 27 `File Sn` citations renumbered with a sheet name (M8). §5.9, caveat S6 |
| D-d | Supplementary tables: how many files | **no more than 5 Excel workbooks, one sheet per table** | §5.9 in full; long-format gene-set sheets replace the 8/10/44-sheet files |
| D-e | Which manuscript file to edit | **`…_revision_260804_manual.docx`** | `13_…` becomes in-place with per-edit skip; `15_…` follows; the Reject-All guarantee weakens. §2.5, §5.10, §5.11, caveat S10 |

### 8.1 The one item left open

**Supplementary Figure S11's waived collision check.** It is the only panel where the check is waived.
With the compaction rungs it can be made denser but will still overlap at 30 terms per group over 44
families. **Proceeding on the recommendation: keep the waiver, denser**, recorded in
`network_qc.json`. Say so if you would rather reduce `top_n` to roughly 20–24 and lose some of the
structure the supplementary version exists to show.

---

## 9. Verification commands

```bash
source ~/venvs/Retroelements_3_11/bin/activate
cd /home/jovyan/Projects/Retroelements/T2T_genes_article/T2T_transposons_genes
WORKING=revision_G3/Revised_manuscript/T2T_genes_article_G3_revision_260804_manual.docx

# --- the frozen four are still untouched -------------------------------------
md5sum TEs_mapped_on_TSS_analysis.ipynb Gene_ontology_analysis.ipynb \
       download_and_process_files_UCSC_genes.ipynb GO_subfamilies.py

# --- 07a reproduces the published 10 kb gene table exactly -------------------
python revision_G3/07a_build_gene_tables.py --verify-10kb
# expect: 38,704 rows / 28,738 genes / 343 empty windows / 582,540 elements
#         LINE 169,930 LTR 51,103 SINE 302,480 DNA 57,684 Retroposon 1,170 RC 173

# --- the grid is complete and nothing above 0.05 leaked ----------------------
python revision_G3/07b_go_grid.py --summary
ls revision_G3/output/GO_grid/GO_*_fdr005.csv | wc -l          # expect 18
python - <<'EOF'
import glob, pandas as pd
for f in sorted(glob.glob("revision_G3/output/GO_grid/GO_*_fdr005.csv")):
    d = pd.read_csv(f, low_memory=False)
    assert d.FDR.max() < 0.05, f
print("all grid cells filtered at FDR < 0.05")
EOF
python revision_G3/07b_go_grid.py --check-reuse       # expect "row-identical" x 6

# --- Figures 2D-2F reproduce the three published sentences ------------------
jupyter nbconvert --execute --inplace revision_G3/nb03_relabelled_figures.ipynb
# the notebook's own assertion prints:
#   Figures 2D-2F reproduce the three published significance statements

# --- fonts: every SVG carries 12 pt (= 16 px in Figma) as its base ----------
# NOTE (corrected 2026-08-05): do NOT grep for "16px". Matplotlib writes font sizes in SVG
# USER UNITS, and its root element is width="...pt" with a matching viewBox, so the user
# unit IS a point: a 12 pt font is written as `font-size: 12px` and Figma renders it at
# 12 x 96/72 = 16.0 px. The check that means something is that every size present is
# exactly 1.2x one of the original point values.
python - <<'EOF'
import glob, re
expected = {round(o * 1.2, 2) for o in (6, 6.5, 7, 7.5, 8, 8.5, 9, 9.6, 10, 16)}
seen = set()
for f in glob.glob("revision_G3/svg/*.svg"):
    seen.update(float(v) for v in re.findall(r"font-size:\s*([0-9.]+)px", open(f).read()))
print("sizes:", sorted(seen))
assert all(round(v, 2) in expected for v in seen), sorted(
    v for v in seen if round(v, 2) not in expected)
print("every size is exactly 1.2x an original point size")
EOF
# measured: [7.2, 7.8, 8.4, 9.0, 9.6, 10.2, 10.8, 12.0], base 12 pt -> 16.0 px in Figma

# --- packing actually got denser --------------------------------------------
python -c "
import json; qc=json.load(open('revision_G3/output/network_qc.json'))
for p,r in qc.items():
    print(f\"{p:28s} area x{r.get('canvas_area_vs_baseline','?')} top_n={r.get('top_n')} \"
          f\"collisions={r.get('label_collisions','waived')}\")"
# expect Fig6A, S9, S10, S11 at <= 0.70 area, collisions 0 (S11 waived)

# --- the subfamily extras are gone ------------------------------------------
ls revision_G3/svg/ | grep S_extra          # expect no output
grep -rn "S_extra" revision_G3/*.ipynb revision_G3/svg/PLACEMENT.md  # expect no hits

# --- the Sankey filter is real and quantified -------------------------------
python -c "import json;print(json.load(open('revision_G3/output/sankey_ribbon_filter.json')))"
# expect nonzero dropped ribbon counts for Fig7H and zero for S8C_sankey_full

# --- the supplementary package: five workbooks, sheets as documented --------
python revision_G3/08_build_supplementary.py --verify
ls revision_G3/supplementary/File_S*.xlsx                  # expect exactly 5
sha256sum -c revision_G3/supplementary/CHECKSUMS.sha256
python - <<'EOF'
import glob, pandas as pd
for f in sorted(glob.glob("revision_G3/supplementary/File_S*.xlsx")):
    x = pd.ExcelFile(f)
    print(f.split('/')[-1], len(x.sheet_names), "sheets:", x.sheet_names)
    assert all(len(s) <= 31 for s in x.sheet_names), f
EOF

# --- the manuscript names every workbook, sheet and figure ------------------
pandoc -t plain --wrap=none "$WORKING" -o /tmp/ms.txt
for n in 1 2 3 4 5; do grep -q "File S$n\." /tmp/ms.txt || echo "File S$n not described"; done
for n in $(seq 1 14); do grep -q "^Figure S$n\." /tmp/ms.txt || echo "Figure S$n caption missing"; done
grep -oP "File S\d+" /tmp/ms.txt | sort -u        # expect only S1..S5, no S6+
grep -c "Figure S12" /tmp/ms.txt                  # expect >= 2 (citation + caption)
grep -c "accompanying tables" /tmp/ms.txt         # expect 0 (M4 replaced it)

# --- manuscript integrity after the in-place tracked edits ------------------
~/venvs/collagen_3_11/bin/python revision_G3/13_manuscript_tracked_edits.py
~/venvs/collagen_3_11/bin/python revision_G3/15_house_style.py
python -c "
import json; r=json.load(open('revision_G3/output/manuscript_edit_report.json'))
print({k: v for k, v in r.items() if 'not_found' in str(v).lower()} or 'no missing targets')"
~/venvs/collagen_3_11/bin/python - <<'EOF'
import zipfile, re
p = ("revision_G3/Revised_manuscript/"
     "T2T_genes_article_G3_revision_260804_manual.docx")
z = zipfile.ZipFile(p); doc = z.read("word/document.xml").decode("utf8", "ignore")
print("sdt:", len(re.findall(r"<w:sdt>", doc)),          # expect 129
      "bibliography:", doc.count("MENDELEY_BIBLIOGRAPHY"),
      "webextension:", any(n.startswith("word/webextensions/") for n in z.namelist()),
      "ins:", len(re.findall(r"<w:ins ", doc)),          # expect > 380
      "del:", len(re.findall(r"<w:del ", doc)))          # expect >= 530
EOF
python revision_G3/11_results_numbers.py                 # every quoted number re-derived
```

---

## 10. TODO

### Phase A — foundations (do these first; everything else depends on them)

- [x] `revision_lib.py` — add `FONT_SCALE = 1.2`, `BASE_FONT_SIZE`, `GLOBAL_FONT_SIZE = 12.0`, `fs()`
- [x] `revision_lib.py` — add `shade()` and `WINDOW_SHADES`
- [x] `revision_lib.py` — add `GO_FDR_CMAP` and `go_fdr_norm()`; use it in both network functions
- [x] `revision_lib.py` — move `CLASS_PALETTE_DIVERGENCE` in from `nb06` cell 3 (keep the provenance comment)
- [x] `revision_lib.py` — add `layout_k`, `layout_iterations`, `label_expand`, `label_forces` to `save_go_network_svg` and `save_go_network_svg_families_by_classes`; thread to `spring_layout` and `adjust_text`
- [x] `revision_lib.py` — record `canvas_area_vs_baseline` and `layout_k` in the `qc` dict

### Phase B — the missing analysis (long-running; launch early)

- [x] `07a_build_gene_tables.py` — write it; `--verify-10kb` must pass all six invariants of §5.6
- [x] Run `07a` for 5 kb, 10 kb, 20 kb
- [x] `07b_go_grid.py` — write it; reuse the two 10 kb cells, run the four new ones
- [x] Run `07b` in the background; confirm 18 `fdr005` files and `INDEX.csv`
- [x] `07b --check-reuse` — the reused 10 kb cells are row-identical to the published tables (plus a full re-run spot check of classes_count: 425/425 (group, term) pairs identical)

### Phase C — figure code changes (after Phase A, before re-execution)

- [x] `nb03` cell 3 — derive `Status` on the families table
- [x] `nb03` cell 7 — switch to `enrichment_df_families`; **swap D/E to the caption's order**; family wording in the legend; singleton classes as points; `n/a` for untestable pairs
- [x] `nb03` cell 17 — add the three caption-conformance assertions
- [x] `nb03` — all 38 font literals through `rl.fs()` (31 measured, not 38; the plan's count was an estimate)
- [x] `nb05` cell 1 — load `family_to_class`
- [x] `nb05` cell 3 — S13A: `width_ratios [1, 3.6]`, class dots `s=110`, colour by class × window shade, two-block legend
- [x] `nb05` cell 5 — S13B: colour by class + legend
- [x] `nb05` cell 9 — S13C: left row-annotation strip via `row_annotation_colour()`
- [x] `nb05` — all 16 font literals through `rl.fs()` (20 rl.fs() calls: 8 pre-existing literals + 12 written into the new panel code)
- [x] `nb06` cell 5 — replace the ladder with the compaction rungs; record density
- [x] `nb06` cell 19 — Fig6B `figsize=(9.8, 6.3)`
- [x] `nb06` cell 21 — S8B `figsize=(12.6, 8.1)`
- [x] `nb06` cell 27 — `go_count_box` `figsize=(2.1, 2.0)`, dot size unchanged, in-panel title and `n=` labels removed, caption sentences printed for M9
- [x] `nb06` cell 35 — `min_ribbon_count`; skip sub-threshold ribbons in both loops; footer sentence; write `sankey_ribbon_filter.json`; two calls (`Fig7H` at 5, `S8C_sankey_full` at 1)
- [x] `nb06` — delete cells 22, 36, 37 and the now-unused `subfamily_classification` load; update cell 20 markdown and cell 41's prefix list
- [x] `nb06` — all 18 font literals through `rl.fs()` (12 measured, not 18)
- [x] `nb01`, `nb02` — font literals through `rl.fs()`
- [x] Delete `svg/S_extra_subfamily_functional_groups_fdr005.svg`, `S_extra_subfamily_go_by_class_fdr005.svg`, `S_extra_subfamily_go_by_length_fdr005.svg`, `S_extra_subfamily_go_by_significance_fdr005.svg`

### Phase D — re-execute, then tune

- [x] Execute `nb01`, `nb02`, `nb03`, `nb05`, `nb06` — zero error outputs, `nb03`'s assertions pass (plus nb07; run order is nb01/02/03/06 -> nb07 -> nb05)
- [x] Inspect `network_qc.json`: is every network panel at ≤ 0.70 baseline area with zero collisions? If not, tune `layout_k` / `label_expand` and re-run that panel only — **no: the target is not reachable.** Measured and documented instead: 1.2x fonts inflate every label box 1.44x in area, so 30 % less area asks for the same labels in 0.49x the space they already needed. Only S11 reaches 0.49 (check waived). Raising the adjust_text forces was measured to make packing WORSE; more spring iterations changes nothing. The ladder now records `canvas_area_vs_baseline`, `compaction_target_met` and the best compaction rung reached, per panel
- [x] Check Figure 6B and the four half-size Figure 7 panels for legibility; apply the §7 S4 fallbacks if needed and record which — Fig6B at (9.8, 6.3) and the four Fig7 panels at (2.1, 2.0) render without error and with no in-panel statistics; **visual legibility is Daniil's call at Figma placement**, and the §7 S4 fallbacks were not needed to produce them
- [x] Verify the font check: 16 px present in every SVG, no size below 1.2× its old value — the check as written was wrong (matplotlib writes user units, not px); corrected in §9. Every size present is exactly 1.2x an original: 7.2, 7.8, 8.4, 9.0, 9.6, 10.2, 10.8, 12.0, base 12 pt = 16.0 px in Figma
- [x] `nb07_go_grid_robustness.ipynb` — write and execute; three `S12*` SVGs + four output CSVs including `go_grid_headline_by_condition.csv`
- [x] `nb05` cell 11 — S13D as the **6-condition** panel: `s=700`, `coolwarm` on −log10(FDR) + colourbar, row pitch 0.30, reading the grid table; re-execute `nb05`
- [x] Fill the `not evaluated` cells in `robustness_headline_claims_by_condition.csv` from the grid (nb05 is the single writer; nb07 supplies the grid verdicts)

### Phase E — the supplementary package

- [x] `08_build_supplementary.py` — write it; assert FDR < 0.05 on every GO sheet, sheet names ≤ 31 chars, long-format gene-set sheets
- [x] Run it; check five workbooks, their sheet inventories, `CHECKSUMS.sha256`, `README.md`, `figures/PLACEHOLDERS.md`
- [x] Confirm the workbook sizes open cleanly in Excel/LibreOffice (`TSS_TE_intersections` is the wide one) — every workbook round-trips through `pd.ExcelFile`/openpyxl with all sheets readable, and `08 --verify` asserts every sheet is inside Excel's 31-character sheet-name, 32,767-character cell and 1,048,576-row limits. A visual check in Excel/LibreOffice is still Daniil's to make

### Phase F — documentation and manuscript (last, per caveat C12)

- [x] `13_manuscript_tracked_edits.py` — switch to the working file, in-place mode, three-way per-edit outcome (§5.10a) — verified by `--dry-run`: 49 edits, 2 applied, 47 already present, **0 not found**. Structural inserts (the tables, the Discussion subsections, the URL paragraph) each gained an `already present` marker so a re-run cannot duplicate them, and the search strings now tolerate 15's `Supplementary File 4` -> `File S4` rename
- [x] `13_…` — add edits M1–M10; write the old→new File S map into the script as data, not prose — Stage M. The map is `FILE_S_LABEL_MAP` (one source of truth: Stage M renumbers from it and the already-applied detection reads it back). M11 and M12 confirmed as no-change.
- [x] Run `13_…`; confirm no `not_found` outcomes in `manuscript_edit_report.json` — **23 applied, 49 already present, 0 not found**, and a second run is a byte-identical no-op (document.xml identical with revision ids normalised).
- [x] `15_house_style.py` — retarget the working file; extend the G2 sweep to `File S1`–`S5` and `Figure S9`–`S14`; run it **after** 13 — retargeted and the sweep extended to Figures S1-S14 with panel letters A-E. It correctly **refuses to run** on the working file (running-title sentinel), because the house style is already applied there; that refusal is now right rather than an obstacle, but M1-M10's new text will need a way to re-run the G2 sweep over newly inserted paragraphs only
- [x] Read the final `top_n` for Figure 6A off `network_qc.json` and set the caption (M10) — `top_n = 5`; the caption reads "Up to five terms", generated from the QC file rather than typed.
- [x] `11_results_numbers.py` — add the `go_grid` and `supplementary` blocks; re-run
- [x] `svg/PLACEMENT.md` — Figure 2 mapping corrected; `S_extra` rows dropped; `S8C` and `S12*` added; compaction and font notes
- [x] `../G3_figure_pvalue_labels_260803.md` — per-family 2D/2E/2F with the corrected letters; add `S8C`, `S12A–C`; Figure 7B/E/F/G statistics move to "caption"; disambiguate `8C` / `S8C`
- [x] `CLAUDE.md`, `project_overview.md`, `README.md`, `../REPRODUCE.md`, `../.gitignore` per §5.13 — including the rewritten idempotency and Reject-All statements — all five updated. `.gitignore` also **un-ignores the three April gene-set workbooks**: `08` reads the published gene sets from them, so without that a clean checkout could not rebuild workbook S2
- [x] Response-letter notes: GO sheets re-emitted at FDR 0.05, and the supplementary files regrouped into five workbooks with a numbering map (caveat S6) — written to `revision_G3/response_letter_notes_260805.md`: the two supplementary changes with the old->new numbering map, the completed window sensitivity with its untidy result, the three caption-conformance corrections, and the legibility/packing trade-off. The letter itself is Daniil's to write
- [x] Re-verify the four frozen MD5s and run the full §9 block — the four MD5s match; no `nbimporter`/`import_ipynb` anywhere; 07a's gate, 07b's reuse check (incl. a full re-run spot check), the 18-cell FDR check, the S_extra removal, the Sankey filter record, and `08 --verify` all pass. The §9 font check was **wrong as written** and is corrected in place. The manuscript greps in §9 depend on M1-M10, which are not written yet


---

## 11. State at the end of the 2026-08-05 implementation pass

46 of 50 TODO items are complete. **Everything computational, every figure, the supplementary
package and all six documentation files are done and verified.** What remains is one coherent
block of authorial work on the manuscript text, deliberately not attempted in this pass.

### What remained after that pass — all of it done on 2026-08-05, see §12 below

- [x] **`13_…` edits M1–M10** — six new figure captions (S9–S14), 27 `File Sn` citations to
  renumber against the map in §5.10b, and several paragraph rewrites (M1, M4, M5, M7).
- [x] **Run `13_…` for real** (it has only been dry-run) and confirm no `not_found` outcomes.
- [x] **M10: Figure 6A's caption** — the achieved value is **`top_n = 5`**, read off
  `output/network_qc.json`. Was 9.
- [x] **Re-running 15's G2 sweep over the newly inserted text.** `15_…` is retargeted and its
  sweep extended to Figures S1–S14, but it refuses to run twice and the house style is already
  applied to the working file. M1–M10's new paragraphs will introduce fresh
  `Supplementary Figure`/`File` strings that need the insertion-aware pass only.

The reason for stopping *at the time* was that this block is writing, not mechanics: it is 30-plus
authored sentences going into a document that carries Daniil's own edits and his acceptance of
594 tracked deletions, where a botched tracked edit is not cheaply reversible. The mechanics it
depends on are finished and proven — `13_…` is in-place, idempotent, guarded against duplicating
its structural inserts, and dry-runs at **49 edits / 2 applied / 47 already present / 0 not
found**.

### The two prerequisites are satisfied

The two things §4b said had to be settled before any caption could be written are now measured
rather than assumed:

| Panel | terms per group | canvas vs nominal |
|---|---|---|
| 4A | 10 | × 1.56 |
| 5A | 10 | × 2.25 |
| **6A** | **5** (was 9) | × 2.25 |
| S9 | 30 | × 2.25 |
| S10 | 29 | × 2.25 |
| S11 | 30, collision check waived | **× 0.49** |

And the Figure 7 legend sentences M9 needs are generated by `nb06`, not to be retyped:

```
(B) Mann-Whitney U, raw p = 0.113 (single test); n = 11 (Significant) vs 2 (Non-Significant).
(E) Mann-Whitney U, raw p = 0.208 (single test); n = 13 (Yes) vs 31 (No).
(F) Mann-Whitney U, raw p = 0.029 (single test); n = 13 (Yes) vs 31 (No).
(G) Mann-Whitney U, p < 0.001 (single test); n = 13 (Yes) vs 31 (No).
```

### Decisions this pass surfaced for Daniil

1. **Figure 6A loses four GO terms per family (9 → 5)** as the price of 16 px text. The panel is
   clean and the caption can say five, but if nine terms matter more than the font size, the
   alternative is to exempt this one panel from `FONT_SCALE`.
2. **"30–50 % denser packing" is not achievable** for the network panels — the geometry is in
   §4b and in `PLACEMENT.md` §4b. Only S11 reaches it, because its collision check is waived.
3. **`hAT-Charlie / MHC class I protein complex` survives 1 of the 6 grid conditions.** It was
   already flagged as percentile-fragile; the grid strengthens the case for softening or
   dropping it.
4. **The GO window effect is real and not tidy** — preservation of published terms drops to
   0.440 at worst across windows, against 0.85–0.93 across percentiles. M4 has to report that
   honestly rather than describe the results as window-robust.


---

## 12. Stage M — completed 2026-08-05

All of M1–M10 are applied to `Revised_manuscript/T2T_genes_article_G3_revision_260804_manual.docx`
as tracked changes. Final state: **129 `<w:sdt>` (128 citations + 1 bibliography), 395 `<w:ins>`,
547 `<w:del>`, zero `<w:t>` inside a `<w:del>`**, webextension part intact. Reject All restores
all eight original `Supplementary File n` captions with their descriptive text and reverts the
corrected cross-references to the published form.

| Edit | What was done |
|---|---|
| M1 | The eight `File Sn` captions tracked-deleted, replaced by five workbook captions listing their sheets |
| M2 | Six captions written for the first time: **Figures S9–S14**, quoting the achieved term counts (S9 30, S10 29, S11 30) |
| M3 | Figure S8's panel C says every ribbon is drawn and names what 7H hides (36 + 50 ribbons, 146 terms). Also **two wrong cross-references corrected**: the published text pointed at `Supplementary Figure 5C`/`5B` for the Sankey and the clustermap, but S5 has only panels A and B and shows TSS distributions — they are panels C and B of S8 |
| M4 | Results sensitivity paragraph reports the GO grid, quoting 0.44 / 0.68 preservation, ρ = 0.61, p ≤ 0.022, 3 of 9 claims, and states the grid is a gene-set effect not a background effect |
| M5 | Methods states six conditions, the shared construction, the constant 28,738 background, and that no permutations were repeated |
| M6 | The Table 1/2 editor note now cites `File S1, sheet enrichment_classes` instead of a bare filename |
| M7 | Data availability rewritten for five workbooks; the raw-p sentence now names File S1 and File S3 |
| M8 | **All 27 `File Sn` citations resolved**: 11 body citations renumbered with their sheet, 8 removed with the old captions (M1), 8 rewritten in Data availability (M7) |
| M9 | Figure 7's legend carries the group sizes and raw Mann–Whitney p for panels B, E, F and G, generated from `figure7_statistics.csv` |
| M10 | Figure 6A's caption reads "Up to **five** terms", read off `network_qc.json` |
| M11, M12 | No change, as the plan predicted — Figure 2's caption was already correct and Figure 7H's ribbon sentence already described the corrected behaviour |

### Five defects found and fixed while making Stage M idempotent

Each was caught by running 13 twice and diffing, which is now part of the procedure:

1. **A replacement whose result contains its own search string looped 200 times.** Figure 7's
   caption keeps its sentence and appends to it, and `replace_inside_insertions` restarted its
   search from the beginning of the run, so it appended until the loop bound. It now makes one
   left-to-right pass, advancing past each replacement.
2. **`present_in_any_form` silently always returned False when handed a paragraph.**
   `find_all_p` iterates only *direct children*, so scoping it to a paragraph found no
   paragraphs at all — every already-applied check was a no-op. It now reads the element's own
   text.
3. **The label map cascaded.** `File S4` → `File S3, sheet GO_TE_groups` was followed by the
   `File S3` rule rewriting that output into `File S2, sheet by_TE_group, sheet GO_TE_groups`.
   Replaced by a single regex pass with a callback.
4. **`File S1` → `File S1, sheet …` renumbered its own output** on a second run. Fixed by a
   per-pair "result already present" skip.
5. **M1, the Lu et al. paragraph, and the repository-URL edit had no already-applied marker**, so
   a second run duplicated them. A tracked-deleted paragraph is still in the XML and `text_of`
   still reads its `delText`, so M1's guard has to check for its *own output* first rather than
   for the absence of the old captions.

### What `15_house_style.py` needs (nothing)

It still refuses to run — the house style is already in this file — and that is now the correct
outcome rather than an obstacle: Stage M's new text was written in G3 house style throughout, so
the G2 assertion holds over it with zero legacy `Supplementary Figure n` / `Supplementary File n`
strings in the document. The new captions were built by cloning the surrounding captions' run
properties, so their `Figure Sn.` labels carry the same `Strong` character style as S1–S8.

### Still Daniil's

Refreshing Mendeley in both docx files, the Figma placement and export, the Zenodo DOI, publishing
the track hub, the ten browser-verify items, and the response letter (notes in
`response_letter_notes_260805.md`). One judgement call is worth a look: **M4's paragraph now
contains both the earlier percentile-only sentence ("of the nine abstract-level claims tested at
both percentiles, eight survive both") and the new six-condition result ("3 of the 9 claims tested
survive all six conditions")**. Both are accurate at their own scope and the new sentence says
"all six conditions" explicitly, but you may prefer to merge them.
