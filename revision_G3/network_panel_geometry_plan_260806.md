# Implementation plan — reshaping the full GO network panels S9, S10, S11

**Date:** 2026-08-06 · **Author:** Claude Code, for Daniil's approval · **Status:** awaiting approval,
nothing implemented

---

## Overview

Three supplementary network panels have the wrong page geometry. `S9_full.svg` and `S10_full.svg`
are square (1630.7 × 1608.5 pt) and must become **half as wide at the same height**;
`S11_full.svg` must become **1.5× taller at the same width**, and its GO-term circles and labels
must be distributed evenly instead of collapsing into tight clumps joined by an impenetrable mat
of edges. In all three, text labels must not intersect.

The key design decision is that **this is not a canvas-scaling problem and cannot be solved by the
existing compaction ladder.** The ladder in `nb06` cell 5 multiplies *both* figure dimensions by one
factor (`CANVAS_RUNGS`), so it cannot express "narrower at the same height" at all. And halving the
width without changing anything else is arithmetically hopeless: labels are drawn as **single lines
of up to 107 characters** (nothing in this repository wraps label text — `grep -rn textwrap` returns
nothing), and one 107-character label at 16 px is ≈ 642 pt, i.e. **79 % of the target width**.

So the plan makes three changes of substance, in this order of importance:

1. **Wrap label text** to a bounded line width. This is what the reference frames do — every long
   term name in `Figure 5 old` and `Figure 6 old` is broken over 2–4 short lines — and it is the
   only lever that buys width without discarding GO terms.
2. **Pin the canvas** per panel to the exact requested pt geometry, and *assert* it after writing,
   so "half the width, same height" is a checked property rather than an intention.
3. **Stop manufacturing the clumps in S11.** They are not an accident of `spring_layout`: line 807
   of `revision_lib.py` deliberately adds same-class family edges at `weight=10.0` to the layout
   graph, ten times the maximum Jaccard weight, which collapses each class onto a point.

Nothing about the *content* of the three panels changes: FDR 0.05, Jaccard ≥ 0.1, ≥ 0 shared genes,
terms ≤ 1,000 genes, and `top_n = 30` are all quoted verbatim in the Figure S9–S11 captions and stay
exactly as they are. Only geometry, label typesetting, edge rendering and layout spacing change —
plus, if the collision check cannot be satisfied at half width, the achieved terms-per-group, which
is the one content-visible fallback and is called out as a decision below.

---

## Background / reference data

### Measured current state and the targets

| Panel | current SVG (pt) | aspect | target (pt) | target aspect | change |
|---|---|---|---|---|---|
| `S9_full.svg` | 1630.70 × 1608.48 | 1.01 | **815.35 × 1608.48** | 0.507 | width ÷ 2, height fixed |
| `S10_full.svg` | 1630.70 × 1608.48 | 1.01 | **815.35 × 1608.48** | 0.507 | width ÷ 2, height fixed |
| `S11_full.svg` | 803.92 × 845.28 | 0.95 | **803.92 × 1267.92** | 0.634 | height × 1.5, width fixed |

With `bbox_inches="tight"` the written size is the trimmed bbox, not `figsize`, so these targets are
reachable exactly only by fixing the canvas (see change 3.1). At 72 pt/in the pinned figures are
**11.32 × 22.34 in** (S9, S10) and **11.17 × 17.61 in** (S11).

### The reference frames

Located by node ID, not name (caveat C9 — the canvas is shared with the subfamilies paper):

| Frame | Node ID | Size | Aspect |
|---|---|---|---|
| `Figure 5 old` | **`861:28`** | 1222 × 1749 | 0.70 |
| `Figure 6 old` | **`861:33`** | 1238 × 1956 | 0.63 |

Both are portrait because each stacks the network (panel A) over the clustermap (panel B). **S11's
target aspect of 0.634 is Figure 6 old's aspect of 0.633** — the requested change lands the panel
exactly on the reference column format, which is good evidence the target is the right one.

What the reference frames do that the current SVGs do not, read off the two screenshots:

| Reference behaviour | Current code |
|---|---|
| long term names wrapped over 2–4 short lines | single line, up to 107 characters |
| labels sit immediately beside their node | `adjust_text(..., max_move=None)` — a label may travel arbitrarily far |
| colourbar and legends tucked into a corner | `plt.colorbar(shrink=0.5, pad=0.05)` reserves a right-hand strip; the families version reserves 15 % of the width via `tight_layout(rect=[0, 0, 0.85, 1])` |
| thin, pale connection lines | Jaccard edges at `width = weight * 10` (up to 10 pt) and class edges at `Fold Enrichment / 2` |
| GO circles spread evenly around their group node | `spring_layout(k=0.35)` on a graph with `weight=10.0` same-class edges |

### Label geometry, measured from the three SVGs

| Panel | text nodes | labels > 3 ch | mean length | longest | longest as % of target width |
|---|---|---|---|---|---|
| S9 | 85 | 75 | 25.3 ch | **75 ch** (≈ 450 pt) | **55 %** |
| S10 | 83 | 75 | 26.0 ch | **107 ch** (≈ 642 pt) | **79 %** |
| S11 | 100 | 87 | 22.2 ch | 52 ch (≈ 312 pt) | 39 % |

Wrapped at 22 characters, the S10 worst case becomes 5 lines of ≤ 22 ch ≈ 132 pt, **16 %** of the
target width. Total label ink then needs roughly 75 labels × ~3 lines × 16 px ≈ 2,700 pt of stacked
height against 1,608 pt of canvas, i.e. about two label columns — which 815 pt of width can hold at
~6 label widths per row. Tight, but not geometrically excluded. This is the basis for expecting the
half-width target to be reachable; it is not a guarantee, hence the ranked fallbacks in §7.

### Current QC state (`output/network_qc.json`)

| Panel | canvas vs nominal | terms/group | collisions | rung reached |
|---|---|---|---|---|
| S9 | × 2.25 (25.5 × 22.5 in) | 30 | 0 | short labels, ×1.50, expand 1.02 |
| S10 | × 2.25 | 29 | 0 | short labels, ×1.50, expand 1.15, top_n = 29 |
| S11 | × 0.49 (13.3 × 11.9 in) | 30 | **waived** | compact rung ×0.70, k = 0.35 |

### Baseline md5s — the three panels that must NOT change

`nb06` writes six SVGs; re-running it regenerates all six. These three are already approved and
placed, and must come back byte-identical:

```
dee24f6bf5aab2ac7b0b5b7b44f84628  svg/Fig4A_simplified.svg
b9c18242ffb9ac1b9aac009fee56c69e  svg/Fig5A_simplified.svg
2382b91c04852ed936a252760d5f9cd0  svg/Fig6A_simplified.svg
99de0393697454bd98e252fa7ee13085  svg/S9_full.svg      <- will change
a6c1f4e015389e64ebce30c3abaaf2cf  svg/S10_full.svg     <- will change
d72afc2215b4db05ec9c4a8ad5b38cb4  svg/S11_full.svg     <- will change
```

Every new argument in §1 and §2 therefore **defaults to the current behaviour**, so the three
main-text panels are untouched by construction and the md5 check is a real regression gate.

---

## Files to change

### 1. `revision_G3/revision_lib.py` — `save_go_network_svg` (S9, S10; shared with 4A/5A)

#### 1a. Import `textwrap` (top of file, with the stdlib imports)

```python
# before
import ast
import os
```
```python
# after
import ast
import os
import textwrap
```

#### 1b. New keyword arguments — signature, ~line 524

```python
# before
    layout_k=0.6,
    layout_iterations=None,
    label_expand=(1.15, 1.25),
    label_forces=((0.6, 0.8), (0.4, 0.6)),
    baseline_figsize=None,
):
```
```python
# after
    layout_k=0.6,
    layout_iterations=None,
    label_expand=(1.15, 1.25),
    label_forces=((0.6, 0.8), (0.4, 0.6)),
    baseline_figsize=None,
    label_wrap=None,          # wrap label text at N characters (None = single line, as published)
    label_max_move=None,      # cap on how far adjust_text may carry a label from its node
    show_colorbar=True,       # False -> use svg/Fig456A_colorbar_vector.svg in Figma instead (G12)
    show_size_legend=True,
    edge_alpha=0.4,
    edge_width_cap=None,      # pt ceiling on drawn edge width; presentation only
    tight_bbox=True,          # False -> figsize maps 1:1 to the written SVG size
):
```

Docstring gains a paragraph: the four new rendering arguments are presentation-only and change no
statistic; `label_wrap` is the width lever the reference frames use; `tight_bbox=False` is what makes
a pinned canvas exact.

#### 1c. Wrap the label text — ~line 675

```python
# before
    texts = [
        ax.text(x, y, G.nodes[node]["label"], fontsize=font_size, fontweight="normal")
        for node, (x, y) in pos.items()
    ]
```
```python
# after
    def _label(node):
        text = G.nodes[node]["label"]
        return textwrap.fill(text, label_wrap) if label_wrap else text

    texts = [
        ax.text(x, y, _label(node), fontsize=font_size, fontweight="normal",
                ha="center", va="center", linespacing=0.95)
        for node, (x, y) in pos.items()
    ]
```

`ha/va="center"` matters once labels are multi-line: left-baseline anchoring puts a 4-line label
entirely below-right of its node, which is the opposite of "labels closer to their cognate term".
`linespacing=0.95` keeps a wrapped label compact vertically.

#### 1d. Keep labels near their node — ~line 679

```python
# before
    adjust_text(
        texts,
        ax=ax,
        arrowprops=dict(arrowstyle="->", color="gray", lw=0.5),
        force_text=label_forces[0],
        force_static=label_forces[1],
        expand=label_expand,
        max_move=None,
    )
```
```python
# after
    adjust_text(
        texts,
        ax=ax,
        arrowprops=dict(arrowstyle="-", color="gray", lw=0.4, shrinkA=1, shrinkB=1),
        force_text=label_forces[0],
        force_static=label_forces[1],
        expand=label_expand,
        max_move=label_max_move,
    )
```

`adjustText` 1.3.0 accepts `max_move` (verified by signature inspection), and `None` is the current
value, so passing the new argument through changes nothing by default. The arrow style becomes a
plain leader line rather than an arrowhead, matching the reference.

#### 1e. Make the colourbar and size legend optional — ~lines 654–673

Wrap the `plt.colorbar(...)` block in `if show_colorbar:` and the `for count in [10, 100, 1000]` /
`ax.legend(...)` block in `if show_size_legend:`. Both default `True`, so 4A/5A are unaffected.

#### 1f. Cap edge width and expose alpha — ~line 644

```python
# before
    nx.draw_networkx_edges(
        G, pos, edgelist=edges,
        width=[d["weight"] / 2 for _u, _v, d in edges],
        edge_color=[d.get("color", "#E0E0E0") for _u, _v, d in edges],
        alpha=0.4, ax=ax,
    )
```
```python
# after
    def _edge_width(d):
        width = d["weight"] / 2
        return min(width, edge_width_cap) if edge_width_cap else width

    nx.draw_networkx_edges(
        G, pos, edgelist=edges,
        width=[_edge_width(d) for _u, _v, d in edges],
        edge_color=[d.get("color", "#E0E0E0") for _u, _v, d in edges],
        alpha=edge_alpha, ax=ax,
    )
```

#### 1g. Honour `tight_bbox` on both save paths — ~lines 694–701

`bbox_inches="tight" if tight_bbox else None` in the `fig.savefig` call, and the same in
`save_svg_collision_checked` (change 4b).

### 2. `revision_G3/revision_lib.py` — `save_go_network_svg_families_by_classes` (S11; shared with 6A)

#### 2a. Same seven new arguments as 1b, plus two layout arguments

```python
# after (added to the signature)
    label_wrap=None,
    label_max_move=None,
    show_colorbar=True,
    show_size_legend=True,
    show_jaccard_legend=True,
    edge_alpha=0.3,
    edge_width_cap=None,
    tight_bbox=True,
    same_class_weight=10.0,      # was hard-coded; lower values spread families of one class
    legend_rect=(0, 0, 0.85, 1), # tight_layout rect; (0,0,1,1) reclaims the legend strip
):
```

#### 2b. The clumping cause — ~line 807

```python
# before
            if G.nodes[fams[i]]["class_name"] == G.nodes[fams[j]]["class_name"]:
                H.add_edge(fams[i], fams[j], weight=10.0)
```
```python
# after
            if G.nodes[fams[i]]["class_name"] == G.nodes[fams[j]]["class_name"]:
                H.add_edge(fams[i], fams[j], weight=same_class_weight)
```

This single literal is the reason the GO circles "form very compact groups". Jaccard edge weights are
≤ 1 and family→term weights are fold enrichments; a same-class weight of 10.0 dominates both, so
`spring_layout` pulls every family of a class onto one point and drags its 30 terms with it. S11 will
use **2.0**, which keeps classes recognisably grouped — the point of the figure — without collapsing
them. The default stays 10.0 so Figure 6A is unchanged.

#### 2c. Wrapped, centred labels — ~lines 858–873

Same `textwrap.fill` treatment as 1c, preserving the existing per-node `color`/`fontweight` logic
(family labels bold black, term labels coloured by dominant class).

#### 2d. `max_move`, edge cap, alpha, optional furniture, `tight_bbox`

As 1d–1g. The edge width function here needs both branches capped:

```python
# before
    def _width(d):
        return d["weight"] / 2 if d["edge_type"] == "fam_go" else d["weight"] * 10
```
```python
# after
    def _width(d):
        width = d["weight"] / 2 if d["edge_type"] == "fam_go" else d["weight"] * 10
        return min(width, edge_width_cap) if edge_width_cap else width
```

#### 2e. The Jaccard legend must use the same transform — ~line 904

The legend draws `lw=v * 10` for J = 0.1, 0.2, 0.5. If drawn edges are capped, the legend
becomes a lie unless it is capped identically:

```python
# before
    j_lines = [
        Line2D([0], [0], color="gray", lw=v * 10, alpha=0.4, label=f"J={v:.2f}")
        for v in [0.1, 0.2, 0.5]
    ]
```
```python
# after
    j_lines = [
        Line2D([0], [0], color="gray",
               lw=min(v * 10, edge_width_cap) if edge_width_cap else v * 10,
               alpha=edge_alpha + 0.1, label=f"J={v:.2f}")
        for v in [0.1, 0.2, 0.5]
    ]
```

This is the Figure 7H lesson applied early: a presentation-only filter that the legend does not know
about is a defect, not a saving.

### 3. `revision_G3/nb06_go_networks_fdr005.ipynb` — cell 5, `write_network`

#### 3a. New parameters and a pinned-canvas attempt list

```python
# before
def write_network(table, panel, plotter, figsize, settings, extra=None,
                  min_top_n=MIN_TOP_N_MAIN_TEXT, require_clean=True, palette=None):
```
```python
# after
# Wrap widths tried in order, narrowest last: wrapping is the width lever the reference
# frames (Figma 861:28 / 861:33) use, and it is tried BEFORE dropping any GO term.
WRAP_RUNGS = (None, 26, 22, 18)

def write_network(table, panel, plotter, figsize, settings, extra=None,
                  min_top_n=MIN_TOP_N_MAIN_TEXT, require_clean=True, palette=None,
                  pinned_figsize=None, wrap_rungs=(None,), target_pt=None, render=None):
```

`render` carries the presentation-only kwargs (`show_colorbar`, `edge_width_cap`, `label_max_move`,
…) straight through to the plotter. When `pinned_figsize` is given the canvas rungs are **not**
climbed — the attempt list becomes wrap × label kind × expand × descending `top_n`, all at the
pinned size — and `tight_bbox=False` is forced so the written SVG is exactly `figsize × 72` pt.

#### 3b. Record the new facts in QC

Each successful attempt adds `label_wrap`, `pinned` and the measured `svg_pt` / `aspect` to the
panel's `network_qc` record, so PLACEMENT.md and the captions can be updated from measurement
rather than memory — the same discipline as `top_n`.

### 4. `revision_G3/revision_lib.py` — new geometry assertion

#### 4a. `assert_svg_geometry(path, target_pt, tol=0.03, qc=None)`

New function beside `assert_no_label_collisions`: parses `width="…pt"` / `height="…pt"` from the
written SVG, raises `RuntimeError` if either is outside `tol` of its target, and records
`svg_pt` + `aspect` into `qc`. Rationale: the whole request is a geometry request, so it gets the
same treatment as C16 gives label overlap — a panel that silently comes out the wrong shape is the
failure mode to design out.

#### 4b. `save_svg_collision_checked(..., tight_bbox=True)`

Thread the flag through so the checked path and the unchecked path write identically.

### 5. `revision_G3/nb06_go_networks_fdr005.ipynb` — cells 12, 13, 14

```python
# before — cell 12
write_network(go_classes_count, "S9_full", rl.save_go_network_svg, figsize=(17, 15),
              settings=FULL, min_top_n=20, require_clean=False)
```
```python
# after — cell 12
# Half the published width at the same height (815.35 x 1608.48 pt), portrait like the
# reference frames. The colourbar leaves the artwork: at this width it would eat ~10 % of
# it, and Fig456A_colorbar_vector.svg is already the vector replacement (G12).
write_network(go_classes_count, "S9_full", rl.save_go_network_svg,
              figsize=(11.32, 22.34), pinned_figsize=(11.32, 22.34),
              target_pt=(815.35, 1608.48), wrap_rungs=WRAP_RUNGS,
              settings=FULL, min_top_n=20, require_clean=True,
              render=dict(show_colorbar=False, show_size_legend=True,
                          edge_width_cap=2.5, edge_alpha=0.25,
                          label_max_move=(40, 40)))
```

Cell 13 (`S10_full`) is the same with `palette=CLASS_PALETTE_DIVERGENCE`. Cell 14:

```python
# after — cell 14
# Same width, 1.5x the height (803.92 x 1267.92 pt) -> aspect 0.634, which is the aspect of
# reference frame 861:33 ("Figure 6 old"). same_class_weight 10.0 -> 2.0 and k 0.35 -> 0.9
# are what stop the families of a class collapsing onto one point.
write_network(go_families, "S11_full", rl.save_go_network_svg_families_by_classes,
              figsize=(11.17, 17.61), pinned_figsize=(11.17, 17.61),
              target_pt=(803.92, 1267.92), wrap_rungs=WRAP_RUNGS,
              settings=FULL, min_top_n=20, require_clean=False,
              extra={"family_to_class": family_to_class},
              render=dict(same_class_weight=2.0, layout_k=0.9, layout_iterations=300,
                          show_colorbar=False, legend_rect=(0, 0, 1, 1),
                          edge_width_cap=2.0, edge_alpha=0.22,
                          label_max_move=(30, 30)))
```

`require_clean=True` for S9/S10 because "ensure the labels do not intersect" is explicit; S11 keeps
`require_clean=False` so the documented waiver remains the fallback, but the ladder now has wrapping
to try first and the collision count at the waiver must be **reported before and after** so the
improvement is measured rather than claimed.

### 6. `revision_G3/nb06_go_networks_fdr005.ipynb` — markdown cells 4 and 11

Cell 4 holds the ladder narrative and the measured packing table; cell 11 explains the S9–S11
settings and the waiver. Both must state the new pinned-canvas rung, the wrap rungs, and the new
achieved numbers, or the notebook's own prose contradicts its output.

### 7. `revision_G3/svg/PLACEMENT.md`

- §3 note 3: `4A and 5A keep 10, S9 at 30, S10 at 29` → the achieved values after the re-run.
- §4b table: new `canvas vs nominal`, terms/group, collisions and rung for S9, S10, S11, plus a row
  or footnote giving the new pt geometry and the reference-frame aspect it matches.
- §4b prose: the "packing target is not achievable" paragraph stays true for 4A/5A/6A but must now
  record that **label wrapping**, which the first ladder never tried, is what made a 2× width
  reduction possible for S9/S10 — otherwise the file reads as if compaction had been exhausted.
- §3 note 4 (the S11 waiver) → updated to the new collision count, or removed if zero is reached.

### 8. `revision_G3/CLAUDE.md`

Lines 295–303: `4A and 5A keep 10, S9 keeps 30, S10 29`, the S11 waiver sentence, and the
"30–50 % denser packing is not achievable" bullet. The last one needs the wrapping caveat added for
the same reason as §7.

### 9. `revision_G3/project_overview.md`

Line ~565: *"`label_collisions == 0` for every panel except the recorded S11 waiver"* — update if the
waiver disappears; and the §7 verification block gains the geometry assertion.

### 10. `revision_G3/13_manuscript_tracked_edits.py` — the S11 caption

`S9_TOP_N` / `S10_TOP_N` / `S11_TOP_N` come from `qc_top_n(...)`, so the caption numbers track
`network_qc.json` automatically **on a fresh insert**. Two things do not track automatically:

```python
# before (~line 2066, inside NEW_FIGURE_CAPTIONS)
     f"the simplified view. Up to {S11_TOP_N} terms per family are shown, with the same edge "
     f"and term-size filters as Figure S9, and node colour denotes the FDR-corrected "
     f"enrichment p-value. At this term count the label field is saturated and some labels "
     f"overlap; Figure 6A is the legible view of the same network and this panel is provided "
     f"for the full structure. GO terms are counted at FDR < 0.05."),
```
```python
# after — only if S11 reaches zero collisions
     f"the simplified view. Up to {S11_TOP_N} terms per family are shown, with the same edge "
     f"and term-size filters as Figure S9, and node colour denotes the FDR-corrected "
     f"enrichment p-value. Term names are wrapped and connection lines are drawn at a capped "
     f"width for legibility; edge weights are unchanged. GO terms are counted at FDR < 0.05."),
```

**And the caption is already in the working document.** M2 is guarded by
`structural_edit_needed(body, "M", …, marker)`, which skips the insert when the marker sentence is
already present — so re-running `13_…` will report *"already present"* and leave the **old** numbers
and the old "some labels overlap" sentence in
`Revised_manuscript/T2T_genes_article_G3_revision_260804_manual.docx`. Updating them needs either a
new tracked `edit()` for the changed sentences or a hand edit in Word. This is the single easiest
thing in this plan to miss.

---

## Files that do NOT need to change

| File | Why not |
|---|---|
| The four frozen files (`../TEs_mapped_on_TSS_analysis.ipynb`, `../Gene_ontology_analysis.ipynb`, `../download_and_process_files_UCSC_genes.ipynb`, `../GO_subfamilies.py`) | hard rule 1; the plotters were copied into `revision_lib.py` precisely so they can be changed here |
| `06_go_rerun_fdr005.py`, `07b_go_grid.py`, `GO_tables_fdr005/` | the GO statistics are untouched; this is a rendering change only |
| `nb01`, `nb02`, `nb03`, `nb05`, `nb07` | none of them calls the network plotters; `nb05`'s dependency is on `nb07`, not `nb06` |
| `svg/Fig4A_simplified.svg`, `Fig5A_simplified.svg`, `Fig6A_simplified.svg` | every new argument defaults to today's behaviour, so re-running `nb06` must reproduce them byte-for-byte — verified by md5, not assumed |
| `svg/Fig456A_colorbar_vector.svg` | already the vector colourbar; dropping the baked-in colourbar from S9–S11 is what it exists for (G12) |
| `10_tables.py`, `11_results_numbers.py`, `08_build_supplementary.py` | no network geometry anywhere in them |
| `15_house_style.py` | its caption edits quote 4A/5A/6A term counts and "the full network is Figure S9/S10/S11" — the pointers stay valid and the counts are not S9–S11's |
| `output/go_label_shortnames.csv` | the 47 curated short labels stay as they are; wrapping composes with them rather than replacing them |

---

## Side effects and caveats

1. **Terms may be lost on S9/S10, and that is the one content-visible risk.** With
   `require_clean=True` the ladder will drop `top_n` toward 20 if wrapping alone cannot clear the
   collisions at half width. The achieved value flows into `network_qc.json` → the S9/S10 captions,
   but only on a fresh insert (§10). Ranked remedies before any term is dropped: wrap at 26 → 22 →
   18 characters, then drop the colourbar (already planned), then the size legend, then `top_n`.
   **If the ladder ends below `top_n = 26` I will stop and report rather than ship a quieter figure**
   — losing a fifth of the terms to a width change is a decision for you, not for the ladder.
2. **`nb06` runtime grows.** Every attempt renders the whole figure because collision detection needs
   a real canvas. Pinning the canvas removes 6 canvas rungs but adds 4 wrap rungs, so the attempt
   count per panel is comparable (~16 before top_n reduction); expect roughly the current 10–20 min
   per full panel, worst case longer for S9/S10 now that they must reach zero at half width.
3. **The C16 check is not loosened anywhere.** S9/S10 move from "clean on a big canvas" to "clean on
   a small one"; S11 keeps its recorded waiver as a fallback. No `pad_points` change, no
   `check_collisions=False` outside the existing waiver path.
4. **Edge capping is presentation-only and must be disclosed.** Edge weights, Jaccard values and
   fold enrichments are unchanged; only drawn stroke width is clipped. The Jaccard legend is capped
   identically (§2e) and the caption says so (§10) — the Figure 7H defect, avoided in advance.
5. **The Figma side is manual and one frame already holds the old S11.** Frame
   **`935:11876`** ("S11_full 1", 1072 × 1127, aspect 0.95) contains the current S11 and will need
   re-import plus a resize to aspect 0.634; S9 and S10 have no frames yet (`svg/PLACEMENT.md` §3
   lists them as still to be created). Hard rule 4 stands: nothing in `revision_G3/` writes to
   Figma.
6. **`network_qc.json` is read by `13_manuscript_tracked_edits.py` only** (verified by grep), and its
   `qc_top_n("Fig6A_simplified")` feeds M10. Since 6A must regenerate identically, M10's number is
   unaffected — but re-running `nb06` rewrites the whole QC file, so keep a copy of the current one
   for the before/after comparison.
7. **`tight_bbox=False` leaves a hairline margin** that `bbox_inches="tight"` would have trimmed.
   This is deliberate: exact, assertable geometry is worth more than a few pt of whitespace, and
   Figma frames crop. The alternative — keep tight bbox and calibrate `figsize` in a measure-and-
   rescale loop — is available if you prefer zero margin, at the cost of two renders per attempt.
8. **Aspect change redistributes the layout, so the panels will not merely be narrower — they will
   look different.** `spring_layout` returns a square coordinate box that the axes stretch to the
   figure aspect, so a 0.507-aspect canvas compresses every horizontal distance by ~2×. Node
   *positions* are still deterministic (`seed=42`), so the figure is reproducible, but it is not the
   old figure squeezed: it is a new arrangement of the same graph.
9. **No cache to invalidate**, and no S3/output-count consequences — this touches only `svg/` and
   `output/network_qc.json`.

---

## Verification commands

```bash
source ~/venvs/Retroelements_3_11/bin/activate
cd /home/jovyan/Projects/Retroelements/T2T_genes_article/T2T_transposons_genes

# 0. keep the before state
cp revision_G3/output/network_qc.json /tmp/network_qc_before.json

# 1. re-run the only notebook involved
jupyter nbconvert --execute --inplace revision_G3/nb06_go_networks_fdr005.ipynb

# 2. the three main-text panels must be byte-identical
md5sum -c <<'EOF'
dee24f6bf5aab2ac7b0b5b7b44f84628  revision_G3/svg/Fig4A_simplified.svg
b9c18242ffb9ac1b9aac009fee56c69e  revision_G3/svg/Fig5A_simplified.svg
2382b91c04852ed936a252760d5f9cd0  revision_G3/svg/Fig6A_simplified.svg
EOF

# 3. geometry: exactly half the width / 1.5x the height, within 3 %
python - <<'PY'
import re
targets = {"S9_full": (815.35, 1608.48), "S10_full": (815.35, 1608.48),
           "S11_full": (803.92, 1267.92)}
for panel, (tw, th) in targets.items():
    svg = open(f"revision_G3/svg/{panel}.svg").read()
    w = float(re.search(r'width="([\d.]+)pt"', svg).group(1))
    h = float(re.search(r'height="([\d.]+)pt"', svg).group(1))
    ok = abs(w - tw) / tw < 0.03 and abs(h - th) / th < 0.03
    print(f"{panel:10s} {w:8.2f} x {h:8.2f} pt  target {tw:8.2f} x {th:8.2f}  "
          f"aspect {w/h:.3f}  {'OK' if ok else 'FAIL'}")
PY

# 4. labels: zero collisions on S9 and S10, and S11 no worse than before
python - <<'PY'
import json
before = json.load(open("/tmp/network_qc_before.json"))
after = json.load(open("revision_G3/output/network_qc.json"))
for panel in ["S9_full.svg", "S10_full.svg", "S11_full.svg"]:
    b, a = before[panel], after[panel]
    print(f"{panel:16s} collisions {b['label_collisions']} -> {a['label_collisions']}   "
          f"top_n {b['top_n']} -> {a['top_n']}   wrap {a.get('label_wrap')}   "
          f"size {a.get('svg_pt')}")
assert after["S9_full.svg"]["label_collisions"] == 0
assert after["S10_full.svg"]["label_collisions"] == 0
assert after["Fig6A_simplified.svg"]["top_n"] == 5, "M10's number must not move"
PY

# 5. no term silently lost
python -c "
import json; qc=json.load(open('revision_G3/output/network_qc.json'))
for p in ['S9_full.svg','S10_full.svg','S11_full.svg']:
    print(p, 'top_n =', qc[p]['top_n'], '(was 30 / 29 / 30)')"

# 6. the doc numbers match the QC file
grep -n "S9 at 30\|S10 at 29\|S9 keeps 30" revision_G3/svg/PLACEMENT.md revision_G3/CLAUDE.md
```

Visual check before any Figma work: render each new SVG to PNG at the size it will occupy and look at
it beside the reference frames.

```bash
python -c "
import cairosvg
for p,w in [('S9_full',816),('S10_full',816),('S11_full',804)]:
    cairosvg.svg2png(url=f'revision_G3/svg/{p}.svg', write_to=f'/tmp/{p}.png',
                     output_width=w*2, background_color='white')
print('wrote /tmp/S9_full.png /tmp/S10_full.png /tmp/S11_full.png')"
```

---

## TODO

- [x] `revision_lib.py` — add `import textwrap` (and `import re`, needed by `assert_svg_geometry`)
- [x] `revision_lib.py` — `save_go_network_svg`: add the new kwargs, all defaulting to current
      behaviour (**8, not 7**: `leader_lines` was added so the plain leader line is opt-in — see
      the deviation note below)
- [x] `revision_lib.py` — `save_go_network_svg`: wrap label text, centre-anchor, `linespacing=0.95`
      (anchor and line spacing applied **only when wrapping**, so an unwrapped panel keeps
      matplotlib's defaults)
- [x] `revision_lib.py` — `save_go_network_svg`: pass `max_move=label_max_move`, plain leader lines
      **behind `leader_lines=False`**
- [x] `revision_lib.py` — `save_go_network_svg`: guard colourbar and size legend behind flags
- [x] `revision_lib.py` — `save_go_network_svg`: `_edge_width` cap + `edge_alpha`
- [x] `revision_lib.py` — `save_go_network_svg`: honour `tight_bbox` on both save paths
- [x] `revision_lib.py` — `save_go_network_svg_families_by_classes`: add the new kwargs (10 with
      `leader_lines`)
- [x] `revision_lib.py` — families: replace the hard-coded `weight=10.0` with `same_class_weight`
- [x] `revision_lib.py` — families: wrap labels, keeping per-node colour and bold family names
- [x] `revision_lib.py` — families: `max_move`, edge cap, `edge_alpha`, `legend_rect`, `tight_bbox`
      (+ when `legend_rect` reclaims the strip, both legends move inside the axes — otherwise they
      are drawn off-canvas and silently cropped by the pinned save)
- [x] `revision_lib.py` — families: cap the Jaccard legend line widths identically (§2e)
- [x] `revision_lib.py` — new `assert_svg_geometry(path, target_pt, tol, qc)`
- [x] `revision_lib.py` — `save_svg_collision_checked`: thread `tight_bbox`
- [x] `revision_lib.py` — docstrings for every new argument; note which are presentation-only
- [x] **Not in the plan, added:** `find_offpage_labels` / `assert_labels_on_page` /
      `record_offpage_labels`. A pinned page crops what `bbox_inches="tight"` would have grown to
      include, so without this "half the width" could be achieved by pushing a term name off the
      page. Enforced on every checked pinned write, recorded (not raised) on the waived one.
      Measured: **0 off-page labels in all three panels**
- [x] `nb06` cell 5 — add `WRAP_RUNGS`, `pinned_figsize`, `wrap_rungs`, `target_pt`, `render`
- [x] `nb06` cell 5 — pinned attempt list, no canvas rungs (**wrap × label kind × expand × k ×
      top_n**: `K_RUNGS = (0.6, 0.9, 1.2)` was added because at half width no wrap rung alone is
      clean at the published k = 0.6, and raising k costs no information — see the deviation note)
- [x] `nb06` cell 5 — record `label_wrap`, `pinned`, `svg_pt`, `aspect` in `network_qc`
- [x] `nb06` cell 5 — call `rl.assert_svg_geometry` after each pinned write (including the waived one)
- [x] `nb06` cell 12 — S9 at (11.32, 22.34) in, `require_clean=True`, render kwargs
- [x] `nb06` cell 13 — S10 likewise, keeping `CLASS_PALETTE_DIVERGENCE`
- [x] `nb06` cell 14 — S11 at (11.17, 17.61) in, `same_class_weight=2.0`, k laddered from 0.9,
      300 iterations, `min_top_n = 30`, no legends (deviation notes 3 and 4)
- [x] `nb06` markdown cell 4 — ladder narrative: pinned ladder, wrap rungs, k rungs, why both come
      before `top_n`
- [x] `nb06` markdown cell 11 — S9–S11 settings and the waiver's new status
- [x] Re-run `nb06`; 4A/5A md5s reproduce their **content** (page size + exact label set, via the
      new `svg/compare_panels.py`) but **no panel is byte-reproducible** — `adjust_text` is
      stochastic and matplotlib stamps a `<dc:date>`. 6A came out at 7 terms per family instead of
      5 on both runs, so all three main-text SVGs **and their QC records** were restored from the
      approved run; their md5s are now exactly the plan's three values
- [x] Confirm geometry within 3 % for all three panels — S9/S10 815.04 × 1608.48 pt (aspect 0.507),
      S11 804.24 × 1267.92 pt (aspect 0.634); worst deviation 0.04 %
- [x] Confirm zero collisions on S9 and S10 (**both, at 30 terms per group**); S11 waived at
      **7 overlapping pairs, down from 11**, with all 30 terms per family kept
- [x] `top_n` never fell: S9 **30**, S10 **30** (up from the published 29), S11 **30** — no term lost
- [x] Render all three to PNG and inspect — S9/S10 read as the reference column format; S11's clumps
      are gone and its terms are spread across the canvas
- [x] `svg/PLACEMENT.md` §3 notes 3 and 4, §3 frame table (new aspects + frame `935:11876`),
      §4b table and prose, and a new **§4c** for the pinned canvas
- [x] `revision_G3/CLAUDE.md` — term counts, waiver, packing bullet, the pinned-canvas bullet, the
      non-reproducibility bullet, and the new `revision_lib` functions/arguments
- [x] `revision_G3/project_overview.md` §7 — waiver count, geometry assertion, the
      non-reproducibility note
- [x] `13_manuscript_tracked_edits.py` — S9/S10/S11 captions: the edge-width cap disclosed, wrapping
      stated **from `network_qc.json`** via the new `qc_wrap_sentence()` (S9/S10 turned out not to
      need wrapping, so the caption must not claim it), S11 pointed at Figure 6A for the node-size
      and Jaccard scales its legends no longer carry
- [ ] **Hand-check the working docx**: M2 is already inserted, so the S9–S11 caption numbers and the
      S11 overlap sentence will NOT update on a re-run (caveat §10). **S10's number moved 29 → 30**,
      so this now matters for two panels, not one
- [ ] Figma (Daniil): re-import S11 into frame `935:11876` at aspect 0.634; create the S9/S10 frames
      at aspect 0.507; place `Fig456A_colorbar_vector.svg` on all three (none carries a colourbar now)
```

---

## Deviations from this plan, and why

1. **The md5 gate cannot pass and never could.** Verification step 2 asks for byte-identical
   4A/5A/6A. Measured: matplotlib writes a `<dc:date>` into every SVG, and `adjust_text` places
   labels stochastically — two renders in one process with the same `PYTHONHASHSEED` differ in
   label coordinates. `svg/compare_panels.py` was written to compare what *is* reproducible (page
   size, mark count, exact label set); 4A and 5A pass it. **Figure 6A does not**: it came out at 7
   terms per family rather than the approved 5, on both runs. Since 6A is approved and placed and
   M10 quotes 5, all three main-text panels and their `network_qc.json` records were restored from
   the approved run. The extra two terms on 6A are available if you want them, at the cost of
   re-placing that panel and updating M10.
2. **`leader_lines` is a flag, not a change.** §1d changes the arrow style unconditionally, which
   would alter 4A/5A/6A. Since §"Files that do NOT need to change" requires the opposite, the plain
   leader line became opt-in and only S9–S11 use it.
3. **`K_RUNGS` was added to the pinned ladder.** With wrapping alone, S9 and S10 still collide at
   half width at every wrap width (S9: 4 pairs at wrap 22, 1 pair at wrap 18 even after dropping
   `top_n` to 20). Raising `spring_layout`'s k spreads nodes at constant page area and costs no
   information, so it belongs in the ladder before any term is dropped — and it is what got both
   panels clean **with every term kept**. S9 needed k = 0.9 with no wrapping at all; S10 k = 1.2.
4. **S11: `min_top_n = 30`, and no legends.** The waiver replays the closest attempt the ladder
   reached; with `top_n` free to fall to 20 the closest attempt is always the one with the fewest
   labels, so the first run silently traded ten terms per family for two fewer overlapping pairs —
   in the one panel whose purpose is the complete structure. Pinning `min_top_n` to the full 30
   removes that trade. Separately, `legend_rect = (0, 0, 1, 1)` leaves nowhere to hang a legend:
   both legends landed on node labels, and the C16 check cannot see legend text. They were dropped
   from the artwork exactly as the colourbar was, and the caption now points at Figure 6A, which
   carries the same node-size and Jaccard scales.
5. **The open question was not acted on.** S9 and S10 still draw every label in black; term labels
   are not tinted by dominant class. It is a two-line change if you want it.

---

## Open question I could not resolve from the repository

The reference frames show **wrapped labels and corner-tucked legends**, which is what this plan
adopts. What they do not settle is whether you want the two supplementary panels to also adopt the
reference's *label colouring* (term labels tinted by their dominant class, as `Figure 6 old` does and
as `save_go_network_svg_families_by_classes` already does for S11) in the **class-level** panels S9
and S10, where `save_go_network_svg` currently draws every label black. It is a two-line change and
would make S9/S10 read like the reference, but it is a visual-identity decision rather than a
geometry one, so it is not in the TODO above.
