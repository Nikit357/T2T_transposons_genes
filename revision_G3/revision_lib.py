"""Shared helpers for the G3 revision of manuscript G3-2026-406828.

This module is AUTHORITATIVE for everything the G3 revision touches. The three
analysis notebooks in the parent directory (`TEs_mapped_on_TSS_analysis.ipynb`,
`Gene_ontology_analysis.ipynb`, `download_and_process_files_UCSC_genes.ipynb`)
and `GO_subfamilies.py` are FROZEN for this revision and are never edited. The
helper functions they define in-cell were therefore COPIED here and modified in
the copies (plan §2, §7.2-7.4, caveat C19).

Provenance of every copied function
-----------------------------------
| Function                                  | Copied from                              |
|-------------------------------------------|------------------------------------------|
| run_goatools_enrichment                   | Gene_ontology_analysis.ipynb cell 6      |
| run_goatools_ordered_enrichment           | Gene_ontology_analysis.ipynb cell 6      |
| save_go_network_svg                       | Gene_ontology_analysis.ipynb cell 6      |
| visualize_go_class_network                | Gene_ontology_analysis.ipynb cell 36     |
| save_go_network_svg_families_by_classes   | Gene_ontology_analysis.ipynb cell 175    |
| CLASS_PALETTE (== class_names_p)          | Gene_ontology_analysis.ipynb cell 3      |
| CLASS_PALETTE_DIVERGENCE                  | Gene_ontology_analysis.ipynb cell 111    |

What changed relative to the frozen originals
---------------------------------------------
1. `FDR_THRESHOLD = 0.05` is the default everywhere (decision D2, reviewer minor
   comment 2). The frozen originals stay at 0.1 and keep working; do NOT import
   from them (caveat C19).
2. The three network plotters gained `min_shared_genes`, `max_term_genes` and an
   enforced label-collision check (WP12, caveat C16).
3. The GO DAG and the 190 MB GAF are loaded ONCE per process and cached. The
   frozen originals reload both inside every call, which is why the frozen
   `GO_subfamilies.py` takes hours over 1,143 subfamilies. The revision's GO
   re-run makes ~1,200 calls, so caching is a prerequisite rather than a
   nicety. The GO results are unaffected: the same DAG and the same association
   dict are handed to `GOEnrichmentStudy` either way.
4. `load_permutation_counts()` reads the compacted permutation store produced by
   `01b_compact_permutation_results.py`, replacing the 6.37 GB
   `consolidated_random_data.csv` the frozen notebook read (WP1b).
5. Font sizes are scaled by `FONT_SCALE = 1.2` so imported SVG text lands at
   16 px in Figma rather than 13.33 px; every hard-coded size in the notebooks
   goes through `fs()` (review item 13).
6. The two network plotters expose their layout parameters (`layout_k`,
   `layout_iterations`, `label_expand`, `label_forces`) so the collision fallback
   ladder can try *compaction* before canvas enlargement, and record what they
   were drawn at in `network_qc.json`.

Never update `go-basic.obo` / `goa_human.gaf` — they are frozen at the
December 31 2025 snapshot cited in the Methods (caveat C5).
"""

from __future__ import annotations

import ast
import functools
import glob
import io
import json
import os
import re
import subprocess
import textwrap
from itertools import combinations

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from adjustText import adjust_text
from goatools.anno.gaf_reader import GafReader
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.obo_parser import GODag
from matplotlib.lines import Line2D

# --------------------------------------------------------------------------- #
# Plotting conventions (project CLAUDE.md; identical across both manuscripts)
# --------------------------------------------------------------------------- #

# Matplotlib sizes are in points; SVG consumers render them in pixels at 96 dpi,
# so 1 pt = 96/72 = 1.3333 px. A 10 pt label therefore arrives in Figma as
# 13.33 px, against the 16 px the frames use. FONT_SCALE = 1.2 makes the base
# 12 pt = exactly 16 px. Every hard-coded size in the notebooks goes through
# fs() so the whole figure scales together, not only the base size — otherwise
# 8 pt axis labels would sit at 10.67 px against 16 px titles.
FONT_SCALE = 1.2
BASE_FONT_SIZE = 10
GLOBAL_FONT_SIZE = BASE_FONT_SIZE * FONT_SCALE  # 12.0 pt -> 16.0 px in Figma


def fs(points):
    """Scale a font size expressed at the original 10 pt baseline."""
    return round(points * FONT_SCALE, 2)


plt.rcParams["svg.fonttype"] = "none"  # keeps SVG text editable in Illustrator
plt.rcParams.update({"font.size": GLOBAL_FONT_SIZE})

# Copied verbatim from Gene_ontology_analysis.ipynb cell 3 (`class_names_p`).
CLASS_PALETTE = {
    "LINE": "#cc660b",
    "LTR": "#70453c",
    "SINE": "#ab1f20",
    "DNA": "#195f90",
    "Retroposon": "#765297",
    "RC": "#238023",
    "SVA": "#765297",
    "Helitron": "#238023",
    "TE_top": "green",
    "TE_bottom": "red",
}

# The six real classes, without the divergence-group and alias keys above.
CLASS_PALETTE_SHORT = {
    "LINE": "#cc660b",
    "LTR": "#70453c",
    "SINE": "#ab1f20",
    "DNA": "#195f90",
    "Retroposon": "#765297",
    "RC": "#238023",
}

CLASS_NAMES = ["LINE", "LTR", "SINE", "DNA", "Retroposon", "RC"]

# Figures 5A and 5B do NOT group by class: the frozen notebook (cell 108, executed
# before the 5A network in cell 112, and again in cell 114 before the 5B
# clustermap) rewrites `class_name` as "{class}_{highest|lowest}" so the two
# divergence tails are separate groups. Collapsing them would silently merge the
# highest- and lowest-divergence gene sets, which is the entire contrast Figure 5
# exists to show. Copied verbatim from frozen GO cell 111
# (`class_names_p_modified`). Lives here rather than in nb06 because nb05's S13C
# row annotation needs the same colours (plan §5.1e).
CLASS_PALETTE_DIVERGENCE = {
    "LINE": "#cc660b", "LINE_lowest": "#db9454", "LINE_highest": "#8f4708",
    "LTR": "#70453c", "LTR_lowest": "#9b7d77", "LTR_highest": "#4e302a",
    "SINE": "#ab1f20", "SINE_lowest": "#c46263", "SINE_highest": "#771516",
    "DNA": "#195f90", "DNA_lowest": "#5e8fb1", "DNA_highest": "#114264",
    "Retroposon": "#765297", "Retroposon_lowest": "#9f86b6",
    "Retroposon_highest": "#523969",
    "RC": "#238023", "RC_lowest": "#65a665", "RC_highest": "#185918",
    "TE_all": "grey", "TE_all_lowest": "#a6a6a6", "TE_all_highest": "#595959",
}


def shade(colour, amount):
    """Lighten (amount < 0) or darken (amount > 0) a hex colour by blending.

    S13A needs three tints of each TE class colour, one per TSS window: 5 kb
    faint, 10 kb the palette colour itself, 20 kb dark. Blending rather than an
    HSV shift keeps the hue exactly on the palette value so the class is still
    identifiable at every tint.
    """
    r, g, b = mcolors.to_rgb(colour)
    if amount < 0:  # blend toward white
        f = -amount
        return (r + (1 - r) * f, g + (1 - g) * f, b + (1 - b) * f)
    return (r * (1 - amount), g * (1 - amount), b * (1 - amount))


WINDOW_SHADES = {"5kb": -0.55, "10kb": 0.0, "20kb": 0.35}

# The GO networks colour nodes by -log10(FDR) on `coolwarm`. S13D must share that
# scale rather than approximate it by eye, so both the colormap name and the
# normalizer are defined here and used by every consumer.
GO_FDR_CMAP = "coolwarm"


def go_fdr_norm(fdr_values, floor=None):
    """Normalizer for -log10(FDR) with the publication threshold as vmin.

    `floor` defaults to FDR_THRESHOLD, so the bottom of the colour scale is
    always the significance cut and not whatever the least significant retained
    term happened to be.
    """
    floor = FDR_THRESHOLD if floor is None else floor
    values = -np.log10(np.asarray(fdr_values, dtype=float))
    return mcolors.Normalize(
        vmin=-np.log10(floor), vmax=float(np.nanmax(values))
    )


# --------------------------------------------------------------------------- #
# Revision-wide parameters
# --------------------------------------------------------------------------- #

# Decision D2: 0.05 everywhere, with no "suggestive 0.05-0.1" band anywhere in
# the main text, the supplementary figures or the supplementary tables.
FDR_THRESHOLD = 0.05

# Decision D1: N=500 is what was actually run. The empirical p floor is
# 2/(500+1) = 0.003992, which is the 0.004 quoted in Table 1.
N_PERMUTATIONS = 500
EMPIRICAL_P_FLOOR = 2.0 / (N_PERMUTATIONS + 1)

# §3.3: the interferon-alpha domain, T2T-CHM13 coordinates.
IFNA_DOMAIN = ("chr9", 21_150_692, 21_370_055)

# §3.0: the canonical inputs. These supersede the broken
# `../../T2T_article/T2T_repeat_masker_processed.csv` path and the missing
# `T2T_genes_sorted.bed` documented in the parent CLAUDE.md. Resolved relative
# to the repository directory that holds this package, so a notebook running in
# revision_G3/ and a script running from the parent both find them.
REPO_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REVISION_DIR = os.path.join(REPO_DIR, "revision_G3")
OUTPUT_DIR = os.path.join(REVISION_DIR, "output")
SVG_DIR = os.path.join(REVISION_DIR, "svg")

REPEATS_BED = os.path.join(REPO_DIR, "T2T_repeat_masker_processed_sorted.bed")
GENES_BED = os.path.join(REPO_DIR, "T2T_genes.bed")
OBO_PATH = os.path.join(REPO_DIR, "go-basic.obo")
GAF_PATH = os.path.join(REPO_DIR, "goa_human.gaf")

# Column order of both the compacted legacy store and the new streaming runs.
PERMUTATION_COUNTS_COLUMNS = [
    "seed",
    "score",
    "subfamily_name",
    "family_name",
    "class_name",
    "n",
]


# --------------------------------------------------------------------------- #
# GO enrichment  (copied from Gene_ontology_analysis.ipynb cell 6, FDR -> 0.05)
# --------------------------------------------------------------------------- #


@functools.lru_cache(maxsize=2)
def _load_godag(obo_path: str) -> GODag:
    """Parse the GO DAG once per process. See module docstring, change 3."""
    print(f"Loading GO DAG from {obo_path} ...")
    return GODag(obo_path)


@functools.lru_cache(maxsize=2)
def _load_associations(gaf_path: str) -> dict:
    """Build the {gene symbol -> set of GO IDs} map once per process.

    Same construction as the frozen original: one entry per DB_Symbol, with no
    evidence-code or aspect filtering, so term memberships are identical.
    """
    print(f"Loading human GO annotations from {gaf_path} ...")
    ogaf = GafReader(gaf_path)
    full_assoc: dict[str, set] = {}
    for ntf in ogaf.associations:
        full_assoc.setdefault(ntf.DB_Symbol, set()).add(ntf.GO_ID)
    return full_assoc


def run_goatools_enrichment(
    gene_list,
    background_list,
    obo_path=OBO_PATH,
    gaf_path=GAF_PATH,
    fdr_threshold=FDR_THRESHOLD,
    verbose=False,
):
    """GO enrichment of `gene_list` against `background_list`.

    Copy of `run_goatools_enrichment` from the frozen Gene_ontology_analysis.ipynb
    cell 6, with the default FDR threshold lowered from 0.1 to 0.05 (D2) and the
    DAG/GAF loads hoisted into process-level caches.

    Returns a DataFrame of the terms with `p_fdr_bh < fdr_threshold`, sorted by
    FDR, carrying BOTH the raw `P-value` and the corrected `FDR` column — which
    is what G3's statistics policy and reviewer minor comment 3 require (WP9).
    Returns an empty DataFrame when nothing is significant.
    """
    godag = _load_godag(obo_path)
    full_assoc = _load_associations(gaf_path)

    goeaobj = GOEnrichmentStudy(
        background_list,
        full_assoc,
        godag,
        propagate_counts=False,
        alpha=fdr_threshold,
        methods=["fdr_bh"],
    )

    if verbose:
        print(f"Running analysis on {len(gene_list)} genes ...")
    results_all = goeaobj.run_study(gene_list, prt=None)
    results_sig = [r for r in results_all if r.p_fdr_bh < fdr_threshold]

    parsed_results = []
    for res in results_sig:
        pop_count = res.pop_count
        study_count = res.study_count
        fold_enrichment = (study_count / res.study_n) / (pop_count / res.pop_n)
        parsed_results.append(
            {
                "Term ID": res.GO,
                "Term Name": res.name,
                "Term Database": res.NS,
                "P-value": res.p_uncorrected,
                "FDR": res.p_fdr_bh,
                "Fold Enrichment": fold_enrichment,
                "Overlap Count": study_count,
                "Total Term Genes (Human)": pop_count,
                "Overlapping Genes": sorted(res.study_items),
                "Full Term Gene List": sorted(res.pop_items),
            }
        )

    df = pd.DataFrame(parsed_results)
    return df.sort_values("FDR") if not df.empty else pd.DataFrame()


def run_goatools_ordered_enrichment(
    ordered_gene_list,
    obo_path=OBO_PATH,
    gaf_path=GAF_PATH,
    fdr_threshold=FDR_THRESHOLD,
    verbose=False,
):
    """GO enrichment of a ranked gene list against itself as the population.

    Copy of `run_goatools_ordered_enrichment` from the frozen
    Gene_ontology_analysis.ipynb cell 6, with FDR 0.1 -> 0.05 (D2) and cached
    DAG/GAF loads. `Overlapping Genes` is comma-joined here, as in the original,
    rather than a list — the two functions differ in that respect and the
    downstream readers depend on it.
    """
    godag = _load_godag(obo_path)
    full_assoc = _load_associations(gaf_path)

    goeaobj = GOEnrichmentStudy(
        ordered_gene_list,
        full_assoc,
        godag,
        propagate_counts=False,
        alpha=fdr_threshold,
        methods=["fdr_bh"],
    )

    if verbose:
        print(f"Running analysis on ordered list of {len(ordered_gene_list)} genes ...")
    results_all = goeaobj.run_study(ordered_gene_list, prt=None)
    results_sig = [r for r in results_all if r.p_fdr_bh < fdr_threshold]

    parsed_results = []
    for res in results_sig:
        pop_count = res.pop_count
        study_count = res.study_count
        fold_enrichment = (study_count / res.study_n) / (pop_count / res.pop_n)
        parsed_results.append(
            {
                "Term ID": res.GO,
                "Term Name": res.name,
                "P-value": res.p_uncorrected,
                "FDR": res.p_fdr_bh,
                "Fold Enrichment": fold_enrichment,
                "Overlap Count": study_count,
                "Total Term Genes (Human)": pop_count,
                "Overlapping Genes": ", ".join(sorted(res.study_items)),
                "Full Term Gene List": ", ".join(sorted(res.pop_items)),
            }
        )

    df = pd.DataFrame(parsed_results)
    return df.sort_values("FDR") if not df.empty else pd.DataFrame()


# --------------------------------------------------------------------------- #
# Label-collision enforcement (WP12, caveat C16) — new in this module
# --------------------------------------------------------------------------- #


def find_label_collisions(fig, texts=None, pad_points=0.5):
    """Return the list of overlapping label pairs as (text_a, text_b, area).

    Renders the figure to a canvas so every text artist has a real bounding box,
    then tests all pairs. `pad_points` inflates each box so labels that merely
    touch still count as colliding — at print size, touching labels read as
    overlapping.
    """
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

    if texts is None:
        texts = [t for ax in fig.axes for t in ax.texts]

    boxes = []
    for t in texts:
        if not t.get_visible() or not t.get_text().strip():
            continue
        bbox = t.get_window_extent(renderer=renderer)
        pad_px = pad_points * fig.dpi / 72.0
        boxes.append((t, bbox.padded(pad_px)))

    collisions = []
    for (ta, ba), (tb, bb) in combinations(boxes, 2):
        overlap = ba.intersection(ba, bb)
        if overlap is not None and overlap.width > 0 and overlap.height > 0:
            collisions.append((ta.get_text(), tb.get_text(), overlap.width * overlap.height))
    return collisions


def assert_no_label_collisions(fig, texts=None, pad_points=0.5, label=""):
    """Raise RuntimeError if any two labels in `fig` overlap.

    Enforced rather than attempted: reviewer minor comment 6 is specifically
    about unreadable overlapping text, so a colliding figure must not be
    saveable. If this fires, use the documented fallback ladder (shorten labels
    via output/go_label_shortnames.csv, then leader lines, then reduce top_n) —
    do not loosen the check (caveat C16).
    """
    collisions = find_label_collisions(fig, texts=texts, pad_points=pad_points)
    if collisions:
        worst = sorted(collisions, key=lambda c: -c[2])[:10]
        detail = "\n".join(f"  {a!r} <-> {b!r} ({area:.0f} px^2)" for a, b, area in worst)
        raise RuntimeError(
            f"{len(collisions)} label collision(s) in {label or 'figure'}:\n{detail}"
        )
    return 0


def find_offpage_labels(fig, texts=None, margin_points=0.0):
    """Return (text, overhang_px) for labels that leave the figure canvas.

    Only meaningful when the figure is saved without `bbox_inches="tight"`: a tight
    bounding box grows the page to contain every label, whereas a pinned page crops it.
    Half the width is not achieved by pushing text off the edge, so this is checked.
    """
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    if texts is None:
        texts = [t for ax in fig.axes for t in ax.texts]

    page = fig.bbox
    margin_px = margin_points * fig.dpi / 72.0
    offpage = []
    for text in texts:
        if not text.get_visible() or not text.get_text().strip():
            continue
        box = text.get_window_extent(renderer=renderer)
        overhang = max(page.x0 + margin_px - box.x0, box.x1 - (page.x1 - margin_px),
                       page.y0 + margin_px - box.y0, box.y1 - (page.y1 - margin_px))
        if overhang > 0:
            offpage.append((text.get_text(), overhang))
    return offpage


def assert_labels_on_page(fig, texts=None, margin_points=0.0, label=""):
    """Raise RuntimeError if any label crosses the figure edge.

    The message opens with a count, in the same shape as
    `assert_no_label_collisions`, so the fallback ladder can read it and try the next
    rung rather than treating an off-page label as a fatal error.
    """
    offpage = find_offpage_labels(fig, texts=texts, margin_points=margin_points)
    if offpage:
        worst = sorted(offpage, key=lambda item: -item[1])[:10]
        detail = "\n".join(f"  {text!r} overhangs by {px:.0f} px" for text, px in worst)
        raise RuntimeError(
            f"{len(offpage)} label(s) outside the page in {label or 'figure'}:\n{detail}"
        )
    return 0


def record_offpage_labels(fig, output_file, texts=None, qc=None):
    """Count labels crossing the page edge into `qc`, without raising.

    The waiver path writes with the collision check off, so it is the one place a pinned
    page could crop a term name unnoticed. Counting keeps the waiver's cost measured
    rather than assumed, in the same spirit as `collisions_at_waived_rung`.
    """
    n_offpage = len(find_offpage_labels(fig, texts=texts))
    if qc is not None:
        qc.setdefault(os.path.basename(output_file), {})["labels_off_page"] = n_offpage
    return n_offpage


def assert_svg_geometry(path, target_pt, tol=0.03, qc=None):
    """Raise RuntimeError unless the written SVG matches `target_pt` within `tol`.

    Page geometry is a request in its own right for the supplementary networks — half
    the published width at the same height (S9, S10) and 1.5x the height at the same
    width (S11) — so it gets the same enforced treatment C16 gives label overlap: a
    panel that silently comes out the wrong shape is the failure mode to design out.
    Records the measured size and aspect into `qc` so PLACEMENT.md and the captions can
    be updated from measurement rather than memory.
    """
    with open(path) as handle:
        header = handle.read(4000)
    width = float(re.search(r'width="([\d.]+)pt"', header).group(1))
    height = float(re.search(r'height="([\d.]+)pt"', header).group(1))
    if qc is not None:
        record = qc.setdefault(os.path.basename(path), {})
        record["svg_pt"] = [round(width, 2), round(height, 2)]
        record["aspect"] = round(width / height, 3)
    target_width, target_height = target_pt
    off_width = abs(width - target_width) / target_width
    off_height = abs(height - target_height) / target_height
    if off_width > tol or off_height > tol:
        raise RuntimeError(
            f"{os.path.basename(path)} geometry {width:.2f} x {height:.2f} pt is outside "
            f"{tol:.0%} of the target {target_width:.2f} x {target_height:.2f} pt "
            f"(off by {off_width:.1%} x {off_height:.1%})"
        )
    return width, height


def save_svg_collision_checked(fig, output_file, texts=None, label="", qc=None,
                               tight_bbox=True):
    """Write `fig` to SVG only after asserting zero label collisions.

    `qc`, when given a dict, records {panel: {"label_collisions": 0, "n_texts": N}}
    for `output/network_qc.json` (§10 verification block). `tight_bbox=False` writes at
    `figsize` exactly instead of at the trimmed bounding box, which is what makes a
    pinned canvas assertable against a pt target.
    """
    n_collisions = assert_no_label_collisions(fig, texts=texts, label=label or output_file)
    if not tight_bbox:
        # A pinned page crops what a tight bounding box would have grown to include, so
        # an off-page label is a silent loss of a term name rather than a wider figure.
        assert_labels_on_page(fig, texts=texts, label=label or output_file)
    os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)
    fig.savefig(output_file, format="svg",
                bbox_inches="tight" if tight_bbox else None, transparent=True)
    if qc is not None:
        panel = os.path.basename(output_file)
        qc[panel] = {
            "label_collisions": n_collisions,
            "n_texts": len([t for ax in fig.axes for t in ax.texts]),
        }
    print(f"SVG saved to {output_file} (label collisions: {n_collisions})")
    return output_file


def _record_layout_qc(qc, output_file, figsize, layout_k, label_expand, norm,
                      baseline_figsize=None):
    """Add the layout settings a network panel was actually drawn at to `qc`.

    `canvas_area_vs_baseline` is what makes the "30-50 % denser packing" claim a
    measured number rather than an assertion: it is the drawn canvas area over the
    panel's nominal area, so 0.49 means the panel was laid out at 70 % per side.
    Recorded only when the caller supplies `baseline_figsize`, because only the
    caller knows what the nominal size for that panel is.
    """
    if qc is None:
        return
    record = qc.setdefault(os.path.basename(output_file), {})
    record["figsize_in"] = [round(figsize[0], 2), round(figsize[1], 2)]
    record["layout_k"] = layout_k
    record["label_expand"] = list(label_expand)
    record["fdr_colour_range"] = [round(float(norm.vmin), 3), round(float(norm.vmax), 3)]
    if baseline_figsize is not None:
        baseline_area = baseline_figsize[0] * baseline_figsize[1]
        record["canvas_area_vs_baseline"] = round(
            (figsize[0] * figsize[1]) / baseline_area, 3
        )


# --------------------------------------------------------------------------- #
# GO network plotters
# (copied from cells 6, 36, 175; + min_shared_genes / max_term_genes / QC)
# --------------------------------------------------------------------------- #


def _parse_gene_list(x):
    """Normalise the `Overlapping Genes` column, which is a stringified list in
    the CSVs and a real list in memory. Copied from the frozen cells."""
    if isinstance(x, str):
        try:
            return ast.literal_eval(x)
        except (ValueError, SyntaxError):
            return [g.strip() for g in x.split(",")]
    return x


def _filter_go_frame(
    df,
    group_column,
    fdr_threshold=FDR_THRESHOLD,
    top_n=5,
    max_term_genes=None,
):
    """Apply the FDR cut, the term-size cut and the per-group top-N cut.

    `max_term_genes` is new (WP12): terms annotated to more than that many human
    genes ("protein binding", "cytoplasm", "nucleus", "cytosol", "nucleoplasm")
    carry no interpretation and dominate the layout, so the main-text networks
    exclude them at 500.
    """
    out = df[df["FDR"] < fdr_threshold].copy()
    if max_term_genes is not None:
        out = out[out["Total Term Genes (Human)"] <= max_term_genes]
    if top_n is not None:
        out = out.sort_values([group_column, "FDR"]).groupby(group_column).head(top_n)
    out = out.copy()
    out["Overlapping Genes"] = out["Overlapping Genes"].apply(_parse_gene_list)
    return out


def _jaccard_edges(df_filtered, jaccard_threshold, min_shared_genes):
    """Yield (term_i, term_j, jaccard, n_shared) for term pairs that pass BOTH
    the Jaccard threshold and the new absolute shared-gene floor.

    `min_shared_genes` is new (WP12): a Jaccard of 0.2 between two 5-gene terms
    is one shared gene, which is not a meaningful link but does add an edge that
    the layout has to satisfy.
    """
    unique_terms = df_filtered.drop_duplicates(subset=["Term ID"])
    term_ids = unique_terms["Term ID"].tolist()
    term_genes = [set(g) for g in unique_terms["Overlapping Genes"].tolist()]

    for i in range(len(term_ids)):
        for j in range(i + 1, len(term_ids)):
            set1, set2 = term_genes[i], term_genes[j]
            union_len = len(set1 | set2)
            if union_len == 0:
                continue
            n_shared = len(set1 & set2)
            jaccard = n_shared / union_len
            if jaccard >= jaccard_threshold and n_shared >= min_shared_genes:
                yield term_ids[i], term_ids[j], jaccard, n_shared


def save_go_network_svg(
    df,
    output_file=None,
    jaccard_threshold=0.1,
    top_n=5,
    font_size=GLOBAL_FONT_SIZE,
    palette=None,
    fdr_threshold=FDR_THRESHOLD,
    min_shared_genes=0,
    max_term_genes=None,
    figsize=(16, 14),
    title=None,
    check_collisions=True,
    qc=None,
    layout_k=0.6,
    layout_iterations=None,
    label_expand=(1.15, 1.25),
    label_forces=((0.6, 0.8), (0.4, 0.6)),
    baseline_figsize=None,
    label_wrap=None,
    label_max_move=None,
    leader_lines=False,
    show_colorbar=True,
    show_size_legend=True,
    edge_alpha=0.4,
    edge_width_cap=None,
    tight_bbox=True,
):
    """GO-term network grouped by TE class; Figures 4A / 5A of the manuscript.

    Copy of `save_go_network_svg` from the frozen Gene_ontology_analysis.ipynb
    cell 6. Added: `fdr_threshold` (0.05, D2), `min_shared_genes` and
    `max_term_genes` (WP12), the enforced collision check (C16), and `title`
    so the baked-in panel title can be suppressed (G11 — figure titles belong in
    the manuscript legend, not in the image).

    `layout_k`, `layout_iterations`, `label_expand` and `label_forces` expose the
    two layout stages the compaction ladder tunes; they default to the values the
    published panels used, so a caller that does not pass them gets the previous
    behaviour unchanged.

    The remaining arguments are **presentation only** and change no statistic:
    `label_wrap` breaks a term name over lines of at most that many characters, which is
    the width lever the reference frames use and the only one that buys page width
    without discarding terms; `label_max_move` caps how far `adjust_text` may carry a
    label from its node; `leader_lines` draws a plain capped leader instead of an
    arrowhead; `show_colorbar` / `show_size_legend` drop furniture that would eat a
    narrow canvas (`svg/Fig456A_colorbar_vector.svg` is the vector colourbar, G12);
    `edge_alpha` and `edge_width_cap` clip drawn stroke width only, leaving edge weights
    untouched; `tight_bbox=False` maps `figsize` 1:1 onto the written SVG size, which is
    what makes a pinned canvas exactly assertable. Every one of them defaults to the
    published behaviour, so Figures 4A and 5A are byte-identical unless a caller opts in.

    Node area is proportional to log10(term size); node colour is -log10(FDR);
    class nodes are diamonds in the shared TE class palette.
    """
    output_file = output_file or os.path.join(SVG_DIR, "GO_network_by_classes.svg")
    palette = CLASS_PALETTE if palette is None else palette

    df_filtered = _filter_go_frame(
        df, "class_name", fdr_threshold, top_n, max_term_genes
    )
    if df_filtered.empty:
        raise ValueError(f"No terms survive FDR < {fdr_threshold} with these filters")

    df_filtered["-log10FDR"] = -np.log10(df_filtered["FDR"].replace(0, 1e-15))
    df_filtered["log_pop_size"] = np.log10(
        df_filtered["Total Term Genes (Human)"].replace(0, 1)
    )

    G = nx.Graph()
    class_color_map = dict(palette)
    class_color_map.setdefault("TE_bottom", "red")
    class_color_map.setdefault("TE_top", "green")

    for _, row in df_filtered.iterrows():
        cls = row["class_name"]
        term_id = row["Term ID"]
        if not G.has_node(cls):
            G.add_node(
                cls,
                label=cls,
                type="class",
                color=class_color_map.get(cls, "gray"),
            )
        G.add_node(
            term_id,
            label=row["Term Name"],
            type="go",
            neg_log_fdr=row["-log10FDR"],
            size=row["log_pop_size"] * 150,
            raw_pop=row["Total Term Genes (Human)"],
        )
        G.add_edge(
            cls,
            term_id,
            weight=row["Fold Enrichment"],
            color=class_color_map.get(cls, "#E0E0E0"),
        )

    for t_i, t_j, jaccard, _n in _jaccard_edges(
        df_filtered, jaccard_threshold, min_shared_genes
    ):
        G.add_edge(t_i, t_j, weight=jaccard * 5, color="#E0E0E0", alpha=0.3)

    fig, ax = plt.subplots(figsize=figsize)
    pos = nx.spring_layout(G, k=layout_k, iterations=layout_iterations or 50, seed=42)

    class_nodes = [n for n, d in G.nodes(data=True) if d["type"] == "class"]
    nx.draw_networkx_nodes(
        G,
        pos,
        nodelist=class_nodes,
        node_shape="d",
        node_size=1200,
        node_color=[G.nodes[n]["color"] for n in class_nodes],
        ax=ax,
    )

    go_nodes = [n for n, d in G.nodes(data=True) if d["type"] == "go"]
    go_colors = [G.nodes[n]["neg_log_fdr"] for n in go_nodes]
    go_sizes = [G.nodes[n]["size"] for n in go_nodes]
    # Shared with S13D through go_fdr_norm: vmin is the FDR threshold, not the
    # least significant retained term, so the two panels read one colourbar.
    norm = go_fdr_norm(df_filtered["FDR"], floor=fdr_threshold)
    nx.draw_networkx_nodes(
        G,
        pos,
        nodelist=go_nodes,
        node_shape="o",
        node_size=go_sizes,
        node_color=go_colors,
        cmap=GO_FDR_CMAP,
        vmin=norm.vmin,
        vmax=norm.vmax,
        alpha=0.8,
        ax=ax,
    )

    edges = G.edges(data=True)

    def _edge_width(d):
        width = d["weight"] / 2
        return min(width, edge_width_cap) if edge_width_cap else width

    nx.draw_networkx_edges(
        G,
        pos,
        edgelist=edges,
        width=[_edge_width(d) for _u, _v, d in edges],
        edge_color=[d.get("color", "#E0E0E0") for _u, _v, d in edges],
        alpha=edge_alpha,
        ax=ax,
    )

    if show_colorbar:
        sm = plt.cm.ScalarMappable(cmap=GO_FDR_CMAP, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, shrink=0.5, pad=0.05)
        # WP9/minor comment 3: the statistic's provenance must be on the artist.
        cbar.set_label("-log10(FDR-corrected GO enrichment p-value)", fontsize=font_size)
        cbar.ax.tick_params(labelsize=font_size)

    if show_size_legend:
        for count in [10, 100, 1000]:
            ax.scatter(
                [], [], c="gray", alpha=0.5, s=np.log10(count) * 150,
                label=f"{count} genes",
            )
        ax.legend(
            scatterpoints=1,
            labelspacing=1.2,
            title="Total genes in term",
            loc="lower left",
            frameon=True,
            fontsize=font_size,
            title_fontsize=font_size,
        )

    def _label(node):
        text = G.nodes[node]["label"]
        return textwrap.fill(text, label_wrap) if label_wrap else text

    # Centre anchoring only matters once labels are multi-line: with the default
    # left-baseline anchor a 4-line label sits entirely below-right of its node, which is
    # the opposite of putting a label beside its cognate term.
    # Anchoring and line spacing are only overridden when wrapping is on, so an
    # unwrapped panel keeps matplotlib's defaults and stays byte-identical.
    wrap_kwargs = (dict(ha="center", va="center", linespacing=0.95) if label_wrap
                   else {})
    texts = [
        ax.text(x, y, _label(node), fontsize=font_size, fontweight="normal",
                **wrap_kwargs)
        for node, (x, y) in pos.items()
    ]
    adjust_text(
        texts,
        ax=ax,
        arrowprops=(dict(arrowstyle="-", color="gray", lw=0.4, shrinkA=1, shrinkB=1)
                    if leader_lines
                    else dict(arrowstyle="->", color="gray", lw=0.5)),
        force_text=label_forces[0],
        force_static=label_forces[1],
        expand=label_expand,
        max_move=label_max_move,
    )

    if title:
        ax.set_title(title, fontsize=font_size)
    ax.axis("off")
    plt.tight_layout()

    try:
        if check_collisions:
            save_svg_collision_checked(fig, output_file, texts=texts, qc=qc,
                                       tight_bbox=tight_bbox)
        else:
            if not tight_bbox:
                record_offpage_labels(fig, output_file, texts=texts, qc=qc)
            os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)
            fig.savefig(output_file, format="svg",
                        bbox_inches="tight" if tight_bbox else None, transparent=True)
        _record_layout_qc(qc, output_file, figsize, layout_k, label_expand, norm,
                          baseline_figsize)
    finally:
        # `finally`, because the collision check RAISES and the fallback ladder catches
        # it to try the next rung. Closing only on success leaves every rejected attempt
        # open, so the inline backend renders all of them into the notebook (113 MB, over
        # GitHub's limit) and matplotlib warns about too many open figures.
        plt.close(fig)
    return G


def save_go_network_svg_families_by_classes(
    df,
    family_to_class,
    output_file=None,
    jaccard_threshold=0.1,
    top_n=5,
    font_size=GLOBAL_FONT_SIZE,
    palette=None,
    fdr_threshold=FDR_THRESHOLD,
    min_shared_genes=0,
    max_term_genes=None,
    figsize=(18, 16),
    title=None,
    check_collisions=True,
    qc=None,
    layout_k=0.6,
    layout_iterations=None,
    label_expand=(1.15, 1.25),
    label_forces=((0.6, 0.8), (0.4, 0.6)),
    baseline_figsize=None,
    label_wrap=None,
    label_max_move=None,
    leader_lines=False,
    show_colorbar=True,
    show_size_legend=True,
    show_jaccard_legend=True,
    edge_alpha=0.3,
    edge_width_cap=None,
    tight_bbox=True,
    same_class_weight=10.0,
    legend_rect=(0, 0, 0.85, 1),
):
    """GO-term network grouped by TE family, coloured by class; Figure 6A.

    Copy of `save_go_network_svg_families_by_classes` from the frozen
    Gene_ontology_analysis.ipynb cell 175, with the same four additions as
    `save_go_network_svg`. `family_to_class` is now an explicit argument — the
    frozen version read a notebook-global of that name, which is exactly the
    coupling the freeze is meant to remove.

    Families of the same class are pulled together by high-weight invisible
    edges added to a layout-only copy of the graph, as in the original.
    `same_class_weight` is the weight of those edges: at the original 10.0 it dominates
    both Jaccard weights (≤ 1) and family→term fold enrichments, so every family of a
    class collapses onto one point and drags its terms with it. Lower values keep the
    classes recognisably grouped without the clumping; the default reproduces Figure 6A.

    The presentation-only arguments are the same as in `save_go_network_svg`, plus
    `show_jaccard_legend` and `legend_rect` — the `tight_layout` rect, which at
    (0, 0, 1, 1) reclaims the 15 % right-hand strip the legends are normally hung in and
    therefore moves both legends inside the axes. All of them default to the published
    behaviour, so Figure 6A is byte-identical unless a caller opts in.
    """
    output_file = output_file or os.path.join(
        SVG_DIR, "GO_network_families_by_classes.svg"
    )
    palette = CLASS_PALETTE if palette is None else palette

    df_filtered = _filter_go_frame(
        df, "family_name", fdr_threshold, top_n, max_term_genes
    )
    if df_filtered.empty:
        raise ValueError(f"No terms survive FDR < {fdr_threshold} with these filters")

    df_filtered["-log10FDR"] = -np.log10(df_filtered["FDR"].replace(0, 1e-15))
    df_filtered["log_pop_size"] = np.log10(
        df_filtered["Total Term Genes (Human)"].replace(0, 1)
    )

    G = nx.Graph()
    term_class_contribution: dict = {}

    for _, row in df_filtered.iterrows():
        fam = row["family_name"]
        term_id = row["Term ID"]
        cls = family_to_class.get(fam, "Other")
        n_overlap = len(row["Overlapping Genes"])

        contrib = term_class_contribution.setdefault(term_id, {})
        contrib[cls] = contrib.get(cls, 0) + n_overlap

        if not G.has_node(fam):
            G.add_node(
                fam,
                label=fam,
                type="class",
                color=palette.get(cls, "gray"),
                class_name=cls,
            )
        G.add_node(
            term_id,
            label=row["Term Name"],
            type="go",
            neg_log_fdr=row["-log10FDR"],
            size=row["log_pop_size"] * 150,
            raw_pop=row["Total Term Genes (Human)"],
        )
        G.add_edge(
            fam,
            term_id,
            weight=row["Fold Enrichment"],
            edge_type="fam_go",
            color=palette.get(cls, "#E0E0E0"),
        )

    for t_i, t_j, jaccard, _n in _jaccard_edges(
        df_filtered, jaccard_threshold, min_shared_genes
    ):
        G.add_edge(t_i, t_j, weight=jaccard, edge_type="jaccard", color="#E0E0E0")

    # Layout-only graph: same-class families get a high-weight edge so the
    # spring layout clusters them, without drawing those edges.
    H = G.copy()
    fams = [n for n, d in G.nodes(data=True) if d["type"] == "class"]
    for i in range(len(fams)):
        for j in range(i + 1, len(fams)):
            if G.nodes[fams[i]]["class_name"] == G.nodes[fams[j]]["class_name"]:
                H.add_edge(fams[i], fams[j], weight=same_class_weight)

    fig, ax = plt.subplots(figsize=figsize)
    pos = nx.spring_layout(H, k=layout_k, iterations=layout_iterations or 80, seed=42,
                           weight="weight")

    class_nodes = [n for n, d in G.nodes(data=True) if d["type"] == "class"]
    nx.draw_networkx_nodes(
        G,
        pos,
        nodelist=class_nodes,
        node_shape="d",
        node_size=1200,
        node_color=[G.nodes[n]["color"] for n in class_nodes],
        ax=ax,
    )

    go_nodes = [n for n, d in G.nodes(data=True) if d["type"] == "go"]
    go_colors = [G.nodes[n]["neg_log_fdr"] for n in go_nodes]
    go_sizes = [G.nodes[n]["size"] for n in go_nodes]
    # Shared with S13D and with the class-level networks through go_fdr_norm.
    norm = go_fdr_norm(df_filtered["FDR"], floor=fdr_threshold)
    nx.draw_networkx_nodes(
        G,
        pos,
        nodelist=go_nodes,
        node_shape="o",
        node_size=go_sizes,
        node_color=go_colors,
        cmap=GO_FDR_CMAP,
        vmin=norm.vmin,
        vmax=norm.vmax,
        alpha=0.8,
        ax=ax,
    )

    edges = G.edges(data=True)

    def _width(d):
        width = d["weight"] / 2 if d["edge_type"] == "fam_go" else d["weight"] * 10
        return min(width, edge_width_cap) if edge_width_cap else width

    nx.draw_networkx_edges(
        G,
        pos,
        edgelist=edges,
        width=[_width(d) for _u, _v, d in edges],
        edge_color=[d.get("color", "#E0E0E0") for _u, _v, d in edges],
        alpha=edge_alpha,
        ax=ax,
    )

    # As in save_go_network_svg: centre anchoring and tighter line spacing only when
    # wrapping is on, so an unwrapped panel keeps matplotlib's defaults.
    wrap_kwargs = (dict(ha="center", va="center", linespacing=0.95) if label_wrap
                   else {})
    texts = []
    for node, (x, y) in pos.items():
        data = G.nodes[node]
        if data["type"] == "class":
            t_color, t_weight = "black", "bold"
        else:
            t_weight = "normal"
            contribs = term_class_contribution.get(node, {})
            dom_class = max(contribs, key=contribs.get) if contribs else None
            t_color = palette.get(dom_class, "black") if dom_class else "black"
        label_text = (textwrap.fill(data["label"], label_wrap) if label_wrap
                      else data["label"])
        texts.append(
            ax.text(
                x, y, label_text, fontsize=font_size, color=t_color,
                fontweight=t_weight, **wrap_kwargs,
            )
        )
    adjust_text(
        texts,
        ax=ax,
        arrowprops=(dict(arrowstyle="-", color="gray", lw=0.4, shrinkA=1, shrinkB=1)
                    if leader_lines
                    else dict(arrowstyle="->", color="gray", lw=0.5)),
        force_text=label_forces[0],
        force_static=label_forces[1],
        expand=label_expand,
        max_move=label_max_move,
    )

    if show_colorbar:
        sm = plt.cm.ScalarMappable(cmap=GO_FDR_CMAP, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, shrink=0.4, pad=0.02)
        cbar.set_label("-log10(FDR-corrected GO enrichment p-value)", fontsize=font_size)
        cbar.ax.tick_params(labelsize=font_size)

    # With the legend strip reclaimed there is no right-hand margin left to hang a legend
    # in, so both legends move into the axes corner — anchored outside they would be
    # drawn off-canvas and silently dropped by the pinned (untrimmed) save.
    strip_reclaimed = legend_rect[2] >= 0.99
    size_anchor = (0.0, 0.78) if strip_reclaimed else (1, 0.7)
    jaccard_anchor = (0.0, 0.55) if strip_reclaimed else (1, 0.4)

    if show_size_legend:
        for count in [10, 100, 1000]:
            ax.scatter(
                [], [], c="gray", alpha=0.5, s=np.log10(count) * 150,
                label=f"{count} genes",
            )
        leg1 = ax.legend(
            title="Term population size",
            loc="lower left",
            bbox_to_anchor=size_anchor,
            frameon=False,
            fontsize=font_size,
            title_fontsize=font_size,
        )
        ax.add_artist(leg1)

    if show_jaccard_legend:
        # Capped identically to the drawn edges: a presentation-only width cap that the
        # legend does not know about would make the legend a lie (the Figure 7H lesson).
        j_lines = [
            Line2D([0], [0], color="gray",
                   lw=min(v * 10, edge_width_cap) if edge_width_cap else v * 10,
                   alpha=edge_alpha + 0.1, label=f"J={v:.2f}")
            for v in [0.1, 0.2, 0.5]
        ]
        ax.legend(
            handles=j_lines,
            title="Jaccard similarity",
            loc="lower left",
            bbox_to_anchor=jaccard_anchor,
            frameon=False,
            fontsize=font_size,
            title_fontsize=font_size,
        )

    if title:
        ax.set_title(title, fontsize=font_size + 4, pad=20)
    ax.axis("off")
    plt.tight_layout(rect=list(legend_rect))

    try:
        if check_collisions:
            save_svg_collision_checked(fig, output_file, texts=texts, qc=qc,
                                       tight_bbox=tight_bbox)
        else:
            if not tight_bbox:
                record_offpage_labels(fig, output_file, texts=texts, qc=qc)
            os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)
            fig.savefig(output_file, format="svg",
                        bbox_inches="tight" if tight_bbox else None, transparent=True)
        _record_layout_qc(qc, output_file, figsize, layout_k, label_expand, norm,
                          baseline_figsize)
    finally:
        # See the note in save_go_network_svg: the collision check raises, so closing
        # here rather than on the success path only.
        plt.close(fig)
    return G


def visualize_go_class_network(
    df,
    output_file="go_network.html",
    jaccard_threshold=0.1,
    top_n=5,
    fdr_threshold=FDR_THRESHOLD,
    min_shared_genes=0,
    max_term_genes=None,
    palette=None,
):
    """Interactive pyvis version of the class-level GO network.

    Copy of `visualize_go_class_network` from the frozen
    Gene_ontology_analysis.ipynb cell 36, with the FDR / min_shared_genes /
    max_term_genes filters added and the class palette switched from Tableau to
    the shared TE class palette so the interactive view matches the figures.
    Kept for inspection only — no manuscript figure is produced from it.
    """
    from pyvis.network import Network

    palette = CLASS_PALETTE if palette is None else palette
    df_filtered = _filter_go_frame(
        df, "class_name", fdr_threshold, top_n, max_term_genes
    )

    G = nx.Graph()
    for _, row in df_filtered.iterrows():
        cls = row["class_name"]
        term_id = row["Term ID"]
        term_name = row["Term Name"]

        if not G.has_node(cls):
            G.add_node(
                cls,
                label=cls,
                color=palette.get(cls, "gray"),
                size=40,
                title=f"Group: {cls}",
                shape="diamond",
                borderWidth=2,
            )

        fdr = row["FDR"] if row["FDR"] > 0 else 1e-15
        G.add_node(
            term_id,
            label=term_name[:30] + "..",
            color="#97c2fc",
            size=min(50, -6 * np.log10(fdr)),
            title=(
                f"Full Term: {term_name}<br>ID: {term_id}<br>"
                f"FDR-corrected p: {row['FDR']:.2e}"
            ),
            shape="dot",
        )
        G.add_edge(
            cls,
            term_id,
            width=max(1, row["Fold Enrichment"] / 3),
            color=palette.get(cls, "#E0E0E0"),
            alpha=0.5,
        )

    for t_i, t_j, jaccard, n_shared in _jaccard_edges(
        df_filtered, jaccard_threshold, min_shared_genes
    ):
        G.add_edge(
            t_i,
            t_j,
            width=jaccard * 15,
            color="#e0e0e0",
            title=f"Shared similarity: {jaccard:.2f} ({n_shared} genes)",
        )

    net = Network(height="850px", width="100%", bgcolor="#ffffff", font_color="black")
    net.from_nx(G)
    net.set_options(
        """
    var options = {
      "physics": {
        "barnesHut": {
          "gravitationalConstant": -50000,
          "centralGravity": 0.2,
          "springLength": 220,
          "springConstant": 0.04,
          "damping": 0.09,
          "avoidOverlap": 1
        },
        "minVelocity": 0.75
      },
      "nodes": { "font": { "size": 16, "face": "arial" } }
    }
    """
    )
    net.save_graph(output_file)
    print(f"Interactive network saved to {output_file}")


# --------------------------------------------------------------------------- #
# Compacted permutation store reader (WP1b) — new in this module
# --------------------------------------------------------------------------- #


def permutation_store_dir(window="10kb"):
    """Directory of the compacted permutation counts for a TSS window size."""
    return os.path.join(OUTPUT_DIR, f"permutation_counts_{window}")


def load_permutation_counts(window="10kb", seeds=None, columns=None):
    """Read the compacted permutation store into one DataFrame.

    Replaces the frozen notebook's read of the 6.37 GB
    `consolidated_random_data.csv` (WP1b). Each row is one unique
    (score, subfamily, family, class) tuple within one permutation, with `n`
    giving how many intersecting elements had that tuple — so every downstream
    statistic must treat `n` as a WEIGHT, not as a row count.

    `window` is one of "5kb", "10kb", "20kb". `seeds` restricts to an iterable
    of permutation indices, which is what the convergence figure uses to walk
    1..500 without holding everything in memory at once.
    """
    store = permutation_store_dir(window)
    if not os.path.isdir(store):
        raise FileNotFoundError(
            f"No compacted store at {store}. Run 01b_compact_permutation_results.py "
            f"(10kb) or 01_permutations_stream.sh (5kb / 20kb) first."
        )

    if seeds is None:
        paths = sorted(
            glob.glob(os.path.join(store, "counts_seed_*.tsv.zst")),
            key=lambda p: int(p.rsplit("_", 1)[1].split(".")[0]),
        )
    else:
        paths = [os.path.join(store, f"counts_seed_{s}.tsv.zst") for s in seeds]

    if not paths:
        raise FileNotFoundError(f"No counts_seed_*.tsv.zst files in {store}")

    frames = [read_counts_file(p, columns=columns) for p in paths]
    return pd.concat(frames, ignore_index=True)


def read_counts_file(path, columns=None):
    """Read one `counts_seed_N.tsv.zst` (or plain `.tsv`) into a DataFrame."""
    dtypes = {
        "seed": "int32",
        "score": "int32",
        "subfamily_name": "category",
        "family_name": "category",
        "class_name": "category",
        "n": "int32",
    }
    if path.endswith(".zst"):
        raw = subprocess.run(
            ["zstd", "-dc", path], check=True, capture_output=True
        ).stdout
        handle = io.BytesIO(raw)
    else:
        handle = path

    # `keep_default_na=False` is REQUIRED, not cosmetic. A small number of
    # subfamilies have an empty `family_name` in the RepeatMasker annotation
    # (2,730 rows in seed 1). Without this, pandas reads them as NaN, and every
    # later `groupby(..., observed=True)` silently DROPS those rows, so
    # family-level totals would not sum to the class-level totals.
    df = pd.read_csv(
        handle,
        sep="\t",
        names=PERMUTATION_COUNTS_COLUMNS,
        header=None,
        dtype=dtypes,
        keep_default_na=False,
    )
    return df[columns] if columns else df


def permutation_totals(window="10kb", by="class_name"):
    """Per-permutation totals of intersecting elements, grouped by `by`.

    This is the quantity the enrichment background needs: for each permutation,
    how many shuffled elements of each class (or family) fell in a TSS window.
    Streams one seed at a time, so it runs in a few hundred MB of RAM.
    """
    store = permutation_store_dir(window)
    paths = sorted(
        glob.glob(os.path.join(store, "counts_seed_*.tsv.zst")),
        key=lambda p: int(p.rsplit("_", 1)[1].split(".")[0]),
    )
    if not paths:
        raise FileNotFoundError(f"No counts_seed_*.tsv.zst files in {store}")

    rows = []
    for path in paths:
        df = read_counts_file(path, columns=["seed", by, "n"])
        agg = df.groupby(["seed", by], observed=True)["n"].sum().reset_index()
        rows.append(agg)
    return pd.concat(rows, ignore_index=True)


def load_manifest(window="10kb"):
    """Read the compaction provenance record written by WP1b."""
    with open(os.path.join(permutation_store_dir(window), "MANIFEST.json")) as fh:
        return json.load(fh)
