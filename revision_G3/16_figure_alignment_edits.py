#!/usr/bin/env python
"""Align the manuscript with the figures as exported from Figma on 2026-08-07.

Reads the 260804 working file and writes a NEW file, so the input is never modified:

    Revised_manuscript/T2T_genes_article_G3_revision_260804_manual.docx  (input, read-only)
        -> Revised_manuscript/T2T_genes_article_G3_revision_260807.docx  (output)

Plan: figures_text_alignment_plan_260807.md. Stages, in the order they must run:

    R  the supplementary renumbering S7-S14 -> S6-S13, after deleting the old S6 caption
       whose figure became Figure 8A
    T  Results, Methods, Abstract and Discussion text
    C  figure captions, main and supplementary
    M  the structural move: all interferon-alpha material to a closing Results subsection,
       plus the new Figure 8 caption

Why the helpers below are copied from 13_manuscript_tracked_edits.py rather than imported
    `13_manuscript_tracked_edits` is a script, not a module: importing it would execute all
    2,157 lines of Phase 5 edits against whatever document it opens. Its module name also
    begins with a digit, so a plain `import` cannot name it. The four helpers this script
    needs are therefore copied verbatim, which is the same choice `revision_lib.py` made for
    the frozen notebooks.

Usage
    ~/venvs/collagen_3_11/bin/python 16_figure_alignment_edits.py [--report] [--force]
"""

import argparse
import json
import os
import re
import shutil
import sys

sys.path.insert(0, os.path.expanduser("~/.claude/skills"))

import docx  # noqa: E402
from docx.oxml import OxmlElement  # noqa: E402
from docx.oxml.ns import qn  # noqa: E402
from lxml import etree  # noqa: E402
from word_rewrite_trackchanges import (  # noqa: E402
    delete_paragraph, find_all_p, heading, ins_paragraph_runs, insert_before, is_p,
    make_del_run, make_run, set_revision_identity, style_of, text_of, wrap_del, wrap_ins,
    _rev_attrs,
)

HERE = os.path.dirname(os.path.abspath(__file__))
DOCX_DIR = os.path.join(HERE, "Revised_manuscript")
DEFAULT_IN = os.path.join(DOCX_DIR, "T2T_genes_article_G3_revision_260804_manual.docx")
DEFAULT_OUT = os.path.join(DOCX_DIR, "T2T_genes_article_G3_revision_260807.docx")
# Neither of these may ever be an input: the first is the read-only submitted baseline, the
# second the record of what the scripts produced before Daniil's 2026-08-04 acceptance pass.
FORBIDDEN_INPUTS = ("T2T_genes_article_G3_submitted_baseline_260418.docx",
                    "T2T_genes_article_G3_revision_260803.docx")

REVISION_AUTHOR = "Claude (G3 revision)"
REVISION_DATE = "2026-08-07T00:00:00Z"

# 128 content controls sit inside body paragraphs and one more, the bibliography, is a
# block-level sibling; `body.findall(".//w:sdt")` therefore counts 129.
EXPECTED_SDT = 129

parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--in", dest="source", default=DEFAULT_IN)
parser.add_argument("--out", dest="target", default=DEFAULT_OUT)
parser.add_argument("--report", action="store_true",
                    help="dry run: resolve every edit, write nothing")
parser.add_argument("--force", action="store_true", help="overwrite an existing --out")
args = parser.parse_args()

if os.path.basename(args.source) in FORBIDDEN_INPUTS:
    sys.exit(f"refusing to read {os.path.basename(args.source)} as the input")
if os.path.exists(args.target) and not (args.force or args.report):
    sys.exit(f"{os.path.basename(args.target)} exists; pass --force to overwrite")

report = []


def record(stage, label, ok, outcome=None):
    """Log one edit's three-way outcome.

    On a fresh copy of the input every edit should report `applied`. `already present`
    means the input already carried the edit and is a warning, not the normal case as it
    is for the in-place scripts; `not found` is fatal, so a silent skip is impossible.
    """
    if outcome is None:
        outcome = "applied" if ok else "not found"
    report.append({"stage": stage, "edit": label, "applied": bool(ok), "outcome": outcome})
    return ok


# =======================================================================================
# Citation-safe tracked editing — copied from 13_manuscript_tracked_edits.py
# =======================================================================================
SOFT_MAP = {"\xa0": " ", "‑": "-", "–": "-", "—": "-", "−": "-",
            "‘": "'", "’": "'", "“": '"', "”": '"', "′": "'"}
SOFT_DROP = {"​", "‌", "‍", "﻿"}


def soften(text):
    """(normalised text, index map back into `text`) for tolerant substring matching."""
    characters, indices = [], []
    for position, character in enumerate(text):
        if character in SOFT_DROP:
            continue
        characters.append(SOFT_MAP.get(character, character))
        indices.append(position)
    return "".join(characters), indices


def soft(text):
    return soften(text)[0]


def strip_proof_errors(paragraph):
    """Drop <w:proofErr> markers, which otherwise fragment contiguous runs."""
    for marker in paragraph.findall(qn("w:proofErr")):
        paragraph.remove(marker)


def run_blocks(paragraph):
    """Maximal contiguous groups of <w:r> children, split by anything else.

    A `<w:sdt>` citation control therefore always terminates a block, which is what keeps
    the rewrite below from ever reaching inside one.
    """
    blocks, current = [], []
    for child in list(paragraph):
        if child.tag == qn("w:r"):
            current.append(child)
        elif current:
            blocks.append(current)
            current = []
    if current:
        blocks.append(current)
    return blocks


def block_text(block):
    return "".join(t.text or "" for r in block for t in r.findall(qn("w:t")))


def tracked_replace_safe(paragraph, replacements):
    """Apply (old, new) substring edits as tracked del+ins without touching citations.

    Each replacement is applied inside the first run-block that contains the whole target
    string; `new=""` is a pure deletion. Returns the targets that were not found.
    """
    strip_proof_errors(paragraph)
    missed = []
    for old, new in replacements:
        placed = False
        target = soft(old)
        for block in run_blocks(paragraph):
            text = block_text(block)
            normalised, index_map = soften(text)
            position = normalised.find(target)
            if position < 0:
                continue
            start = index_map[position]
            end = index_map[position + len(target) - 1] + 1
            before, matched, after = text[:start], text[start:end], text[end:]
            base = block[0].find(qn("w:rPr"))
            pieces = []
            if before:
                pieces.append(make_run(before, base=base))
            pieces.append(wrap_del(make_del_run(matched, base=base)))
            if new:
                pieces.append(wrap_ins(make_run(new, base=base)))
            if after:
                pieces.append(make_run(after, base=base))
            anchor = block[-1]
            for piece in pieces:
                anchor.addprevious(piece)
            for run in block:
                paragraph.remove(run)
            placed = True
            break
        if not placed:
            missed.append(old)
    return missed


def delete_paragraph_fully(paragraph):
    """Tracked-delete a paragraph that also contains this revision's own `<w:ins>` runs.

    Needed for the old Figure S6 caption: the G2 sweep bolded its label as an insertion, and
    `delete_paragraph` only converts plain runs, so a bare tracked delete would leave the
    bold "Figure S6." visible and the document would show two S6 captions. An insertion
    removed by the same revision is simply dropped — it was never in the baseline, so there
    is nothing to accept or reject — while the original runs are tracked-deleted as usual.
    """
    for insertion in list(paragraph.iter(qn("w:ins"))):
        parent = insertion.getparent()
        if parent is not None:
            parent.remove(insertion)
    delete_paragraph_safe(paragraph)


def delete_paragraph_safe(paragraph):
    """Tracked-delete a paragraph, including the runs inside its citation controls."""
    delete_paragraph(paragraph)
    for sdt in paragraph.findall(".//" + qn("w:sdt")):
        for run in sdt.findall(".//" + qn("w:r")):
            if run.getparent().tag == qn("w:del"):
                continue
            for text in run.findall(qn("w:t")):
                text.tag = qn("w:delText")
            marker = _rev_attrs(OxmlElement("w:del"))
            run.getparent().replace(run, marker)
            marker.append(run)


def replace_inside_insertions(scope, pairs):
    """Apply (old, new) to text inside `<w:ins>` runs within `scope`, editing in place.

    Needed because `run_blocks` deliberately cannot see inside an insertion — which is what
    protects the citation controls, but also hides every paragraph Phase 5 added, and most
    of this script's targets are in those. Editing our own inserted text directly is the
    correct tracked semantics: the span is already marked as added.

    Returns the per-pair hit counts. Longest target first, so a short target cannot eat a
    prefix of a longer one.
    """
    counts = {old: 0 for old, _ in pairs}
    for old, new in sorted(pairs, key=lambda pair: -len(pair[0])):
        target = soft(old)
        for insertion in list(scope.iter(qn("w:ins"))):
            for run in insertion.findall(qn("w:r")):
                node = run.find(qn("w:t"))
                if node is None or not node.text:
                    continue
                # One left-to-right pass per run, advancing past each replacement, because
                # `new` legitimately contains `old` when an edit appends to a sentence.
                text, offset = node.text, 0
                while True:
                    normalised, index_map = soften(text[offset:])
                    position = normalised.find(target)
                    if position < 0:
                        break
                    start = offset + index_map[position]
                    end = offset + index_map[position + len(target) - 1] + 1
                    text = text[:start] + new + text[end:]
                    offset = start + len(new)
                    counts[old] += 1
                node.text = text
    return counts


def paragraph_by(body, needle, occurrence=0):
    """The nth body paragraph whose visible text contains `needle`, softly matched."""
    target = soft(needle)
    hits = find_all_p(body, lambda text: target in soft(text))
    return hits[occurrence] if len(hits) > occurrence else None


# =======================================================================================
# New here: one edit entry point that reaches both plain runs and our own insertions
# =======================================================================================
def merge_insertion_runs(paragraph):
    """Concatenate adjacent identically-formatted runs inside each `<w:ins>`.

    Word splits a long inserted run at arbitrary points — the Figure S14 caption arrives as
    257 characters followed by 336 — and `replace_inside_insertions` edits one `<w:t>` at a
    time, so a target that straddles such a split can never match. Merging is lossless
    because only runs whose `<w:rPr>` is byte-identical are joined, which is what keeps a
    bold caption label from being absorbed into the plain body that follows it.
    """
    for insertion in paragraph.iter(qn("w:ins")):
        previous = None
        for run in list(insertion.findall(qn("w:r"))):
            node = run.find(qn("w:t"))
            if node is None:
                previous = None
                continue
            if previous is not None and _run_format(previous) == _run_format(run):
                head = previous.find(qn("w:t"))
                head.text = (head.text or "") + (node.text or "")
                insertion.remove(run)
            else:
                previous = run


def _run_format(run):
    properties = run.find(qn("w:rPr"))
    return b"" if properties is None else etree.tostring(properties)


def replace_anywhere(paragraph, pairs):
    """Apply (old, new) whether the target sits in original runs or inside a `<w:ins>`.

    The two mechanisms are disjoint by construction — `run_blocks` stops at an insertion
    boundary, `replace_inside_insertions` only looks inside one — so trying the plain path
    first and the insertion path for whatever it missed covers a paragraph of either kind
    without either pass being able to double-apply.

    Returns per-target replacement counts. The count matters rather than a bare boolean
    because the two paths differ: the plain path replaces one occurrence per call, while
    the insertion path replaces every occurrence in the paragraph in one pass.
    """
    merge_insertion_runs(paragraph)
    counts = {old: 0 for old, _new in pairs}
    missed = tracked_replace_safe(paragraph, pairs)
    for old, _new in pairs:
        if old not in missed:
            counts[old] = 1
    retry = [(old, new) for old, new in pairs if old in missed]
    if retry:
        counts.update(replace_inside_insertions(paragraph, retry))
    return counts


def distinguishing_part(old, new, length=80):
    """The opening of `new` after whatever prefix it shares with `old`.

    Empty when `new` is a prefix of `old` (a pure deletion), in which case there is no
    text whose presence could show the edit had already been made.
    """
    if not new:
        return ""
    shared = 0
    while shared < min(len(old), len(new)) and old[shared] == new[shared]:
        shared += 1
    return new[shared:shared + length]


def edit(body, stage, needle, pairs, occurrence=0):
    """Locate a paragraph by `needle`, then apply every (old, new) pair inside it."""
    paragraph = paragraph_by(body, needle, occurrence)
    if paragraph is None:
        return record(stage, f"locate: {needle[:70]}", False)
    pending, skipped = [], []
    for old, new in pairs:
        # An append-style edit leaves its own search string in place, so a target whose
        # result is already there must be skipped rather than applied a second time. The
        # test has to be made on the part of `new` that differs from `old`: most of these
        # edits extend a sentence, so `new` opens with `old` and any prefix of `new` would
        # match the unedited text and report a false "already present".
        marker = distinguishing_part(old, new)
        if marker and soft(marker) in soft(text_of(paragraph)):
            skipped.append(old)
        else:
            pending.append((old, new))
    for old in skipped:
        record(stage, f"replace: {old[:70]}", True, "already present")
    if not pending:
        return True
    counts = replace_anywhere(paragraph, pending)
    unresolved = False
    for old, _new in pending:
        if counts[old]:
            record(stage, f"replace: {old[:70]}", True)
        else:
            record(stage, f"replace: {old[:70]}", False)
            unresolved = True
    return not unresolved


def delete_children_tracked(paragraph, children):
    """Tracked-delete a run of paragraph children of mixed revision state.

    Needed for the tail of the Figure 5A paragraph, which G13's gene-symbol italicisation
    fragmented into eight `<w:del>`/`<w:ins>` pairs interleaved with plain runs. The three
    cases are not the same edit:
      - `<w:del>` is already deleted; leaving it is what lets Reject All restore the
        original sentence;
      - `<w:ins>` is this revision's own insertion, so it is dropped outright — there is
        nothing for a reviewer to accept or reject in text that was never in the baseline;
      - a plain `<w:r>` is baseline text and must be wrapped in a tracked deletion.
    """
    for child in children:
        if child.tag == qn("w:ins"):
            paragraph.remove(child)
        elif child.tag == qn("w:r"):
            for node in child.findall(qn("w:t")):
                node.tag = qn("w:delText")
            marker = _rev_attrs(OxmlElement("w:del"))
            child.addprevious(marker)
            paragraph.remove(child)
            marker.append(child)


# =======================================================================================
# Open the copy
# =======================================================================================
if not args.report:
    shutil.copyfile(args.source, args.target)
document = docx.Document(args.source if args.report else args.target)
body = document.element.body
set_revision_identity(REVISION_AUTHOR, REVISION_DATE)

sdt_before = len(body.findall(".//" + qn("w:sdt")))
if sdt_before != EXPECTED_SDT:
    sys.exit(f"input has {sdt_before} content controls, expected {EXPECTED_SDT}")

# =======================================================================================
# Stage R — the supplementary renumbering
#
# Old Supplementary Figure 6, the UCSC view of the interferon-alpha domain, became Figure 8A
# during layout, so every later supplementary figure moves down one. Its caption is deleted
# first: a tracked deletion turns the text into <w:delText>, which the renumbering below
# cannot see, so the vacated number cannot collide with the S7 that is about to take it.
# =======================================================================================
old_s6 = paragraph_by(body, "Figure S6. UCSC Genome Browser visualization")
if old_s6 is None:
    record("R", "C-7 locate the old Figure S6 caption", False)
else:
    delete_paragraph_fully(old_s6)
    record("R", "C-7 delete the old Figure S6 caption (its figure is now Figure 8A)", True)

# Counts measured from the input, so a missed reference is a hard failure rather than a
# silent partial renumber.
RENUMBER = [
    ("Figures S12 and S13", "Figures S11 and S12", 1),
    ("Figure S8A", "Figure S7A", 2),
    ("Figure S8B", "Figure S7B", 2),
    ("Figure S8C", "Figure S7C", 1),
    ("Figure S10", "Figure S9", 2),
    ("Figure S11", "Figure S10", 2),
    ("Figure S12", "Figure S11", 1),
    ("Figure S13", "Figure S12", 1),
    ("Figure S14", "Figure S13", 2),
    ("Figure S7", "Figure S6", 2),
    ("Figure S8", "Figure S7", 1),
    ("Figure S9", "Figure S8", 4),
]


def renumber_supplementary(body, mapping):
    """Renumber every supplementary reference in one collision-proof two-phase pass.

    A straight search-and-replace cannot do this in any order. Descending fails because
    `Figure S9 -> Figure S8` writes text that the next pair, `Figure S8 -> Figure S7`,
    then rewrites again; ascending fails symmetrically. So phase one replaces each old
    label with a unique sentinel and phase two expands the sentinels, which makes the
    result independent of order. The sentinel is plain ASCII because `soften()` normalises
    the text it searches and would strip anything more exotic.

    Longest target first within each phase, so `Figure S8A` is consumed before the
    `Figure S8` that is its prefix.
    """
    ordered = sorted(mapping, key=lambda row: -len(row[0]))
    totals = {old: 0 for old, _new, _n in ordered}
    for index, (old, _new, _expected) in enumerate(ordered):
        sentinel = f"@@SF{index}@@"
        for paragraph in list(body.iter(qn("w:p"))):
            # Loop because the plain path replaces one occurrence per call; the insertion
            # path replaces every occurrence in the paragraph at once, which is why the
            # returned count is added rather than incremented by one.
            while True:
                replaced = replace_anywhere(paragraph, [(old, sentinel)])[old]
                if not replaced:
                    break
                totals[old] += replaced
    for index, (old, new, _expected) in enumerate(ordered):
        sentinel = f"@@SF{index}@@"
        for paragraph in list(body.iter(qn("w:p"))):
            replace_inside_insertions(paragraph, [(sentinel, new)])
    return totals


hits = renumber_supplementary(body, RENUMBER)
for old, new, expected in RENUMBER:
    record("R", f"renumber {old} -> {new} ({hits[old]} of {expected})",
           hits[old] == expected)

# =======================================================================================
# Stage T — Results, Methods, Abstract, Discussion
# =======================================================================================
# T-1  Figure 4A carries 10 terms per group, not 30 (output/network_qc.json).
#
# No "the full network is Figure S8" pointer is added here, nor in T-2 and T-3. It would be
# the first citation of S8, S9 and S10 anywhere in the body and would therefore place them
# ahead of S6 and S7, breaking the supplementary citation order. The Figure 4, 5 and 6
# captions already carry those cross-references.
edit(body, "T", "The integrative network visualization of top 30", [(
    "The integrative network visualization of top 30 the most significant terms per each "
    "of the remaining group (Figure 4A)",
    "The integrative network visualisation of up to ten of the most significant terms per "
    "remaining group (Figure 4A)")])

# T-2  Figure 5A carries 10 terms and a 500-gene term-size limit; 1,000 is the
#      supplementary setting. The two halves sit in different runs, either side of the
#      inserted File S3 citation, so they are two replacements rather than one.
edit(body, "T", "We visualized top 30 GO terms of each TE-divergence groups", [
    ("We visualized top 30 GO terms of each TE-divergence groups",
     "We visualised up to ten GO terms for each TE-divergence group"),
    ("filtered by no more than 1000 genes per term to avoid too general classification "
     "(Figure 5A).",
     "filtered by no more than 500 genes per term to avoid overly general classification "
     "(Figure 5A)."),
])

# T-3  Figure 6A carries five terms per family after the 1.2x font increase (edit M10's
#      body half, applied to the caption in Phase 5 but never here).
edit(body, "T", "Visualization of top 30 GO terms by family", [(
    "Visualization of top 30 GO terms by family in Figure 6A showed a high degree of "
    "functional distinction between processes by families.",
    "Visualisation of up to five GO terms per family in Figure 6A showed a high degree of "
    "functional distinction between processes by families.")])

# T-4  Figure 8 was built with four panels, not the three the plan anticipated.
edit(body, "T", "The null distributions and the subfamily composition", [(
    "The null distributions and the subfamily composition are shown in Figure 8.",
    "The domain is shown in Figure 8A, the null distributions of the three matched tests "
    "in Figure 8B, the L1 subfamily composition with the divergence of every individual "
    "element in Figure 8C, and the leave-one-out means in Figure 8D.")])

# T-6  The convergence panels are B-D of the renumbered figure; A is the stability panel
#      moved in from the old S13E.
edit(body, "T", "The convergence trajectories are shown in Figure S13", [(
    "The convergence trajectories are shown in Figure S13 and the checkpoint values",
    "The convergence trajectories are shown in Figure S13B-D and the checkpoint values")])

# T-7  Point the gene-set and rank stability sentences at the panel that now carries them.
edit(body, "T", "the gene rankings agree with Kendall tau between 0.48 and 0.79", [
    ("the gene rankings agree with Kendall tau between 0.48 and 0.79.",
     "the gene rankings agree with Kendall tau between 0.48 and 0.79 (Figure S13A)."),
    ("All comparisons are given in Figures S11 and S12 and in File S5.",
     "All comparisons are given in Figures S11, S12 and S13A and in File S5."),
])

# T-8  hAT-Charlie MHC class I: the grid is harsher than the percentile arm alone
#      (output/go_grid_headline_by_condition.csv).
edit(body, "T", "the exception is the hAT-Charlie MHC class I association", [(
    "which is significant only in the 5 % set.",
    "which is significant only in the 5 % set and survives only one of the six window x "
    "cut-off conditions (Figure S12D); it is reported here as the weakest of the "
    "family-level associations rather than as an established one.")])

# T-9  Figure 3B shows Helitrons at K-S p = 0.222; the sentence read as though the pattern
#      were uniform across all six classes.
edit(body, "T", "The same pattern was observed for individual TE classes", [(
    "The same pattern was observed for individual TE classes (Figure 3B), with two peaks "
    "found in SINEs and SVA elements.",
    "The same pattern was observed for individual TE classes (Figure 3B), with two peaks "
    "found in SINEs and SVA elements; the difference is significant for five of the six "
    "classes, Helitrons being the exception (K-S p = 0.222, raw).")])

# T-10 The Abstract rounded IFNA test 3 to 0.002 where the Results, Figure 8B and File S4
#      all give 0.0017.
edit(body, "T", "matched for L1 count (empirical p = 0.006) or gene density",
     [("or gene density (p = 0.002).", "or gene density (p = 0.0017).")])

# T-11 The only Figure 8 citation in the Discussion.
edit(body, "T", "The interferon-alpha domain on chromosome 9 is the strongest single locus",
     [("matched for L1 count and for gene density.",
       "matched for L1 count and for gene density (Figure 8).")])

# =======================================================================================
# Stage C — figure captions
# =======================================================================================
# C-1  Figure 2D-2F are family-level, Helitron and SVA are single points, and the placed
#      panels use a `* (ns)` bracket that the published star key contradicts.
edit(body, "C", "For all group comparisons in boxplots", [(
    "Pearson correlation p-values on scatterplots are not corrected for multiple "
    "hypotheses.",
    "Pearson correlation p-values on scatterplots are not corrected for multiple "
    "hypotheses. Panels (D), (E) and (F) compare the 44 TE families grouped by class. The "
    "Helitron and SVA classes contain a single family each and are therefore drawn as "
    "individual points rather than boxes. Mann-Whitney U p-values are Benjamini-Hochberg "
    "corrected across the testable class pairs; a bracket marked \"* (ns)\" denotes a pair "
    "significant on the raw p-value that does not survive the correction.")])

# C-2  Figure 5's caption omits the edge filter that Figure 4's carries, though 4A and 5A
#      use the same settings (nb06 SIMPLE).
edit(body, "C", "Figure 5. Functional analysis of genes whose TSS are enriched", [(
    "GO terms with more than 500 genes were excluded to avoid overly general terms.",
    "GO terms with more than 500 genes were excluded to avoid overly general terms, and "
    "edges require a Jaccard index of at least 0.2 and at least 5 shared genes.")])

# C-3  Figure 6's caption states neither filter.
edit(body, "C", "Figure 6. Functional analysis of genes with TSS enriched", [(
    "the full network is Figure S10.",
    "the full network is Figure S10. GO terms with more than 500 genes were excluded to "
    "avoid overly general terms, and edges require a Jaccard index of at least 0.2 and at "
    "least 5 shared genes.")])

# C-4  Panel 7G: the Results quote the exact value, the caption rounded it to < 0.001.
edit(body, "C", "Figure 7. Analysis of main factors impacting",
     [("(G) Mann-Whitney U, raw p < 0.001;", "(G) Mann-Whitney U, raw p = 1.2 x 10-6;")])

# C-6  Supplementary S1 names the K-S test on the artwork but not in its caption, and one
#      of the six class tests is not significant.
edit(body, "C", "Figure S1. Ridge plots for length distribution", [(
    "(A) for all classes and (B) for individual classes.",
    "(A) for all classes and (B) for individual classes. Distributions were compared by "
    "two-sample Kolmogorov-Smirnov test. The p-values shown are raw, as each panel is a "
    "single test; panel (B) reports one test per class, six in total, of which five are "
    "significant.")])

# C-8  State the FULL network settings completely: min_shared_genes = 0 (nb06 FULL).
edit(body, "C", "Figure S8. Complete connection map of GO terms enriched", [(
    "with edges drawn at a Jaccard index of at least 0.1 and terms of at most 1,000",
    "with edges drawn at a Jaccard index of at least 0.1, with no shared-gene floor, and "
    "terms of at most 1,000")])

# C-9  The panel gained a term when it was re-laid out on the pinned canvas on 2026-08-07
#      (network_qc.json, top_n = 30), so 29 is simply wrong.
edit(body, "C", "Figure S9. Complete connection map of GO terms for the divergence",
     [("Up to 29 terms per group are shown", "Up to 30 terms per group are shown")])

# C-10 Replace the qualitative overlap statement with the measured one.
edit(body, "C", "Figure S10. Complete connection map of GO terms per TE family", [(
    "the label field is saturated and some labels overlap",
    "the label field is saturated and seven pairs of labels overlap")])

# C-11 Panel E moved out of this figure and into S13 as its panel A.
edit(body, "C", "Figure S12. Concordance of the enrichment analysis", [(
    " (E) Overlap coefficient of the top 5 % gene sets and Kendall correlation of the gene "
    "rankings between window pairs, with 95 % bootstrap confidence intervals.", "")])

# C-12 Retitle and re-letter for four panels: A is the stability panel moved in from the
#      old S13E, D the 44-family drift panel Daniil is placing manually (decision D-2).
#      Split into three replacements because each must lie within a single run: the caption
#      label is a separate bold run after the G2 sweep.
edit(body, "C", "Figure S13. Convergence of the permutation background", [
    ("Convergence of the permutation background at 500 permutations.",
     "Stability of the gene sets across TSS window sizes, and convergence of the "
     "permutation background at 500 permutations."),
    ("(A) Running mean of the random odds ratio for each TE class as permutations "
     "accumulate. (B) Running standard deviation of the same quantity. (C) Drift of the "
     "running mean for all 44 families, expressed in units of the final standard "
     "deviation, so that families with different absolute odds ratios are comparable.",
     "(A) Overlap coefficient of the top 5 % gene sets and Kendall correlation of the gene "
     "rankings with 95 % bootstrap confidence intervals, between each pair of window "
     "sizes, for each TE class and for all TEs. The Kendall correlations are raw, as each "
     "is a single test. (B) Running mean of the random odds ratio for each TE class as "
     "permutations accumulate, as a fraction of its value at N = 500; the shaded band is "
     "+/- 1 %. (C) Running standard deviation of the same quantity, on the same scale. "
     "(D) Drift of the running mean for all 44 families, expressed in units of the final "
     "standard deviation, so that families with different absolute odds ratios are "
     "comparable."),
    ("for the worst class and within 0.10 for the worst family",
     "for the worst-behaved class and within 0.10 for the worst-behaved family"),
])

# C-13 The extended discussion is submitted but nothing names it. Note it is cited ONCE in
#      the body, not twice: the second mention the plan listed is editor note [160], which
#      is deleted before submission anyway.
edit(body, "C", "The wider argument for building compound TE-derived gene expression", [(
    "is developed in the extended discussion provided as supplementary material.",
    "is developed in the extended discussion provided as File S6.")])

file_s5 = paragraph_by(body, "File S5. Sensitivity and robustness analyses")
if file_s5 is None:
    record("C", "C-13 locate the File S5 entry", False)
else:
    insert_after_target = file_s5
    caption_runs = [make_run(
        "File S6. Extended discussion. The window-size literature review, the per-class and "
        "per-family mechanistic review and the TE-derived biomarker material, moved out of "
        "the Discussion at the reviewers' request and provided in full.",
        base=None)]
    new_entry = ins_paragraph_runs(caption_runs, style=style_of(file_s5))
    insert_after_target.addnext(new_entry)
    record("C", "C-13 add the File S6 entry", True)

# =======================================================================================
# Stage M — all interferon-alpha material to a closing Results subsection (decision D-3)
#
# Figure 8 was first cited between Figures 5 and 6, so the main-text citation order ran
# 1-2-3-4-5-8-6-7. Moving only the "in detail" block would have left the Figure 8A citation
# stranded in the Figure 5A paragraph and the ordering unfixed, so the trailing sentence of
# that paragraph moves too.
# =======================================================================================
GENE_SYMBOLS = ["IFNA10", "IFNA16", "IFNA17", "IFNA21", "IFNA4", "IFNA6", "IFNA7", "IFNW1"]

fig5a_paragraph = paragraph_by(body, "These immune system GO terms have been sharing")
# Matched by prefix, not equality: `text_of` concatenates deleted and inserted text, so
# this heading reads "DiscussionDISCUSSION" — Phase 5 replaced the all-caps original and
# the tracked deletion is still carried in the paragraph.
discussion = None
for candidate in find_all_p(body, lambda t: t.strip().startswith("Discussion")):
    if style_of(candidate) in ("Heading1", "Heading 1"):
        discussion = candidate
        break

lead = paragraph_by(body, "The interferon alpha domain in detail. Because this locus")
outliers = paragraph_by(body, "The signal is not carried by a handful of outliers")
testing = paragraph_by(body, "Formal testing against matched random windows")
fig7_caption = paragraph_by(body, "Figure 7. Analysis of main factors impacting")

# `is None` rather than a truth test: an lxml element with no children is falsy, so
# `all([...])` on paragraph elements reports a perfectly good paragraph as missing.
_anchors = {"Figure 5A paragraph": fig5a_paragraph, "Discussion heading": discussion,
            "interferon lead": lead, "outliers": outliers, "formal testing": testing,
            "Figure 7 caption": fig7_caption}
_absent = [name for name, element in _anchors.items() if element is None]
if _absent:
    record("M", f"locate the interferon-alpha block ({', '.join(_absent)})", False)
else:
    # T-5b, part one: everything from the italicised IFNA10 onwards is the sentence being
    # moved. Capture the children by identity BEFORE the replacement below re-splits the
    # run that precedes them.
    children = list(fig5a_paragraph)
    first = next((i for i, child in enumerate(children)
                  if "IFNA10" in "".join(t.text or "" for t in child.iter(qn("w:delText")))),
                 None)
    if first is None:
        record("M", "T-5b locate the moved sentence's tail", False)
    else:
        tail = children[first:]
        missed = tracked_replace_safe(fig5a_paragraph, [(
            "These immune system GO terms have been sharing the same core interferon gene "
            "set, namely ",
            "These immune system GO terms all share a single core interferon gene set, "
            "which is characterised in detail at the end of this section.")])
        if missed:
            record("M", "T-5b replace the lead-in of the moved sentence", False)
        else:
            delete_children_tracked(fig5a_paragraph, tail)
            record("M", f"T-5b split the Figure 5A paragraph ({len(tail)} tail elements "
                        f"tracked-deleted)", True)

    # T-5, part one: the run-in bold lead becomes the subsection heading, so it goes.
    if replace_inside_insertions(
            lead, [("The interferon alpha domain in detail. Because this locus",
                    "Because this locus")]).get(
                        "The interferon alpha domain in detail. Because this locus"):
        record("M", "T-5 remove the run-in lead now that it is a heading", True)
    else:
        record("M", "T-5 remove the run-in lead", False)

    # The opening paragraph is fresh text, so its gene symbols do not inherit G13's
    # italics and are built italic run by run.
    opening = ["These immune system GO terms share the same core interferon gene set, "
               "namely "]
    for position, symbol in enumerate(GENE_SYMBOLS):
        opening.append(symbol)
        if position < len(GENE_SYMBOLS) - 2:
            opening.append(", ")
        elif position == len(GENE_SYMBOLS) - 2:
            opening.append(" and ")
    opening.append(
        ", whose TSS neighbourhoods lie in the interferon alpha domain of chromosome 9 "
        "(coordinates 21,150,692 to 21,370,055, a 220 kb region), with per-gene average "
        "divergence of the intersecting LINE elements between 95.0 and 161.7, while the "
        "mean divergence of the individual L1 elements across the whole domain is 135.7 "
        "against a genome-wide L1 mean of 188.2 (File S1, sheet TSS_TE_intersections, "
        "Figure 8A).")
    opening_runs = [make_run(piece, italic=piece in GENE_SYMBOLS) for piece in opening]
    opening_paragraph = ins_paragraph_runs(opening_runs, style=style_of(outliers))

    figure8_caption = ins_paragraph_runs([make_run(
        "Figure 8. The interferon-alpha domain of chromosome 9 and the evolutionary age of "
        "its L1 elements. (A) UCSC Genome Browser view of the 220 kb domain "
        "(chr9:21,150,692-21,370,055, T2T-CHM13v2.0/hs1) showing the curated NCBI RefSeq "
        "gene set and the RepeatMasker annotation. (B) Null distributions of the mean L1 "
        "divergence in random 220 kb autosomal windows under three progressively "
        "better-matched null models: windows containing at least one L1 element (top, "
        "10,000 windows), at least 40 L1 elements (middle, 10,000 windows, controlling for "
        "local L1 density) and at least 10 annotated genes and at least one L1 (bottom, "
        "3,582 windows, controlling for gene density). The orange line marks the observed "
        "domain mean of 135.7 and the dashed line the null mean. The empirical p-values "
        "shown are two-sided and raw, as each is a single test. (C) L1 subfamily "
        "composition of the domain: the left bars give the number of elements in each of "
        "the 36 subfamilies and the right strip the divergence of every individual L1 "
        "element, both coloured by whether the subfamily is young primate-specific (L1HS "
        "and L1P*) or older mammalian (L1M*). The dashed line marks the genome-wide L1 "
        "mean of 188.2 and the orange line the domain mean of 135.7. Young "
        "primate-specific copies are enriched in the domain relative to the rest of the "
        "genome (38 of 77 against 133,450 of 545,659; Fisher exact odds ratio 3.01, raw "
        "p = 3.2 x 10-6). (D) Mean L1 divergence recomputed with each of the 77 elements "
        "removed in turn, with the mean after removing the five youngest elements marked "
        "at 143.8. Every subset remains below the L1-count-matched null mean of 189.2, so "
        "the deficit is not carried by a small number of outliers.")],
        style=style_of(fig7_caption))

    subsection = [heading("The interferon-alpha domain in detail", level=2),
                  opening_paragraph, lead, outliers, testing, figure8_caption]
    for element in (lead, outliers, testing):
        if element.getparent() is not body:
            record("M", "T-5 the interferon block is not a direct child of the body", False)
            break
    else:
        for element in (lead, outliers, testing):
            body.remove(element)
        insert_before(discussion, subsection)
        record("M", "T-5 move the interferon-alpha block to a closing Results subsection",
               True)

# =======================================================================================
# Report
# =======================================================================================
sdt_after = len(body.findall(".//" + qn("w:sdt")))
record("verify", f"Mendeley content controls {sdt_before} -> {sdt_after}",
       sdt_after == EXPECTED_SDT)

if not args.report:
    document.save(args.target)

outcomes = {"applied": 0, "already present": 0, "not found": 0}
for row in report:
    outcomes[row["outcome"]] += 1
failed = [row for row in report if not row["applied"]]

print(f"{len(report)} edits: {outcomes['applied']} applied, "
      f"{outcomes['already present']} already present, {outcomes['not found']} NOT FOUND")
for row in failed:
    print(f"  NOT FOUND [{row['stage']}] {row['edit']}")
if outcomes["already present"]:
    print("  warning: 'already present' on a fresh copy means the input already had the edit")
if not args.report:
    print(f"wrote {os.path.relpath(args.target, HERE)}")

with open(os.path.join(HERE, "output", "figure_alignment_edit_report.json"), "w") as handle:
    json.dump({"source": os.path.basename(args.source),
               "target": os.path.basename(args.target),
               "dry_run": args.report, "outcomes": outcomes,
               "sdt_before": sdt_before, "sdt_after": sdt_after,
               "edits": report}, handle, indent=2)
print("report -> output/figure_alignment_edit_report.json")

if failed:
    sys.exit(1)
