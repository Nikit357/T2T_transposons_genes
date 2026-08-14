#!/usr/bin/env python
"""18_review_edits_260809.py — the 2026-08-09 scientific review, applied as tracked changes.

Input :  Revised_manuscript/T2T_genes_article_G3_revision_260807_manual.docx   (read-only)
Output:  Revised_manuscript/T2T_genes_article_G3_revision_260809.docx

Deterministic rebuild from its input, like 16_..., not an accumulation of in-place edits, so
it may be re-run at will.

WHAT IT DOES
  A. Acknowledges and removes the three 'Daniil to Claude' comments by doing what they ask.
  B. Rewrites 'The interferon-alpha domain' as a self-contained analysis, in five passages,
     every number pointed at a File S4 sheet or a Figure 8 panel.
  C. Rewrites the prior-work comparison, citing File S4 sheets one at a time, and adds Table 3.
  D. Rewrites the sensitivity subsection, citing File S5 sheets and individual S11-S13 panels.
  E. Corrects four numbers that do not reproduce from the source tables.
  F. Renumbers the supplementary figures by order of first citation and rewrites all 13 legends.
  G. Removes the remaining process scaffolding and the proofless statements.
  H. House style: US spelling, plain-text exponents.

CITATION SAFETY
  Only one rewritten paragraph carries a Mendeley content control (the Lu et al. paragraph, one
  <w:sdt>). It is never deleted and never rebuilt: the two edits to it lie wholly before and
  wholly after the control. No other paragraph touched here contains an <w:sdt>, which is
  asserted before any mutation. Reference NUMBERS are not touched anywhere - Mendeley owns them.

Run:  ~/venvs/collagen_3_11/bin/python revision_G3/18_review_edits_260809.py
      ~/venvs/collagen_3_11/bin/python revision_G3/18_review_edits_260809.py --report
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

sys.path.insert(0, "/home/jovyan/.claude/skills")

from docx import Document                                    # noqa: E402
from docx.oxml.ns import qn                                  # noqa: E402
from docx.oxml import OxmlElement                            # noqa: E402
from nar_review_tools import safe_tracked_replace            # noqa: E402  NOT tracked_replace
from word_rewrite_trackchanges import (                      # noqa: E402
    _rev_attrs, delete_paragraph, ins_paragraph, insert_after, make_run,
    set_revision_identity, wrap_ins,
)

HERE = Path(__file__).resolve().parent
SRC = HERE / "Revised_manuscript" / "T2T_genes_article_G3_revision_260807_manual.docx"
DST = HERE / "Revised_manuscript" / "T2T_genes_article_G3_revision_260809.docx"

AUTHOR = "Claude (scientific review 2026-08-09)"
DATE = "2026-08-09T00:00:00Z"

W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"

report: list[tuple[str, str, str]] = []          # (edit id, status, detail)


# ----------------------------------------------------------------- locating helpers

def visible(p) -> str:
    """Text as Word displays it with markup hidden: original + insertions, deletions out."""
    parts = []
    for t in p.iter(qn("w:t")):
        anc, dead = t.getparent(), False
        while anc is not None and anc.tag != qn("w:p"):
            if anc.tag == qn("w:del"):
                dead = True
                break
            anc = anc.getparent()
        if not dead:
            parts.append(t.text or "")
    return "".join(parts)


def find_p(body, needle: str, must_be_unique: bool = True):
    """Resolve a text anchor to exactly one paragraph element. Raises if ambiguous."""
    hits = [p for p in body.iter(qn("w:p")) if needle in visible(p)]
    if not hits:
        raise LookupError(f"anchor not found: {needle!r}")
    if must_be_unique and len(hits) > 1:
        raise LookupError(f"anchor matches {len(hits)} paragraphs, lengthen it: {needle!r}")
    return hits[0]


def assert_no_citation(p, where: str):
    if p.find(".//" + qn("w:sdt")) is not None:
        raise AssertionError(f"refusing to touch a paragraph holding a citation control: {where}")


OUR_AUTHORS = ("Claude",)          # insertions we are allowed to retype in place


def _retype_in_own_ins(p, old: str, new: str) -> bool:
    """Edit text that lives inside an unaccepted <w:ins> of OUR OWN earlier pass.

    safe_tracked_replace cannot see such text: those runs are children of the w:ins, not of
    the paragraph. Retyping mirrors what Word does when you edit your own not-yet-accepted
    insertion, and Reject All still restores the original because the whole insertion is
    dropped. Never do this inside another author's revision.
    """
    for ins in p.findall(qn("w:ins")):
        author = ins.get(qn("w:author")) or ""
        if not author.startswith(OUR_AUTHORS):
            continue
        ts = [t for r in ins.findall(qn("w:r")) for t in r.findall(qn("w:t"))]
        joined = "".join(t.text or "" for t in ts)
        if old not in joined:
            continue
        start = joined.index(old)
        end = start + len(old)
        acc, done = 0, False
        for t in ts:
            txt = t.text or ""
            lo, hi = acc, acc + len(txt)
            acc = hi
            if hi <= start or lo >= end:
                continue
            head = txt[: max(0, start - lo)]
            tail = txt[max(0, end - lo):] if hi > end else ""
            t.text = head + ("" if done else new) + tail
            done = True
        return True
    return False


def spans_revision(p, old: str) -> bool:
    """Would replacing `old` straddle an existing <w:ins>/<w:del> inside this paragraph?

    safe_tracked_replace matches across DIRECT-CHILD runs only, and rebuilds the span by
    inserting the replacement before the first matched run and deleting runs first..last.
    Any revision element sitting BETWEEN those runs is not removed, but it ends up after the
    replacement instead of inside it - so no text is lost, yet Reject All returns the
    already-deleted fragment in a different position. Refuse such an anchor and pick a
    shorter one that lies inside a single uninterrupted stretch.
    """
    kids = [c for c in p if c.tag in (qn("w:r"), qn("w:ins"), qn("w:del"))]
    runs, texts = [], []
    for i, c in enumerate(kids):
        if c.tag == qn("w:r"):
            runs.append(i)
            texts.append("".join(t.text or "" for t in c.findall(qn("w:t"))))
    joined = "".join(texts)
    start = joined.find(old)
    if start < 0:
        return False
    end, acc, first, last = start + len(old), 0, None, None
    for j, t in enumerate(texts):
        if acc < end and acc + len(t) > start:
            first = j if first is None else first
            last = j
        acc += len(t)
    if first is None or first == last:
        return False
    return any(kids[i].tag != qn("w:r") for i in range(runs[first], runs[last]))


def edit(p, pairs, eid: str):
    """safe_tracked_replace with a per-replacement report; a miss is fatal.

    Falls back to retyping inside our own unaccepted insertion, which is the one case where
    a 'not-found' is expected rather than wrong.
    """
    for old, new in pairs:
        if spans_revision(p, old):
            raise AssertionError(
                f"{eid}: anchor {old[:60]!r} straddles an existing revision; "
                f"shorten it so it lies inside one uninterrupted run")
        (status,) = [s for _, s in safe_tracked_replace(p, [(old, new)])]
        if status == "not-found" and _retype_in_own_ins(p, old, new):
            status = "ok-in-ins"
        report.append((eid, status, old[:70]))
        if status not in ("ok", "ok-in-ins"):
            raise AssertionError(f"{eid}: {status} for {old[:80]!r}")


def delete_paragraph_smart(p) -> str:
    """Delete a paragraph, choosing the representation its content actually needs.

    A paragraph whose entire content is OUR OWN unaccepted <w:ins> was never in the original
    document, so the correct deletion is to drop the element outright - which is what Word
    does when you delete your own not-yet-accepted insertion, and it leaves Reject All
    untouched. word_rewrite_trackchanges.delete_paragraph only wraps DIRECT-CHILD runs, so on
    such a paragraph it is a silent no-op: the note stays fully visible in the Accept All view.
    Anything with original text is deleted the tracked way instead.
    """
    runs = [c for c in p if c.tag == qn("w:r")]
    inses = [c for c in p if c.tag == qn("w:ins")]
    foreign = [i for i in inses
               if not (i.get(qn("w:author")) or "").startswith(OUR_AUTHORS)]
    if foreign:
        raise AssertionError(
            "refusing to delete a paragraph carrying another author's unaccepted insertion")
    has_own_text = any((t.text or "").strip() for r in runs for t in r.findall(qn("w:t")))

    # Our own insertions were never in the original, so dropping them outright is both the
    # correct deletion and invisible to Reject All. Leaving them is NOT an option: a tracked
    # deletion only wraps direct-child runs, so an insertion left behind stays fully visible
    # in the Accept All view - which is how two stale legend stubs survived the first attempt.
    for i in inses:
        p.remove(i)

    if not has_own_text:
        p.getparent().remove(p)
        return "removed (was our own unaccepted insertion)"
    delete_paragraph(p)
    return f"deleted (tracked; {len(inses)} of our insertions dropped)" if inses \
        else "deleted (tracked)"


def replace_block(body, anchors: list[str], new_paragraphs: list[str], eid: str):
    """Delete a run of paragraphs (tracked) and insert replacements after the last one."""
    ps = [find_p(body, a) for a in anchors]           # resolve ALL refs before mutating
    for p, a in zip(ps, anchors):
        assert_no_citation(p, f"{eid}: {a[:50]}")
    tail = ps[-1]
    new_els = [ins_paragraph(t) for t in new_paragraphs]
    insert_after(tail, new_els)
    hows = [delete_paragraph_smart(p) for p in ps]
    report.append((eid, "ok", f"{len(ps)} paragraphs replaced by {len(new_els)}"
                              f" ({'; '.join(sorted(set(hows)))})"))


def ins_table(rows: list[list[str]], header_bold: bool = True):
    """A w:tbl whose every row is marked as a tracked insertion."""
    tbl = OxmlElement("w:tbl")
    pr = OxmlElement("w:tblPr")
    style = OxmlElement("w:tblStyle"); style.set(qn("w:val"), "TableGrid"); pr.append(style)
    wd = OxmlElement("w:tblW"); wd.set(qn("w:w"), "0"); wd.set(qn("w:type"), "auto"); pr.append(wd)
    borders = OxmlElement("w:tblBorders")
    for edge in ("top", "left", "bottom", "right", "insideH", "insideV"):
        b = OxmlElement(f"w:{edge}")
        b.set(qn("w:val"), "single"); b.set(qn("w:sz"), "4"); b.set(qn("w:color"), "auto")
        borders.append(b)
    pr.append(borders)
    tbl.append(pr)
    grid = OxmlElement("w:tblGrid")
    for _ in rows[0]:
        gc = OxmlElement("w:gridCol"); grid.append(gc)
    tbl.append(grid)
    for i, row in enumerate(rows):
        tr = OxmlElement("w:tr")
        trpr = OxmlElement("w:trPr")
        trpr.append(_rev_attrs(OxmlElement("w:ins")))          # the row itself is inserted
        tr.append(trpr)
        for cell in row:
            tc = OxmlElement("w:tc")
            tcpr = OxmlElement("w:tcPr")
            wd = OxmlElement("w:tcW"); wd.set(qn("w:w"), "0"); wd.set(qn("w:type"), "auto")
            tcpr.append(wd); tc.append(tcpr)
            p = OxmlElement("w:p")
            ppr = OxmlElement("w:pPr"); rpr = OxmlElement("w:rPr")
            rpr.append(_rev_attrs(OxmlElement("w:ins"))); ppr.append(rpr); p.append(ppr)
            p.append(wrap_ins(make_run(cell, bold=(i == 0 and header_bold))))
            tc.append(p)
            tr.append(tc)
        tbl.append(tr)
    return tbl


# ================================================================= the new prose

IFNA = [
    "The divergence-stratified analysis placed low-divergence LINE elements next to B, T and NK "
    "cell activation genes (Figure 5A). Those three GO terms share a single core gene set - "
    "IFNA10, IFNA16, IFNA17, IFNA21, IFNA4, IFNA6, IFNA7 and IFNW1 - whose TSS neighborhoods all "
    "fall inside one 219,363 bp interval on chromosome 9 (chr9:21,150,692-21,370,055), the "
    "interferon-alpha domain. Because one locus carries three of the class-level associations, we "
    "analyzed it separately and in full.",

    "The domain contains 175 annotated transposable elements (File S4, sheet IFNA_elements; "
    "Figure 8A). L1 is its largest family with 77 elements, 44 % of the domain, followed by Alu "
    "(33), L2 (15), MIR (13), hAT-Charlie (11), ERV1 (10) and ERVL-MaLR (6); the remaining 10 "
    "elements belong to seven further families. L1 density in the domain is 351 elements per Mb "
    "against 181.4 per Mb genome-wide, that is 565,459 L1 copies over 3,117 Mb, a 1.94-fold "
    "excess.",

    "The L1 elements of the domain are younger than the genome-wide L1 population. Their mean "
    "divergence is 135.7 against 188.2 for all L1 elements, and the per-gene average divergence "
    "of the LINE elements around each of the eight interferon TSS ranges from 95.0 for IFNA6 to "
    "161.7 for IFNA21 (File S1, sheet TSS_TE_intersections). The 77 elements span 36 distinct L1 "
    "subfamilies (File S4, sheet IFNA_subfamily_composition; Figure 8C), including the "
    "primate-specific L1PA2-L1PA8, L1PA10, L1PA14, L1PA15, L1P1, L1P3, L1P4e, L1P5, L1PB and "
    "L1PREC2 alongside older mammalian L1M clades. 38 of the 77 elements belong to a young "
    "primate-specific clade, against 133,450 of the 545,659 genome-wide L1 copies assignable to "
    "either the young L1HS and L1P* or the older L1M* clade (Fisher exact odds ratio 3.01, raw "
    "p = 3.2 x 10-6); HAL1 and X9_LINE elements, which belong to neither clade, are excluded from "
    "this test.",

    "The low mean divergence is not a by-product of the high L1 density of the domain, of its "
    "gene density, or of a few extreme elements (File S4, sheet IFNA_tests; Figure 8B). Against "
    "10,000 random 220 kb autosomal windows containing at least one L1, the observed mean of "
    "135.7 lies 2.28 standard deviations below the null mean of 189.4 (two-sided empirical "
    "p = 0.022, raw). Restricting the null to windows with at least 40 L1 elements, which matches "
    "the L1 density of the domain, strengthens the result (10,000 windows, z = -3.07, "
    "p = 0.0061), and restricting it to gene-dense windows with at least 10 annotated genes and "
    "at least one L1 strengthens it further (3,582 windows, z = -3.01, p = 0.0017), where 1 of "
    "the 3,582 null windows is as extreme as the observation. Removing individual elements does "
    "not change the conclusion (Figure 8D): excluding the single element annotated with "
    "divergence 0, implausible for its L1P3 clade and possibly a RepeatMasker artefact, raises "
    "the mean to 137.5, and removing the five youngest elements raises it to 143.8, both below "
    "every null mean.",

    "The domain is not a product of the new assembly: no base of it is newly resolved relative to "
    "hg38 and the whole interval is alignable to the previous reference (File S4, sheet "
    "assembly_bound), so this L1 cluster was visible before the T2T assembly. What the analysis "
    "establishes is positional. Young L1 copies are concentrated in a locus of innate immune "
    "genes to a degree that three matched null models reject, which identifies the domain as a "
    "candidate for functional follow-up. It does not establish that these elements regulate the "
    "interferon genes, which would require the expression and chromatin evidence that the present "
    "design does not carry.",
]

LU_HEAD_NEW = ("To place the current results next to the earlier repeat-based classification of "
               "genes, we tested the overlap between the current top-5 % gene sets and the three "
               "repeat-enriched mouse gene categories of Lu et al ")   # citation control follows

LU_TAIL_NEW = (". Their categories are defined on the mouse genome, so both sides were mapped to "
               "a common universe of human genes with a mouse ortholog through MGI homology "
               "classes (File S4, sheet prior_work_categories). The mapping is not uniform across "
               "their three categories: 2,042 of 2,439 low-complexity-repeat-enriched genes "
               "(83.7 %) and 1,694 of 2,041 SINE-enriched genes (83.0 %) carry an MGI ortholog, "
               "against 652 of 1,480 L1-enriched genes (44.1 %), so every comparison involving "
               "their L1 category rests on roughly half as many genes.")

LU_REST = [
    "The like-for-like pairs agree. Our Alu-enriched genes overlap their SINE-enriched category "
    "2.87-fold above chance, sharing 273 genes (odds ratio 4.02, FDR = 4.4 x 10-61), which is the "
    "strongest association in the matrix and means that a SINE-based categorization built on the "
    "mouse genome recovers substantially the same human genes as an Alu-based one built here. Our "
    "SINE class as a whole behaves the same way at 2.71-fold (270 genes, odds ratio 3.69, "
    "FDR = 6.1 x 10-55), as expected because Alu contributes most of the human SINE content. Our "
    "L1-enriched genes overlap their L1-enriched category 2.22-fold (49 genes, odds ratio 2.44, "
    "FDR = 6.3 x 10-7), a weaker but consistent agreement on the second element class.",

    "The cross pairs are depleted, which is the control the like-for-like agreement needs. Our "
    "SINE set against their L1 category is 0.13-fold (4 genes, FDR = 7.6 x 10-9) and our Alu set "
    "against their L1 category 0.41-fold (12 genes, FDR = 5.8 x 10-4), so the agreement above is "
    "specific to the matching element class rather than a general tendency of TE-rich gene sets "
    "to coincide.",

    "One counterpart pair diverges, and we report it rather than set it aside. Our MIR-enriched "
    "genes are depleted, not enriched, in their SINE-enriched category (45 genes, 0.46-fold, "
    "FDR = 2.8 x 10-9). MIR elements are ancient tRNA-derived SINEs shared across mammals, "
    "whereas the SINE content driving a mouse-based categorization is dominated by the young "
    "rodent-specific B1 and B2 families: the two are the same class but not the same elements, so "
    "the disagreement is a property of SINE evolution rather than of either analysis. The eight "
    "comparisons that carry these conclusions are collected in Table 3.",
]

LU_CLOSE = ("The full matrix of all 33 TE group by category pairs, with fold enrichment, raw and "
            "adjusted Fisher p-values, Jaccard and overlap coefficients, is given in File S4, "
            "sheet prior_work_overlap_matrix, and the individual shared genes in sheet "
            "prior_work_shared_genes.")

TABLE3_TITLE = ("Table 3. Overlap between the current top-5 % gene sets and the repeat-enriched "
                "mouse gene categories of Lu et al. Fold is the observed overlap divided by the "
                "overlap expected by chance in the shared human-mouse ortholog universe; "
                "p-values are Benjamini-Hochberg adjusted across all 33 pairs of the full matrix "
                "(File S4, sheet prior_work_overlap_matrix).")

TABLE3 = [
    ["Current gene set", "Lu et al. category", "Shared genes", "Fold over expected",
     "Odds ratio", "Adjusted Fisher p"],
    ["Alu (family)", "SINE-enriched", "273", "2.87", "4.02", "4.4*10-61"],
    ["SINE (class)", "SINE-enriched", "270", "2.71", "3.69", "6.1*10-55"],
    ["L1 (family)", "L1-enriched", "49", "2.22", "2.44", "6.3*10-7"],
    ["L1 (family)", "SINE-enriched", "48", "0.67", "0.63", "3.2*10-3"],
    ["Alu (family)", "L1-enriched", "12", "0.41", "0.39", "5.8*10-4"],
    ["SINE (class)", "L1-enriched", "4", "0.13", "0.12", "7.6*10-9"],
    ["MIR (family)", "SINE-enriched", "45", "0.46", "0.43", "2.8*10-9"],
    ["MIR (family)", "L1-enriched", "12", "0.40", "0.38", "5.0*10-4"],
]

SENS = [
    "We repeated the whole enrichment analysis at TSS neighborhoods of 5 kb and 20 kb, keeping "
    "the same TSS definition and a separate N = 500 permutation background for each window (File "
    "S5, sheets window_classes and window_families). The observed-to-random odds ratio of every "
    "class and family at all three windows is shown in Figure S11A: the ordering of TE groups is "
    "preserved, with Spearman rho between window pairs from 0.828 to 1.000 and a label-shuffling "
    "permutation test rejecting chance concordance in all six pairs, the largest p being 0.0090 "
    "(File S5, sheet window_concordance). The Bland-Altman comparison of the 44 families (Figure "
    "S11B) puts the mean difference between windows at 0.027 to 0.084 observed-to-random odds "
    "ratio units, so widening the window shifts the family distribution slightly upward without "
    "reordering it.",

    "Significance calls change for 11 of the 50 tested groups (File S5, sheet window_flips): the "
    "RC class, called Helitrons in Table 1, and ten families - ERVL, Gypsy, Helitron, I-Jockey, "
    "MULE-MuDR, PiggyBac, TcMar, TcMar-Mariner, hAT and tRNA-Deu. All eleven are small groups "
    "lying near the significance threshold and no other class flips, so the window size decides "
    "borderline calls for small groups and nothing else.",

    "The gene sets themselves are stable across windows (Figure S12; File S5, sheets "
    "geneset_stability and rank_stability, in which the RepeatMasker class names Retroposon and "
    "RC denote SVA and Helitrons). The overlap coefficient of the top-5 % gene sets between "
    "window pairs ranges from 0.285 to 0.931 and every hypergeometric test against chance is "
    "significant, the weakest being 3.5 x 10-205; the per-gene rankings agree with Kendall tau "
    "between 0.482 and 0.793. The two weakest cells are not the same one: the smallest gene-set "
    "overlap is that of the SVA class between 5 and 20 kb at 0.285, which is expected for the "
    "smallest class, whose 6,274 elements give sparse per-gene counts and heavy tying, whereas "
    "the lowest Kendall tau is that of the DNA class over the same window pair at 0.482.",

    "Widening the gene-set cut-off from the top and bottom 5 % to 10 % preserves a median of "
    "90.3 % of the significant GO terms per group and adds more, losing 144 terms and gaining 902 "
    "across all groups (Figure S11C; File S5, sheets percentile_summary and percentile_terms). Of "
    "the nine functional associations this paper reports, eight remain significant at the 10 % "
    "cut-off (File S5, sheet headline_by_condition). The exception is the association of "
    "hAT-Charlie with the MHC class I protein complex, which holds only in the 5 % set and "
    "survives one of the six window by cut-off conditions (Figure S11D); it is reported as the "
    "weakest of the family-level associations.",

    "Repeating the Gene Ontology analysis at every combination of the three windows and the two "
    "cut-offs gives an 18-cell grid (Figure S13A; File S5, sheet GO_grid_index), in which the "
    "window matters more than the cut-off. Widening the cut-off always finds more terms, as a "
    "larger foreground should, but widening the window does not, and not in the same direction at "
    "every level: the number of significant terms falls with window width for the class sets "
    "defined by element count, peaks at 10 kb for the divergence-stratified sets and rises for "
    "the families. Of the terms this paper reports, the fraction still significant in another "
    "cell is 0.440 at worst and 0.677 at median (Figure S13B; File S5, sheet "
    "GO_grid_preservation), and the per-group term counts of every other cell correlate with the "
    "published cell at a Spearman rho of at least 0.614, every label-shuffling permutation test "
    "giving p <= 0.022 (File S5, sheet GO_grid_concordance). Three of the nine reported claims "
    "survive all six conditions and all nine hold in the published one; the TE groups that keep "
    "at least one significant term in each cell are shown in Figure S13C. Because the grid varies "
    "the gene sets and not the permutation background, these differences are a property of which "
    "genes enter each set, and the enrichment odds ratios of Table 1 and Figure 2 are unchanged "
    "by it.",
]

S1_LEGEND_BODY = (
    " Technical controls on the permutation background and on element length. "
    "(A) Running mean of the random odds ratio for each TE class as permutations accumulate, as a "
    "fraction of its value at N = 500; the shaded band is +/- 1 %. (B) Running standard deviation "
    "of the same quantity, on the same scale. (C) Drift of the running mean for all 44 families, "
    "expressed in units of the final standard deviation, so that families with different absolute "
    "odds ratios are comparable. (D, E) Ridge plots for length distribution comparison between "
    "all TEs (blue) and those mapped on TSS neighborhoods (red), (D) for all classes together and "
    "(E) for individual classes; distributions were compared by two-sample Kolmogorov-Smirnov "
    "test. The p-values shown are raw, as each panel is a single test; panel (E) reports one test "
    "per class, six in total, of which five are significant."
)

S12_LEGEND_BODY = (
    " Stability of the gene sets across TSS window sizes. Overlap coefficient of the "
    "top 5 % gene sets and Kendall correlation of the gene rankings with 95 % bootstrap "
    "confidence intervals, between each pair of window sizes, for each TE class and for all TEs. "
    "The Kendall correlations are raw, as each is a single test. The RepeatMasker class names "
    "Retroposon and RC denote the SVA and Helitron classes."
)


# ================================================================= the edit programme

def apply_edits(doc):
    body = doc.element.body

    # ---- A. the three 'Daniil to Claude' comments -------------------------------------
    for i, key in enumerate(["Daniil to Claude: I rewrote this subsection",
                             "Daniil to Claude: Please rewrite this subsection citing individual",
                             "Daniil to Claude: The same request here"], start=1):
        p = find_p(body, key)
        assert_no_citation(p, f"A{i}")
        how = delete_paragraph_smart(p)
        report.append((f"A{i}", "ok", f"comment resolved and {how}"))

    # ---- B. interferon-alpha domain ---------------------------------------------------
    replace_block(
        body,
        ["As we shown in Figure 5A, LINE elements of low divergence",
         "We characterized LINEs in interferon alpha domain in a greater detail",
         "The 77 L1 elements in the domain span 36 distinct subfamilies",
         "We also compared the observed low divergence to the genomic background"],
        IFNA, "B")

    # subsection heading: it is a sentence, make it a title
    p = find_p(body, "The interferon-alpha domain in depth analysis.")
    edit(p, [("The interferon-alpha domain in depth analysis.",
              "The interferon-alpha domain in depth")], "B-head")

    # ---- C. prior-work comparison -----------------------------------------------------
    lu = find_p(body, "To place the current results results next to the earlier")
    if lu.find(".//" + qn("w:sdt")) is None:
        raise AssertionError("C: expected the Lu paragraph to carry its citation control")
    # Resolve BOTH anchors before mutating: the head edit consumes the text the tail
    # anchor is derived from. Each lies wholly on one side of the <w:sdt>; neither spans it.
    head_old = ("To place the current results results next to the earlier repeat-based "
                "classification of genes, we tested the overlap between the current top-5 % gene "
                "sets and the three published repeat-enriched gene categories of the landmark "
                "study by Lu et al ")
    tail_old = visible(lu).split("Lu et al ", 1)[1]
    tail_old = tail_old[tail_old.index("."):]                 # from the '.' after the citation
    edit(lu, [(head_old, LU_HEAD_NEW)], "C-head")
    edit(lu, [(tail_old, LU_TAIL_NEW)], "C-tail")

    new_els = [ins_paragraph(t) for t in LU_REST]
    new_els.append(ins_paragraph(TABLE3_TITLE))
    new_els.append(ins_table(TABLE3))
    new_els.append(ins_paragraph(LU_CLOSE))
    insert_after(lu, new_els)
    report.append(("C-body", "ok", f"{len(LU_REST)} paragraphs + Table 3 + closing pointer"))

    p = find_p(body, "Comparison with the previous repeat-based gene categorization.")
    edit(p, [("Comparison with the previous repeat-based gene categorization.",
              "Comparison with the previous repeat-based gene categorization")], "C-head2")

    # C-disc: the species difference is the dominant methodological term, and it also explains
    # why the assembly bound is not the interesting one here.
    p = find_p(body, "the length-bias correction by permutation")
    insert_after(p, [ins_paragraph(
        "The largest methodological difference between the two studies is not the window or the "
        "normalization but the species. Their categories were derived on the mouse genome and "
        "ours on the human T2T assembly, so the comparison is necessarily cross-species and has "
        "to pass through orthology, which caps the agreement that could be observed even if both "
        "analyses were otherwise identical: only genes with an MGI homology-class ortholog can "
        "be tested at all, and that is 83.7 % and 83.0 % of their low-complexity- and "
        "SINE-enriched categories but only 44.1 % of their L1-enriched one. The choice of human "
        "assembly is, for this particular comparison, close to irrelevant for the same reason - "
        "their gene categories were never derived from any human assembly, so the 6.70 % of "
        "newly resolved sequence bounded above cannot move the overlap at all except through the "
        "0.55 % of human genes that fall in it. What survives both the species gap and the "
        "design gap is therefore the more notable result: the gene sets still agree 2.87-fold "
        "for Alu against their SINE-enriched category and 2.22-fold for L1 against their "
        "L1-enriched one, and are depleted in every cross pair. Two independent designs, on two "
        "genomes, recover overlapping sets of human genes as the ones most associated with the "
        "same element classes.")])
    report.append(("C-disc", "ok", "species-difference paragraph added to Comparison with prior work"))

    # ---- D. sensitivity ---------------------------------------------------------------
    replace_block(
        body,
        ["In order to test robustness of the current analysis to TSS neighborhood size changes",
         "To further evaluate the gene sets robustness",
         "Repeating the Gene Ontology analysis at every combination"],
        SENS, "D")

    p = find_p(body, "Analysis sensitivity to the TSS window size and the gene-set cut-off.")
    edit(p, [("Analysis sensitivity to the TSS window size and the gene-set cut-off.",
              "Analysis sensitivity to the TSS window size and the gene-set cut-off")], "D-head")

    # ---- E. numbers that do not reproduce ---------------------------------------------
    # E1: MIR borderline band is 10 terms, not 6, and the values live nowhere in the package
    p = find_p(body, "6 GO terms with FDR between 0.05 and 0.1 were rejected as insignificant")
    edit(p, [
        ("6 GO terms with FDR between 0.05 and 0.1 were rejected as insignificant:  cellular "
         "response to cadmium ion and detoxification of copper ion (FDR = 0.078 for both), as are "
         "macrophage activation (FDR = 0.053), phosphatidylinositol-4,5-bisphosphate binding, "
         "negative regulation of response to cytokine stimulus (both FDR = 0.062) and the "
         "complement component C1q complex (FDR = 0.078).",
         "10 GO terms of MIR elements fell between FDR 0.05 and 0.1 and were rejected as "
         "insignificant, among them cellular response to cadmium ion, detoxification of copper "
         "ion, the complement component C1q complex and neuromuscular junction (FDR = 0.078 for "
         "all four), macrophage activation (FDR = 0.053) and actin binding, "
         "phosphatidylinositol-4,5-bisphosphate binding and negative regulation of response to "
         "cytokine stimulus (FDR = 0.062 for all three); the complete band is listed in File S3, "
         "sheet GO_borderline."),
    ], "E1")

    # E2: the LTR class has three terms in that band, not two
    p = find_p(body, "At the 0.05 FDR threshold the LINE class had insignificant term")
    edit(p, [
        ("At the 0.05 FDR threshold the LINE class had insignificant term of flavone metabolic "
         "process (FDR = 0.088) and the LTR had glutamatergic synapse and positive regulation of "
         "lipopolysaccharide-mediated signaling (both FDR = 0.086), that could be deemed "
         "significant under a more relaxed FDR threshold of 0.1 but are rejected here.",
         "At the 0.05 FDR threshold the LINE class had one insignificant term, flavone metabolic "
         "process (FDR = 0.088), and the LTR class three: nucleotide binding (FDR = 0.076), "
         "glutamatergic synapse and positive regulation of lipopolysaccharide-mediated signaling "
         "(FDR = 0.086 for both). All four would be significant under a more relaxed FDR "
         "threshold of 0.1 and are rejected here; every term in that band, for every group, is "
         "listed in File S3, sheet GO_borderline."),
    ], "E2")

    # E3: Figure 8C legend - 545,659 is the two-clade total, not the L1 total
    p = find_p(body, "Young primate-specific copies are enriched in the domain")
    edit(p, [("(38 of 77 against 133,450 of 545,659; Fisher exact odds ratio 3.01, raw "
              "p = 3.2 x 10-6)",
              "(38 of 77 against 133,450 of the 545,659 genome-wide L1 copies assignable to "
              "either clade; Fisher exact odds ratio 3.01, raw p = 3.2 x 10-6)")], "E3")

    # E4: the assembly bound has a source; give it one
    p = find_p(body, "The interferon-alpha domain contains no newly resolved sequence at all")
    edit(p, [("The interferon-alpha domain contains no newly resolved sequence at all and is "
              "fully alignable to hg38.",
              "The interferon-alpha domain contains no newly resolved sequence and is fully "
              "alignable to hg38 (File S4, sheet assembly_bound).")], "E4")

    # ---- F. supplementary figure renumbering, body ------------------------------------
    #      old S8->S6, S9->S7, S6->S8, S7->S9; S1D/1E written out in full
    p = find_p(body, "either for all TEs or for their individual classes (Figure S1D, 1E)")
    edit(p, [("S1D, 1E", "S1D, S1E")], "F1")

    # both legends share the sentence that carries the citation, so anchor on the legend head
    p = find_p(body, "Figure 4. Functional analysis of genes whose TSS are enriched or depleted")
    edit(p, [("Figure S8", "Figure S6")], "F2")

    p = find_p(body, "Figure 5. Functional analysis of genes whose TSS are enriched with TEs")
    edit(p, [("Figure S9", "Figure S7")], "F3")

    p = find_p(body, "then visualized the TE count distributions in the TSS vicinity")
    edit(p, [("Figure S6", "Figure S8")], "F4")

    p = find_p(body, "and their overlapping’s is shown in")
    edit(p, [("Figure S7A", "Figure S9A")], "F5a")
    edit(p, [("Figure S7A", "Figure S9A")], "F5b")

    p = find_p(body, "We also compared family-level GO terms with the respective class-level")
    edit(p, [("Figure S7B", "Figure S9B")], "F6")

    p = find_p(body, "unfiltered representation in Figure S7C")
    edit(p, [("Figure S7C", "Figure S9C"), ("Figure S7B", "Figure S9B")], "F7")

    # ---- F. supplementary legends: renumber, rewrite two, reorder ---------------------
    # A legend is the paragraph that STARTS with its own number. Substring matching would
    # collide with the body citations F2-F7 have just created ("...the full network is
    # Figure S6." reads as "Figure S6. " too).
    def find_legend(n: int):
        hits = [p for p in body.iter(qn("w:p"))
                if visible(p).strip().startswith(f"Figure S{n}.")]
        if len(hits) != 1:
            raise LookupError(f"legend S{n}: {len(hits)} paragraphs start with it")
        return hits[0]

    legends = {old: find_legend(old) for old in range(1, 14)}
    OLD2NEW = {1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 8, 7: 9, 8: 6, 9: 7, 10: 10,
               11: 13, 12: 11, 13: 12}

    # A legend whose NUMBER changes but whose text does not is patched in place: the number
    # token sits in its own <w:ins> from an earlier pass, so a one-token replacement is enough.
    # The two legends whose CONTENT changes (S1 gains the convergence panels, old S13 loses
    # them) are replaced wholesale instead: their text is split across a direct-child run and
    # an <w:ins>, so no single anchor spans it. Neither carries a citation control.
    for old, new in OLD2NEW.items():
        if old in (1, 13):
            continue
        if old != new:
            edit(legends[old], [(f"Figure S{old}.", f"Figure S{new}.")], f"F-leg{old}")

    for old, full_new in ((1, "Figure S1." + S1_LEGEND_BODY),
                          (13, "Figure S12." + S12_LEGEND_BODY)):
        p = legends[old]
        assert_no_citation(p, f"F-leg{old}")
        new_p = ins_paragraph(full_new)
        insert_after(p, [new_p])
        how = delete_paragraph_smart(p)
        legends[old] = new_p                       # the reorder below must move the new one
        report.append((f"F-leg{old}", "ok",
                       f"legend replaced in full -> S{OLD2NEW[old]} ({how})"))

    # physical reorder into S1..S13 (a structural move, not a Word tracked move)
    anchor = legends[1]
    for new in range(2, 14):
        old = [o for o, n in OLD2NEW.items() if n == new][0]
        el = legends[old]
        el.getparent().remove(el)
        anchor.addnext(el)
        anchor = el
    report.append(("F-order", "ok", "13 supplementary legends reordered to S1..S13"))

    # ---- G. process scaffolding and proofless statements ------------------------------
    notes = [p for p in body.iter(qn("w:p")) if "[EDITOR NOTE]" in visible(p)]
    how = []
    for p in notes:
        assert_no_citation(p, "G-note")
        how.append(delete_paragraph_smart(p))
    from collections import Counter as _C
    report.append(("G-notes", "ok",
                   f"{len(notes)} editor notes removed: " + "; ".join(
                       f"{n} {k}" for k, n in _C(how).items())))

    p = find_p(body, "Gemini 2.5 PRO")
    edit(p, [("and Claude Code Opus 5.0 (ref) was", "and Claude Code (Opus 5) were")], "G1")

    # (The 'It should be noted that the 95.0-161.7 range' caveat lived inside the first
    #  interferon paragraph, which block B has already deleted; the rewritten passage 3
    #  states the same distinction without the throat-clearing, so nothing is left to do.)

    p = find_p(body, "Alu and MIR insertions-adjacent GO terms were remarkably different")
    edit(p, [("were remarkably different", "differed")], "G3")

    # (G4: 'which is what justifies reporting the background at 500' lived in the old S13
    #  legend, which F-leg13 replaces in full; the new S1 legend states the convergence
    #  panels without that self-assessment, so there is nothing left to strip.)
    if any("what justifies reporting the background" in visible(p)
           for p in body.iter(qn("w:p"))):
        raise AssertionError("G4: the proofless clause survived the legend rewrite")
    report.append(("G4", "ok", "proofless clause gone with the old S13 legend"))

    p = find_p(body, "Sensitivity and robustness analyses. Sheets window_classes")
    edit(p, [("Sensitivity and robustness analyses.",
              "Sensitivity analyses across TSS window size and gene-set cut-off.")], "G5")

    # G6: the last surviving "what changed since the previous draft" phrase
    p = find_p(body, "serine-type endopeptidase inhibitor activity (all FDR = 0.086)")
    edit(p, [("at the tightened threshold", "at FDR 0.05")], "G6")

    # G7: an unquantified intensifier in the prose this pass inserted
    p = find_p(body, "recovers substantially the same human genes")
    edit(p, [("recovers substantially the same human genes",
              "recovers the same human genes")], "G7")

    # G8: keep the supplementary titles on one term - "sensitivity", not "robustness"
    p = find_p(body, "Robustness of the Gene Ontology results across window size")
    edit(p, [("Robustness of the Gene Ontology results",
              "Sensitivity of the Gene Ontology results")], "G8")

    # G9: two rounded ranges, replaced by the per-class values derived from the published
    #     gene sets (File S2, sheet by_divergence, joined to the per-TSS class divergences)
    # The sentence already cites File S2, sheet by_divergence, so only the two rounded ranges
    # need replacing; each numeric token lies inside a single run, unlike the whole clause.
    p = find_p(body, "The upper limits for lowest divergence were at the level of 10-17%")
    edit(p, [
        ("at the level of 10-17% depending on the TE class",
         "9.8 % for LTR, 11.1 % for SINE, 14.3 % for DNA and 16.8 % for LINE"),
        ("at the level of 28-31%",
         "27.9 % for SINE, 29.6 % for LTR, 31.1 % for DNA and 31.4 % for LINE"),
    ], "G9")

    # G10: the last two British spellings the mechanical pass could not reach
    for anchor, pair, eid in (
            ("Visualisation of up to five GO terms per family",
             ("Visualisation", "Visualization"), "G10a"),
            ("behaviour", ("behaviour", "behavior"), "G10b")):
        try:
            edit(find_p(body, anchor), [pair], eid)
        except LookupError:
            report.append((eid, "already-fixed", anchor[:60] + " (caught by the H-us pass)"))

    # ---- H. house style ---------------------------------------------------------------
    us = {"neighbourhood": "neighborhood", "neighbourhoods": "neighborhoods",
          "characterised": "characterized", "signalling": "signaling",
          "randomised": "randomized", "co-localised": "co-localized",
          "localised": "localized", "normalised": "normalized",
          "normalising": "normalizing", "normalisation": "normalization",
          "localisation": "localization", "visualisation": "visualization",
          "visualised": "visualized", "organised": "organized",
          "analysed": "analyzed", "coloured": "colored", "colour": "color",
          "colours": "colors", "behaviour": "behavior", "categorisation": "categorization",
          "labelled": "labeled", "summarises": "summarizes"}
    n_us = 0
    n_us_skipped: list = []
    # Longest key first, and recount after each pass: 'colour' is a substring of 'coloured'
    # and 'neighbourhood' of 'neighbourhoods', so a naive count queues phantom replacements
    # that then report not-found.
    us_keys = sorted(us, key=len, reverse=True)
    for p in list(body.iter(qn("w:p"))):
        txt = visible(p)
        if not txt or "[EDITOR NOTE]" in txt:
            continue
        pairs = []
        for brit in us_keys:
            n = txt.count(brit)
            if n:
                pairs.extend([(brit, us[brit])] * n)
                txt = txt.replace(brit, us[brit])
        if not pairs:
            continue
        for old, new in pairs:
            (status,) = [s for _, s in safe_tracked_replace(p, [(old, new)])]
            if status == "not-found":
                status = "ok-in-ins" if _retype_in_own_ins(p, old, new) else status
            if status in ("ok", "ok-in-ins"):
                n_us += 1
            else:
                n_us_skipped.append((old, status))
    report.append(("H-us", "ok", f"{n_us} British spellings changed to US"
                                 + (f"; {len(n_us_skipped)} skipped" if n_us_skipped else "")))
    for old, status in n_us_skipped:
        report.append(("H-us", status, old))

    # Unicode superscripts break journal typesetting and make the number unfindable by a
    # plain-text search for '10-91', which is what the verification protocol relies on.
    # Flatten the whole token so the change is a tracked edit in original text and a retype
    # inside our own insertions - never an untracked edit to text the author wrote.
    import re as _re
    SUPMAP = str.maketrans("⁻⁰¹²³⁴⁵⁶⁷⁸⁹", "-0123456789")
    SUPTOK = _re.compile(r"[⁻⁰¹²³⁴⁵⁶⁷⁸⁹]+")
    n_sup, n_sup_skipped = 0, []
    for p in list(body.iter(qn("w:p"))):
        txt = visible(p)
        if not txt or not SUPTOK.search(txt):
            continue
        seen = set()
        pairs = []
        for m in SUPTOK.finditer(txt):
            tok = m.group(0)
            if tok in seen:
                continue
            seen.add(tok)
            pairs.extend([(tok, tok.translate(SUPMAP))] * txt.count(tok))
        for old, new in pairs:
            (status,) = [s for _, s in safe_tracked_replace(p, [(old, new)])]
            if status == "not-found":
                status = "ok-in-ins" if _retype_in_own_ins(p, old, new) else status
            if status in ("ok", "ok-in-ins"):
                n_sup += 1
            else:
                n_sup_skipped.append((old, status))
    report.append(("H-sup", "ok", f"{n_sup} Unicode superscript tokens flattened"
                                  + (f"; {len(n_sup_skipped)} skipped" if n_sup_skipped else "")))
    for old, status in n_sup_skipped:
        report.append(("H-sup", status, repr(old)))

    return doc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", action="store_true", help="print the last run's report only")
    a = ap.parse_args()
    if a.report:
        print(DST.exists() and f"{DST} exists" or "not built yet")
        return

    if not SRC.exists():
        raise SystemExit(f"input not found: {SRC}")
    set_revision_identity(AUTHOR, DATE)

    shutil.copy2(SRC, DST)
    doc = Document(str(DST))
    apply_edits(doc)
    doc.save(str(DST))

    bad = [r for r in report if r[1] not in ("ok", "ok-in-ins", "already-fixed")]
    width = max(len(r[0]) for r in report)
    for eid, status, detail in report:
        print(f"  {eid:<{width}}  {status:<10} {detail}")
    print(f"\n  {len(report)} operations, {len(bad)} not ok")
    print(f"  -> {DST}")
    if bad:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
