#!/usr/bin/env python
"""19_review_edits_260810.py — the final scientific review, applied as tracked changes.

Input :  Revised_manuscript/T2T_genes_article_G3_revision_260809.docx   (read-only from here on)
Output:  Revised_manuscript/T2T_genes_article_G3_revision_260810.docx

A deterministic rebuild from its input, like 16_ and 18_, not an accumulation of in-place
edits, so it may be re-run at will.

WHY THIS PASS EXISTS
  The figures were re-exported to current_figures_260810/ with the Figma renaming of
  Figures_renaming_260809.md applied, and the supplementary package is supplementary_260809/.
  Everything the manuscript cites was re-checked against those two deliverables and against
  the source tables. Twelve numbers did not reproduce; four structural defects were found that
  only show up in the Accept-All view.

WHAT IT DOES
  A. Corrects twelve numbers that do not reproduce from the source tables (P5). The largest is
     the class-level summary sentence, which quotes the ERVK family row for "all TEs".
  B. Corrects five GO term lists in the divergence subsection that name terms which are not
     significant at FDR 0.05, or do not exist (P5).
  C. Removes the duplicate Table 1 / Table 2 pair after Funding, and finishes the tracked
     deletion of the already-struck 7x11 table in Results so Accept All leaves no empty grid.
  D. Completes the truncated Ethical statement and removes its misplaced duplicate (P7).
  E. Repoints three stale cross-references left by the supplementary renumbering (P3).
  F. Adds the two missing body citations of Figures S6 and S7 so every object is cited from
     prose, not only from another figure's legend (P3).
  G. Deletes the remaining proofless and revision-history statements (P1, P6).
  H. House style: US spelling, one multiplication sign (P9).

CITATION SAFETY
  Five touched paragraphs carry Mendeley content controls. No control is deleted and no
  paragraph is rebuilt from its concatenated text: every anchor lies wholly inside runs that
  are direct children of the paragraph, which is all safe_tracked_replace touches. The
  <w:sdt> and MENDELEY_CITATION counts are asserted unchanged at the end. Reference NUMBERS
  are never touched - Mendeley owns them.

WHAT IT DELIBERATELY DOES NOT DO
  Orphan citation-only paragraphs, the [ZENODO DOI] placeholder, the missing citations in the
  Discussion and the two artwork defects in Figure 8 and Figure 3 are reported, not edited.
  See review_report_260810.md beside this script.

Run:  ~/venvs/collagen_3_11/bin/python revision_G3/19_review_edits_260810.py
      ~/venvs/collagen_3_11/bin/python revision_G3/19_review_edits_260810.py --report
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
import word_rewrite_trackchanges as wrt                      # noqa: E402
from word_rewrite_trackchanges import (                      # noqa: E402
    _rev_attrs, delete_paragraph, set_revision_identity,
)

HERE = Path(__file__).resolve().parent
SRC = HERE / "Revised_manuscript" / "T2T_genes_article_G3_revision_260809.docx"
DST = HERE / "Revised_manuscript" / "T2T_genes_article_G3_revision_260810.docx"

AUTHOR = "Claude (scientific review 2026-08-10)"
DATE = "2026-08-10T00:00:00Z"

OUR_AUTHORS = ("Claude",)

report: list[tuple[str, str, str]] = []


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
    hits = [p for p in body.iter(qn("w:p")) if needle in visible(p)]
    if not hits:
        raise LookupError(f"anchor not found: {needle!r}")
    if must_be_unique and len(hits) > 1:
        raise LookupError(f"anchor matches {len(hits)} paragraphs, lengthen it: {needle!r}")
    return hits[0]


def n_sdt(el) -> int:
    return len(el.findall(".//" + qn("w:sdt")))


def _retype_in_own_ins(p, old: str, new: str) -> bool:
    """Edit text living inside an unaccepted <w:ins> of one of our own earlier passes.

    safe_tracked_replace cannot see it: those runs are children of the w:ins, not of the
    paragraph. Retyping mirrors what Word does when you edit your own not-yet-accepted
    insertion, and Reject All still restores the original because the whole insertion is
    dropped. Never done inside another author's revision.
    """
    for ins in p.findall(qn("w:ins")):
        if not (ins.get(qn("w:author")) or "").startswith(OUR_AUTHORS):
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

    safe_tracked_replace matches across DIRECT-CHILD runs only and rebuilds the span from the
    first matched run, so a revision element sitting between the first and last matched run
    survives but moves - Reject All then returns an already-deleted fragment in the wrong
    place. Refuse such an anchor; shorten it until it lies inside one uninterrupted stretch.
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

    The <w:sdt> count of the paragraph is asserted unchanged around every replacement, which
    is the mechanical proof that no Mendeley control was disturbed.
    """
    for old, new in pairs:
        if spans_revision(p, old):
            raise AssertionError(
                f"{eid}: anchor {old[:60]!r} straddles an existing revision; "
                f"shorten it so it lies inside one uninterrupted run")
        before = n_sdt(p)
        (status,) = [s for _, s in safe_tracked_replace(p, [(old, new)])]
        if status == "not-found" and _retype_in_own_ins(p, old, new):
            status = "ok-in-ins"
        if n_sdt(p) != before:
            raise AssertionError(f"{eid}: citation control count changed on {old[:50]!r}")
        report.append((eid, status, old[:74]))
        if status not in ("ok", "ok-in-ins"):
            raise AssertionError(f"{eid}: {status} for {old[:80]!r}")


def drop_own_insertions(p, eid: str) -> None:
    """Remove our own unaccepted <w:ins> children, leaving any original / deleted content.

    Used on the two duplicate table captions, which are an original caption already wrapped
    in <w:del> plus a replacement we inserted. Dropping only the insertion leaves a paragraph
    that Accept All removes and Reject All restores exactly.
    """
    dropped = 0
    for ins in list(p.findall(qn("w:ins"))):
        author = ins.get(qn("w:author")) or ""
        if not author.startswith(OUR_AUTHORS):
            raise AssertionError(f"{eid}: refusing to drop another author's insertion")
        p.remove(ins)
        dropped += 1
    report.append((eid, "ok", f"dropped {dropped} of our insertions from a caption"))


def remove_inserted_table(tbl, eid: str) -> None:
    """Remove a table every one of whose rows is marked as our own tracked insertion."""
    rows = tbl.findall(qn("w:tr"))
    for tr in rows:
        trpr = tr.find(qn("w:trPr"))
        ins = trpr.find(qn("w:ins")) if trpr is not None else None
        if ins is None:
            raise AssertionError(f"{eid}: table row is not an insertion, refusing to remove")
        if not (ins.get(qn("w:author")) or "").startswith(OUR_AUTHORS):
            raise AssertionError(f"{eid}: table row inserted by another author")
    tbl.getparent().remove(tbl)
    report.append((eid, "ok", f"removed an inserted table of {len(rows)} rows"))


def mark_rows_deleted(tbl, eid: str) -> None:
    """Finish a half-done table deletion: strike the ROWS, not only their text.

    A table whose every <w:t> already sits inside a <w:del> still leaves an empty grid after
    Accept All, because a row is only removed when its <w:trPr> carries a <w:del>. This adds
    that marker, which is what Word writes when you delete table rows with tracking on.
    """
    rows = tbl.findall(qn("w:tr"))
    live = [t for t in tbl.iter(qn("w:t")) if _alive(t)]
    if live:
        raise AssertionError(f"{eid}: table still has {len(live)} undeleted text runs")
    for tr in rows:
        trpr = tr.find(qn("w:trPr"))
        if trpr is None:
            trpr = OxmlElement("w:trPr")
            tr.insert(0, trpr)
        if trpr.find(qn("w:del")) is None:
            trpr.append(_rev_attrs(OxmlElement("w:del")))
    report.append((eid, "ok", f"marked {len(rows)} already-struck rows as deleted"))


def _alive(t) -> bool:
    anc = t.getparent()
    while anc is not None and anc.tag != qn("w:p"):
        if anc.tag == qn("w:del"):
            return False
        anc = anc.getparent()
    return True


def delete_paragraph_smart(p, eid: str) -> None:
    """Delete a paragraph in whichever representation its content needs."""
    runs = [c for c in p if c.tag == qn("w:r")]
    inses = [c for c in p if c.tag == qn("w:ins")]
    foreign = [i for i in inses if not (i.get(qn("w:author")) or "").startswith(OUR_AUTHORS)]
    if foreign:
        raise AssertionError(f"{eid}: paragraph carries another author's insertion")
    if n_sdt(p):
        raise AssertionError(f"{eid}: refusing to delete a paragraph holding a citation")
    has_own_text = any((t.text or "").strip() for r in runs for t in r.findall(qn("w:t")))
    for i in inses:
        p.remove(i)
    if not has_own_text:
        p.getparent().remove(p)
        report.append((eid, "ok", "removed (was our own unaccepted insertion)"))
    else:
        delete_paragraph(p)
        report.append((eid, "ok", "deleted (tracked)"))


# ----------------------------------------------------------------- the edits

def apply_edits(doc) -> None:
    body = doc.element.body

    # ============================================================ A. numbers (P5)

    # A1  The class-level summary quotes the ERVK FAMILY row for "all TEs".
    #     ERVK: observed 1.9509 / random 1.7782 / fold 1.0971  (enrichment_families_with_random.csv)
    #     all TEs: 582,540 of 3,709,429 elements in windows against a mean 542,973.1 over the
    #     500-seed background -> observed 1.95, random 1.79, enrichment 1.086. The fold change
    #     is a ratio of count ratios and so does not depend on the window-geometry constant.
    edit(find_p(body, "In general, the observed OR for all TEs"), [
        ("was 1.94 (1.78 the random one)", "was 1.95 (1.79 the random one)"),
        ("an enrichment by a factor of 1.097", "an enrichment by a factor of 1.086"),
    ], "A1 all-TE enrichment")

    # A2  Count thresholds of the top-5 % class gene sets, read off the realised sets.
    #     LINE: 1,257 genes have >= 10 elements, so the 1,436th gene carries 9, not 10.
    #     TE top: 1,279 genes have >= 27, so the boundary count is 26.
    edit(find_p(body, "There were 1436 genes for each of the major classes"), [
        ("genes having at least 5, 10, 18 and 5 elements",
         "genes having at least 5, 9, 18 and 5 elements"),
        ("(starting with 27 till 42 elements per TSS)",
         "(starting with 26 till 42 elements per TSS)"),
    ], "A2 gene-set thresholds")

    # A3  Maximum copies per TSS: L2 is 12 and hAT-Charlie is 10, not 17 and 19; and the
    #     families are stated to be ordered by TSS count, so L2 (30,376) precedes L1 (30,184).
    edit(find_p(body, "we ordered TE families by number of TSS having at least one"), [
        ("30184 (78.0%, again up to 19 copies) had L1 elements, "
         "30376 (78.5%, again up to 17 copies) had L2 elements and "
         "18731 (48.4%, again up to 19 copies) had hAT-Charlie elements",
         "30376 (78.5%, up to 12 copies) had L2 elements, "
         "30184 (78.0%, again up to 19 copies) had L1 elements and "
         "18731 (48.4%, up to 10 copies) had hAT-Charlie elements"),
    ], "A3 family TSS counts")

    # A4  Two functional groups tie at two terms, not one (Figure 6B recomputation).
    edit(find_p(body, "Protein ubiquitination and degradation was the next least present"), [
        ("Protein ubiquitination and degradation was the next least present group, "
         "found in 2 cases.",
         "Protein ubiquitination and degradation and embryogenesis were the next least "
         "present groups, with 2 terms each."),
    ], "A4 functional-group counts")

    # A5  hAT-hAT19: 314.0 against a family mean of 226.2, i.e. 1.39-FOLD.
    edit(find_p(body, "and this single element had divergence"), [
        ("had divergence 1.39 higher than the average all elements in a family",
         "had divergence 1.39-fold higher than the family average"),
    ], "A5 hAT-hAT19 ratio")

    # A6  1,140 is a count of TSS; the 1,140 TSS carry 962 distinct genes.
    edit(find_p(body, "due to the low number of TSS with SVA elements"), [
        ("(1,140 genes, 2.9%)", "(1,140 TSS, 2.9%)"),
    ], "A6 SVA TSS count")

    # A7  TcMar-Tc2 is the one family in the depleted list with no class given.
    edit(find_p(body, "In contrast, 9 families were significantly depleted"), [
        ("TcMar-Tc2 (0.819)", "TcMar-Tc2 (DNA, 0.819)"),
    ], "A7 TcMar-Tc2 class")

    # ============================================================ B. GO term lists (P5)
    #
    # Every list below named terms that fall in the 0.05 < FDR <= 0.1 band, which this paper
    # rejects, or that the 500-gene cap keeps out of Figure 5A, or that do not exist at all.
    # Recomputed from output/GO_classes_divergence_fdr005.csv with the same <= 500-gene filter
    # the panel uses.
    div = find_p(body, "We visualized up to ten GO terms for each TE-divergence group")
    edit(div, [
        # B1  TE_all / lowest. "spermatogenesis" is not a term anywhere; "establishment of
        #     spindle localization" is FDR 0.071 and "fatty acid catabolic process" 0.080.
        ("returning terms about rRNA binding, spermatogenesis, mitotic spindle localization, "
         "subcortical maternal complex and fatty acids catabolism.",
         "returning the ribonucleoprotein complex, rRNA binding, the subcortical maternal "
         "complex and establishment of the Sertoli cell barrier."),
        # B2  GO:0016712 is the term itself; "flavin-dependent oxidoreductase activity" is a
        #     paraphrase a reader cannot find in File S3.
        ("a single term of flavin-dependent oxidoreductase activity (FDR = 0.039)",
         "a single oxidoreductase activity term acting on paired donors with reduced flavin "
         "as one donor (GO:0016712, FDR = 0.039)"),
        # B3  Typo, and all TEs of highest divergence carry no olfactory or smell term at
        #     either threshold - only LTRs of highest divergence do.
        ("Groups of high divergence (all TEs, LINEs, SINEs and LTRs) we co-clustered, with "
         "all TEs and LTRs of high divergence sharing olfactory receptors with LINEs of "
         "lowest divergence.",
         "Groups of high divergence (all TEs, LINEs, SINEs and LTRs) were co-clustered, and "
         "LTRs of high divergence shared olfactory receptor terms with LINEs of lowest "
         "divergence."),
        # B4  LINE / highest at FDR 0.05. cell-cell adhesion (0.066), cell population
        #     proliferation (0.052) and very-low-density lipoprotein particle (0.056) are all
        #     in the rejected band.
        ("LINEs of highest divergence were inserted near genes of voltage-gated potassium "
         "channel complex, cell adhesion, differentiation and proliferation, as well as "
         "genes of very-low-density lipoprotein components.",
         "LINEs of highest divergence were inserted near genes of potassium ion transport "
         "and the voltage-gated potassium channel complex, of protein kinase binding, of "
         "ventricular cardiac muscle cell differentiation and of the extracellular matrix "
         "structural constituent conferring tensile strength."),
        # B5  98 shared terms, so the count replaces "dozens".
        ("shared dozens of GO terms connected with nervous system",
         "shared 98 GO terms connected with nervous system"),
        # B6  The three most significant shared terms under the panel's own 500-gene cap are
        #     postsynaptic membrane, chemical synaptic transmission and axon guidance.
        #     "synaptic chemical transition" is not a GO term, and nervous system development
        #     has 642 annotated genes so the panel excludes it.
        ("Top 3 GO terms with the most significant enrichment were nervous system "
         "development, postsynaptic membrane and synaptic chemical transition",
         "The three most significant of these were postsynaptic membrane, chemical synaptic "
         "transmission and axon guidance"),
        # B7  ... and the range therefore closes at 10-10, not 10-11. The superscript run is
        #     edited on its own so the exponent keeps its formatting.
        ("-11", "-10"),
    ], "B divergence GO lists")

    # B8  SINEs of lowest divergence do have two significant terms; both are above the
    #     500-gene cap, which is why Figure 5A shows none.
    edit(find_p(body, "Finally, DNA elements of highest and lowest divergence"), [
        ("Finally, DNA elements of highest and lowest divergence, and SINEs of lowest "
         "divergence, revealed no significant GO enrichments reflecting",
         "Finally, DNA elements of highest and lowest divergence revealed no significant GO "
         "enrichments, and SINEs of lowest divergence only protein binding and extracellular "
         "exosome, two terms too general to be informative, reflecting"),
    ], "B8 SINE lowest divergence")

    # ============================================================ E. stale cross-references (P3)

    # E1  The SVA class-level association is the count-based one, Figure 4A, not Figure 5A.
    edit(find_p(body, "Finally, SVA elements, as previously for the class-level analysis"), [
        ("as previously for the class-level analysis (Figure 5A)",
         "as previously for the class-level analysis (Figure 4A)"),
    ], "E1 SVA pointer")

    # E2  The renumbering moved the first full network from S8 to S6, so "the same filters as
    #     Figure S8" now points at the log-scaled TSS distributions, which have no filters.
    edit(find_p(body, "Figure S7. Complete connection map of GO terms for the divergence"), [
        ("the same edge and term-size filters as Figure S8",
         "the same edge and term-size filters as Figure S6"),
    ], "E2 S7 legend cross-ref")
    edit(find_p(body, "Figure S10. Complete connection map of GO terms per TE family"), [
        ("the same edge and term-size filters as Figure S8",
         "the same edge and term-size filters as Figure S6"),
    ], "E3 S10 legend cross-ref")

    # ============================================================ F. missing body citations (P3)
    #
    # S6 and S7 were cited only from the legends of Figures 4 and 5. Their position in the
    # sequence is correct either way, but a supplementary figure that prose never mentions is
    # a numbering audit failure and easy for a reader to miss.
    edit(find_p(body, "The integrative network visualization of up to ten"), [
        ("per remaining group (Figure 4A) showed",
         "per remaining group (Figure 4A; the complete network is Figure S6) showed"),
    ], "F1 cite S6 in prose")
    edit(div, [
        ("to avoid overly general classification (Figure 5A)",
         "to avoid overly general classification (Figure 5A; the complete network is "
         "Figure S7)"),
    ], "F2 cite S7 in prose")

    # ============================================================ G. deletions (P1, P6)

    # G1  "retaining / losing ... at FDR 0.05" tells the reader about a previous threshold.
    #     "at FDR 0.05" is an unaccepted insertion of the 2026-08-09 pass and is blanked
    #     rather than struck: striking our own not-yet-accepted text would leave a deletion
    #     of text the original document never had.
    crp = find_p(body, "CR1 elements had surprising connections")
    edit(crp, [
        ("retaining synapse", "with synapse"),
        (" but losing regulation of postsynaptic", ", while regulation of postsynaptic"),
        ("inhibitor activity (all FDR = 0.086) ",
         "inhibitor activity are rejected (all FDR = 0.086)"),
        ("is retained at the L1 family level", "is significant at the L1 family level"),
    ], "G1 CR1 threshold history")
    if not _retype_in_own_ins(crp, "at FDR 0.05", ""):
        raise AssertionError("G1: could not blank the 'at FDR 0.05' insertion")
    report.append(("G1 CR1 threshold history", "ok-in-ins", "blanked 'at FDR 0.05'"))

    # G2  The roadmap listed the subsections in an order the Discussion does not follow, and
    #     ended on a "makes it interpretable" claim about the paper rather than about the data.
    edit(find_p(body, "What follows is organized as the principal findings"), [
        ("the specific hypotheses it generates, its connection with cancer, its use as a "
         "null model for epigenomic association studies, and the mechanistic framework that "
         "makes the positional pattern interpretable.",
         "the specific hypotheses it generates, the mechanistic framework for the positional "
         "pattern, its connection with cancer and its use as a null model for epigenomic "
         "association studies."),
    ], "G2 Discussion roadmap order")

    # G3  The paper does use arms-race vocabulary for its own hypothesis - in the title and
    #     the abstract - so the sentence as written is not true of it.
    edit(find_p(body, "the causal vocabulary of an evolutionary arms race"), [
        ("and the causal vocabulary of an evolutionary arms race is used in this paper only "
         "when reporting other people's conclusions.",
         "and the arms-race reading of the interferon-alpha domain is offered as a hypothesis "
         "rather than as a finding."),
    ], "G3 arms-race vocabulary claim")

    # G4  "they are the framework within which the positional pattern becomes interpretable"
    #     is a claim about the paper, and Figure 9 is cited twice in two adjacent sentences.
    edit(find_p(body, "The enrichment and depletion values reported above"), [
        ("which the schematic in Figure 9 summarizes and which the following four mechanisms "
         "describe. None of them is tested by this design; they are the framework within "
         "which the positional pattern becomes interpretable.",
         "which the schematic in Figure 9 summarizes. None of the four mechanisms below is "
         "tested by this design."),
    ], "G4 mechanistic framework claim")
    edit(find_p(body, "The degree of enrichment or deficiency of TE groups"), [
        ("can be determined by the following factors (Figure 9):",
         "can be determined by the following factors:"),
    ], "G5 duplicate Figure 9 citation")

    # G6  "some ... some ... some", "exactly that baseline", and a defensive disclaimer that
    #     inverts P8 by putting the contribution in the subordinate clause.
    edit(find_p(body, "A large literature reports that some TE family"), [
        ("that some TE family carries some epigenomic mark near genes of some function",
         "that a given TE family carries a given epigenomic mark near genes of a given "
         "function"),
        ("The results reported here supply exactly that baseline",
         "The results reported here supply that baseline"),
        ("Read this way, the absence of epigenomic data in this study is the point rather "
         "than a gap, since the proximity layer has to be measured on its own before "
         "anything can be conditioned on it.",
         "The proximity layer has to be measured on its own before anything can be "
         "conditioned on it."),
    ], "G6 null-model subsection")

    # G7  "landmark" is praise, not description.
    edit(find_p(body, "The previous landmark repeat-based categorization"), [
        ("The previous landmark repeat-based categorization",
         "The previous repeat-based categorization"),
    ], "G7 landmark")

    # G8  "at all" as filler.
    edit(find_p(body, "Thirty-one of the 44 families"), [
        ("yield no significant functional gene groups at all, which suggests",
         "yield no significant functional gene groups, which suggests"),
    ], "G8 at all")

    # G9  "(343 ones)".
    edit(find_p(body, "The mapping showed that only 0.89% of unique TSS"), [
        ("(343 ones)", "(343)"),
    ], "G9 343 ones")

    # G10a  "beyond the scope of this revision" describes the review process, not the study.
    edit(find_p(body, "Whether the residual differences are methodological"), [
        ("which was beyond the scope of this revision",
         "which lies outside the scope of the present study"),
    ], "G10a scope of this revision")

    # G10b  "several groups" is countable: the raw odds ratio is above 1 for all six classes,
    #       while the observed-to-random ratio is below 1 for LINE, LTR, DNA and Helitrons.
    edit(find_p(body, "Second, once the length bias of the odds ratio is removed"), [
        ("The direction of these effects reverses for several groups if the raw odds ratio",
         "The direction of these effects reverses for four of the six classes if the raw odds "
         "ratio"),
    ], "G10b four of six classes")

    # G10c  "The former / the latter" reach back past their own clause (P4c).
    edit(find_p(body, "The former was associated with a single significant GO term"), [
        ("The former was associated with a single significant GO term",
         "hAT-Tip100 was associated with a single significant GO term"),
        ("(FDR = 0.040), the latter with olfactory receptor activity",
         "(FDR = 0.040), hAT-Charlie with olfactory receptor activity"),
    ], "G10c former/latter")

    # G11  The eight interferon genes, in coordinate-free but sorted order.
    edit(find_p(body, "Those three GO terms share a single core gene set"), [
        ("IFNA10, IFNA16, IFNA17, IFNA21, IFNA4, IFNA6, IFNA7 and IFNW1",
         "IFNA4, IFNA6, IFNA7, IFNA10, IFNA16, IFNA17, IFNA21 and IFNW1"),
    ], "G10 interferon gene order")

    # ============================================================ H. house style (P9)

    for eid, anchor, pairs in [
        ("H1 defense", "in a continuous evolutionary arms race with host",
         [("host defence systems", "host defense systems")]),
        ("H2 artifact", "implausible for its L1P3 clade",
         [("possibly a RepeatMasker artefact", "possibly a RepeatMasker artifact")]),
        ("H3 artifact", "known property of this design rather than a correctable",
         [("a correctable artefact of our implementation",
           "a correctable artifact of our implementation")]),
        ("H4 tumor", "Finally, no eQTL, expression or",
         [("expression or tumour data", "expression or tumor data")]),
        ("H5 tumor", "more than 500 L1 insertions per",
         [("per tumour reported in bladder cancer", "per tumor reported in bladder cancer"),
          ("worth examining first in tumour cohorts", "worth examining first in tumor cohorts")]),
        ("H6 networkx", "Network visualizations were",
         [("network version 3.6.1", "networkx version 3.6.1")]),
        ("H7 ChatGPT", "was used for grammar corrections",
         [("Chat GPT was used", "ChatGPT was used")]),
        ("H8 missing space", "The T2T RepeatMasker annotation itself",
         [("annotation itself was derived from", "annotation itself was derived from ")]),
    ]:
        edit(find_p(body, anchor), pairs, eid)

    # H9  One multiplication sign. The manuscript already writes "2.3 × 10-91" in Results;
    #     the interferon-alpha, prior-work and sensitivity subsections write " x 10".
    times = [
        ("Fisher exact odds ratio 3.01, raw p = 3.2 x 10",
         "The L1 elements of the domain are younger"),
        ("Fisher exact odds ratio 3.01, raw p = 3.2 x 10",
         "Figure 8. The interferon-alpha domain of chromosome 9"),
        ("odds ratio 4.02, FDR = 4.4 x 10", "The like-for-like pairs agree"),
        ("odds ratio 3.69, FDR = 6.1 x 10", "The like-for-like pairs agree"),
        ("odds ratio 2.44, FDR = 6.3 x 10", "The like-for-like pairs agree"),
        ("(4 genes, FDR = 7.6 x 10", "The cross pairs are depleted"),
        ("(12 genes, FDR = 5.8 x 10", "The cross pairs are depleted"),
        ("0.46-fold, FDR = 2.8 x 10", "One counterpair pair diverges"),
        ("the weakest being 3.5 x 10", "The gene sets themselves are stable across windows"),
    ]
    for i, (old, anchor) in enumerate(times, start=1):
        try:
            p = find_p(body, anchor)
        except LookupError:
            # the MIR paragraph opens with a different first clause in some builds
            p = find_p(body, "Our MIR-enriched genes are depleted, not enriched")
        edit(p, [(old, old.replace(" x 10", " × 10"))], f"H9.{i} multiplication sign")

    # ============================================================ C. structural (P7)

    # C1  The duplicate Table 1 / Table 2 pair after Funding. It drops the unadjusted Fisher
    #     p and merges mean and SD into one uncopyable cell; the pair in Results carries both
    #     and is placed at first citation, so the trailing pair is the one to go. Its rows are
    #     all our own tracked insertions, so removing them outright leaves Reject All intact.
    cap1 = find_p(body, "Table 1. Enrichment of TE classes in gene TSS 10 kb neighborhoods.")
    cap2 = find_p(body, "Table 2. Statistical support for TE class enrichment in gene TSS 10 kb")
    kids = list(body)
    t1 = kids[kids.index(cap1) + 1]
    t2 = kids[kids.index(cap2) + 1]
    for el, name in ((t1, "trailing Table 1"), (t2, "trailing Table 2")):
        if el.tag != qn("w:tbl"):
            raise AssertionError(f"C1: expected a table after the caption of {name}")
    remove_inserted_table(t1, "C1 duplicate Table 1")
    remove_inserted_table(t2, "C1 duplicate Table 2")
    drop_own_insertions(cap1, "C1 duplicate Table 1 caption")
    drop_own_insertions(cap2, "C1 duplicate Table 2 caption")

    # C2  The 7x11 table in Results whose every text run is already struck. A row is only
    #     removed by Accept All when its <w:trPr> carries a <w:del>, so as it stands the
    #     accepted document keeps an empty grid between Table 1 and the Table 2 caption.
    empty = [t for t in body.findall(qn("w:tbl"))
             if not [x for x in t.iter(qn("w:t")) if _alive(x)]]
    if len(empty) != 1:
        raise AssertionError(f"C2: expected exactly one fully-struck table, found {len(empty)}")
    mark_rows_deleted(empty[0], "C2 struck-through table")

    # ============================================================ D. Ethical statement (P7)

    # D1  The section paragraph stops mid-clause in the submitted baseline; the sentences it
    #     needs are sitting in the supplementary block instead. Finish the section, then
    #     remove the misplaced copy.
    edit(find_p(body, "This study represents a purely computational analysis"), [
        ("obtained from the T2T genome assembly and the RepeatMasker track",
         "obtained from the T2T genome assembly and the RepeatMasker track. No new "
         "biological material was collected and the work did not involve human participants, "
         "animal subjects or experimental interventions, so no institutional review board "
         "approval or informed consent was required."),
    ], "D1 complete Ethical statement")
    delete_paragraph_smart(find_p(body, "This study analyzed only publicly available genomic data"),
                           "D2 misplaced Ethical statement")

    # D3  The five workbooks are already described one by one in Supplementary material, four
    #     paragraphs earlier. Say it once.
    edit(find_p(body, "The supplementary tables are provided as five thematic workbooks"), [
        ("each opening with a README sheet that describes its contents and names the script "
         "that produced it. File S1: the TE-to-TSS map and the enrichment statistics for "
         "classes, families and subfamilies. File S2: the foreground gene sets used for Gene "
         "Ontology analysis, in long format. File S3: the Gene Ontology results at FDR < "
         "0.05 with their functional-group classification. File S4: the interferon-alpha "
         "domain, the newly resolved sequence and the overlap with prior work. File S5: the "
         "window-size, percentile and Gene Ontology grid sensitivity analyses. The "
         "enrichment sheets",
         "described one by one in the Supplementary material section, each opening with a "
         "README sheet that names the script that produced it. The enrichment sheets"),
    ], "D3 duplicated workbook list")


# ----------------------------------------------------------------- verification

def verify(src: Path, dst: Path) -> None:
    import zipfile

    def counts(path):
        with zipfile.ZipFile(path) as z:
            xml = z.read("word/document.xml").decode("utf-8")
        return {
            "sdt": xml.count("<w:sdt>"),
            "MENDELEY_CITATION": xml.count("MENDELEY_CITATION"),
            "MENDELEY_BIBLIOGRAPHY": xml.count("MENDELEY_BIBLIOGRAPHY"),
            "hyperlink": xml.count("<w:hyperlink"),
            "ins": xml.count("<w:ins "),
            "del": xml.count("<w:del "),
            "t_in_del": xml.count("<w:del ") and None,
        }

    a, b = counts(src), counts(dst)
    print("\n--- integrity ---")
    for k in ("sdt", "MENDELEY_CITATION", "MENDELEY_BIBLIOGRAPHY", "hyperlink"):
        flag = "OK" if a[k] == b[k] else "*** CHANGED ***"
        print(f"  {k:24s} {a[k]:>5} -> {b[k]:>5}   {flag}")
        if a[k] != b[k]:
            raise AssertionError(f"{k} count changed: {a[k]} -> {b[k]}")
    print(f"  {'<w:ins>':24s} {a['ins']:>5} -> {b['ins']:>5}")
    print(f"  {'<w:del>':24s} {a['del']:>5} -> {b['del']:>5}")

    doc = Document(dst)
    bad = []
    for d in doc.element.body.iter(qn("w:del")):
        if d.find(".//" + qn("w:t")) is not None:
            bad.append(d)
    print(f"  <w:t> inside <w:del>      {len(bad)}   {'OK' if not bad else '*** BAD ***'}")
    if bad:
        raise AssertionError("a deletion carries live <w:t>; it must be <w:delText>")

    tables = doc.element.body.findall(qn("w:tbl"))
    live = [t for t in tables if [x for x in t.iter(qn("w:t")) if _alive(x)]]
    print(f"  tables total / with live text   {len(tables)} / {len(live)}   "
          f"{'OK' if len(live) == 3 else '*** expected 3 ***'}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", action="store_true", help="print the last run's report only")
    args = ap.parse_args()

    if args.report:
        if not DST.exists():
            sys.exit("no output yet; run without --report first")
        print(f"{DST.name} exists ({DST.stat().st_size:,} bytes)")
        return

    if not SRC.exists():
        sys.exit(f"missing input: {SRC}")
    set_revision_identity(AUTHOR, DATE)
    # word_rewrite_trackchanges numbers revisions from a fixed base of 1000, and three earlier
    # passes have already used that range, so start above the highest id already in the file.
    # Duplicate w:id values are what made Word's revision pane unreliable in the first attempt.
    import re
    import zipfile
    existing = [int(x) for x in re.findall(
        r'<w:(?:ins|del)\s[^>]*w:id="(\d+)"',
        zipfile.ZipFile(SRC).read("word/document.xml").decode("utf-8"))]
    wrt._rid[0] = max(existing or [1000]) + 1
    print(f"revision ids start at {wrt._rid[0] + 1} "
          f"({len(existing)} revisions already in the input)")

    shutil.copy2(SRC, DST)
    doc = Document(DST)
    apply_edits(doc)
    doc.save(DST)

    width = max(len(e) for e, _, _ in report)
    print(f"--- {len(report)} operations ---")
    for eid, status, detail in report:
        print(f"  {eid:{width}}  {status:10}  {detail}")
    bad = [r for r in report if r[1] not in ("ok", "ok-in-ins")]
    print(f"\n{len(report) - len(bad)} of {len(report)} succeeded")
    if bad:
        sys.exit("some edits did not apply")
    verify(SRC, DST)
    print(f"\nwrote {DST}")


if __name__ == "__main__":
    main()
