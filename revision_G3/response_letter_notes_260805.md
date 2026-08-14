# Response-letter notes — things the letter must say explicitly

Written 2026-08-05, from the review-corrections pass. These are the points where a reviewer
comparing the April submission against the revision would otherwise find a silent difference.
Every number here is re-derivable: `revision_G3/output/results_numbers.json`,
`go_grid_summary.json`, `network_qc.json`, `sankey_ribbon_filter.json`, `supplementary/INVENTORY.json`.

---

## 1. The supplementary material changed in two ways at once (caveat S6)

Both are improvements, but a reviewer working from their own notes will find the numbering no
longer matches, so say it plainly rather than leaving it to be discovered.

**a. Every Gene Ontology result moved from FDR 0.1 to FDR 0.05**, in the supplementary tables as
well as the main text, with no "suggestive" band anywhere (this is minor comment 2). The April
supplementary files were at 0.1. Measured effect: classes-by-count 504 → 425 terms,
classes-by-divergence 516 → 414, families 195 → 140, subfamilies 1,231 → 1,003; families with at
least one term 14 → **13 of 44** (Dong-R4 is the one lost).

**b. Fourteen candidate files became five thematic workbooks** (Files S1–S5), each with one sheet
per table and a README sheet describing them. So every `File Sn` citation in the text is
renumbered and now names a sheet as well as a file. The mapping is in
`revision_G3/supplementary/README.md` and should be reproduced in the letter:

| April 2026 | Now |
|---|---|
| File S1 (TSS/TE coordinates) | File S1, sheet `TSS_TE_intersections` |
| File S2 (enrichment) | File S1, sheets `enrichment_families` / `enrichment_subfamilies` |
| File S3 (gene sets by TE group) | File S2, sheet `by_TE_group` |
| File S4 (GO by TE group) | File S3, sheet `GO_TE_groups` |
| File S5 (gene sets by divergence) | File S2, sheet `by_divergence` |
| File S6 (GO by divergence) | File S3, sheet `GO_by_divergence` |
| File S7 (gene sets by family) | File S2, sheet `by_family` |
| File S8 (GO by family) | File S3, sheet `GO_by_family` |
| — (Lu et al. overlap) | File S4, sheet `prior_work_overlap_matrix` |
| — ("the accompanying tables") | File S5 |

Three further changes inside that package, each worth a sentence:

* **The old File S2 caption/content mismatch is resolved by carrying both tables.** The text
  described File S2 as "enrichment statistics of TE subfamilies"; the file contained the 44
  families. Workbook S1 now ships `enrichment_families` **and** `enrichment_subfamilies`, both
  with raw and adjusted p-values, so the caption and the Data availability claim both become true
  without anything being dropped.
* **Two sheets of the April File S3 were empty** — `TE top` and `TE bottom`. Those gene sets drive
  GO results the paper reports, so only their inclusion was missing, not the sets. They are
  reconstructed by the same construction the GO scripts use, verified first against the six
  non-empty sheets: every disagreement sits exactly on the tie boundary (for the LINEs set, 1,257
  genes are strictly above the boundary count of 9 and 1,033 tie at it for 179 remaining places),
  and no gene above the boundary differs in any set.
* **The gene-set sheets are long format**, which also restores the `hAT?` and `hAT-Tip100?` family
  names — Excel forbids `?` in a sheet name, so the April File S7 had shipped them as `hAT...` and
  `hAT-Tip100...`.

---

## 2. The window-size sensitivity analysis is now complete, and the result is not tidy

Major comment 5 asked for alternative window sizes **and** alternative percentiles. The April
revision answered the window half for enrichment and the percentile half for GO; GO had never been
run at 5 kb or 20 kb. It now is — the full 3 windows × 2 percentiles grid, all three GO levels, 18
cells. Report the outcome as measured:

* Widening the **percentile** always finds more terms (9 of 9 level × window combinations).
  Widening the **window** does not, and not even in the same direction at every level:
  classes-by-count falls 510 → 425 → 299, classes-by-divergence peaks at 10 kb (209 / 414 / 303),
  families rises 137 → 140 → 201. A wider window adds elements to *every* gene, so it does not
  simply add power — it changes which genes are in the top 5 % and dilutes the promoter-proximal
  signal with more distal sequence.
* Of the terms the paper reports, the fraction still significant in another cell is at worst
  **0.440** and at median **0.677** — weaker than the percentile-only result (0.85–0.93). **The
  TSS window matters more to GO than the gene-set percentile does.** Say so.
* **3 of the 9 headline claims survive all six conditions**; all 9 hold in the published 10 kb /
  5 % condition. The weakest are `SVA / termination of RNA polymerase II transcription` (2 of 6)
  and `hAT-Charlie / MHC class I protein complex` (1 of 6). The hAT-Charlie term was already
  flagged as percentile-fragile; the grid makes the case for softening or dropping it stronger.
* Concordance of per-group term counts against the published cell: lowest Spearman ρ **0.614**,
  every label-shuffling permutation p ≤ **0.022**.
* **The grid is GO only — no permutations were re-run.** So a difference between cells is a
  gene-set effect, not a background effect, and the enrichment odds ratios of Table 1 and Figure 2
  are unchanged. State this or the two analyses will be conflated.

---

## 3. Three figure corrections that make the artwork agree with captions the journal already has

None of these creates a new claim; each removes a disagreement.

* **Figures 2D–2F** were plotted at subfamily level with D and E swapped. The published caption is
  family-level and assigns D = observed/random OR, E = observed OR, F = random OR. At family level
  the panels reproduce the three published significance sentences exactly, so **no text changes**.
  Two visible consequences: the `n=` per class now reads 1–22 instead of 2–600, and Retroposon and
  RC have one family each, so they appear as single points and their class pairs are marked `n/a`
  rather than `ns`.
* **Figure 7H** now applies the ribbon filter its own caption already promised ("Connecting ribbons
  were filtered by at least 5 GO terms. This filtering was applied to the visualization only").
  The previous version applied the threshold only to the ribbon count *labels*. Filtering hides 36
  class → group and 50 group → family ribbons, 146 GO terms. Because the filter is visual only,
  bar heights are unfiltered, so retained ribbons do not fill their bars — which the caption
  already accounts for. The unfiltered version is now **Supplementary Figure S8C**.
* **Supplementary Figures S9–S14 receive captions for the first time.** Only S1–S8 were captioned;
  S9–S11, S13 and S14 were cited without captions, and S12 was neither cited nor captioned. The
  new GO-grid figure takes S12, so the final inventory is a contiguous S1–S14.

---

## 4. Legibility: what was done, and one thing that could not be

Minor comment 6 asked for legible figures. All panel text was raised by exactly 1.2×, which puts
the base size at 16 px as rendered — matplotlib emits points and SVG consumers convert at 96/72,
so the previous 10 pt arrived as 13.33 px.

The trade-off is worth disclosing if the point is raised: **denser packing of the network panels
was not achievable alongside the larger text.** A fixed number of labels at a fixed size needs a
minimum area; 1.2× text inflates every label box by 1.44× in area, so packing 30 % tighter would
mean the same labels in 0.49× the area they already required. The no-overlap requirement is
enforced programmatically — a figure with overlapping labels cannot be written — so the resolution
was to keep the larger text and let the canvas grow, which costs nothing in a vector figure.

One consequence to state in the Figure 6A legend rather than leave implicit: at the larger text
size that panel shows **five** GO terms per family, not nine.
