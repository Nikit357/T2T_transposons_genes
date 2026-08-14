# Figure statistics labels — raw vs FDR-corrected, per panel

**Manuscript:** G3-2026-406828 · **Deliverable of:** WP9 (reviewer minor comment 3)
**Written:** 2026-08-04 · **Companion to:** `revision_G3/svg/PLACEMENT.md`

Reviewer minor comment 3: *"…clearly indicate in axis labels or legends whether p-values are
raw or FDR-adjusted… In Figure 1D, for example, the third vertical bar plot labeled
'-log10(p-value)' should read '-log10(FDR-corrected p-value)'."*

This is also **journal policy**. G3's statistics guidance (`G3_article_guidelines.md` §5)
requires that it be clear whether reported p-values are raw or corrected, that the type of
correction be named, **and that raw p-values be available in the supporting materials** so a
reader can apply their own correction. Our tables already carry both columns, so the obligation
is met by saying so — see §4.

The **Correct label** column is the verbatim string to put on the artwork. Every panel listed
under "regenerated" already carries it; the panels under §3 need editing by hand in Figma.

---

## 1. Regenerated panels — label already correct in the SVG

| Figure / panel | Frame | Statistic | Source column | Correct label | Raw or corrected | Appears as |
|---|---|---|---|---|---|---|
| 1D bar 1 (OR) | `856:7` | empirical permutation p, BH-adjusted | `p_adjusted_empirical_bh` | `Observed OR vs random control (ns: FDR-corrected empirical p ≥ 0.05)` | **corrected** | panel title + legend entry |
| 1D bar 2 (fold change) | `856:7` | — (descriptive ratio) | `OR_Observed_to_Random` | `Fold change (observed / random)` | n/a | panel title |
| 1D bar 3 (significance) | `856:7` | Fisher exact, BH-adjusted | `Enrichment_p_value_adjusted` | title `Significance -log10(FDR-corrected p)` · x-axis `-log10(FDR-corrected Fisher exact p-value)` | **corrected** | title + axis label |
| 1D bars 4–5 | `856:7` | — (descriptive) | `Total_TEs_number_log`, `Average_divergence_all` | `log10(total elements)`, `Divergence score` | n/a | axis labels |
<!-- NAMING HAZARD: `8C` and `S8C` are two different panels one letter apart. `8C` is the
     interferon-alpha subfamily-composition panel of the new main Figure 8. `S8C` is the
     unfiltered class -> functional group -> family Sankey, panel C of Supplementary Figure
     S8, added 2026-08-05. Every reference must carry the `S`, including in the response
     letter. -->

| 2D / 2E / 2F boxes | `856:8` | per-**family** empirical p, BH-adjusted | `p_adjusted_empirical_bh` | `family significant (FDR-corrected empirical p < 0.05)` | **corrected** | in-panel legend (circle vs cross) |
| 2D / 2E / 2F matrix | `856:8` | Mann–Whitney U, BH across the **testable** class pairs | computed in `nb03` | colourbar `-log10(FDR-corrected p)` · caption strip `Mann–Whitney U, Benjamini–Hochberg corrected across the <n> testable class pairs · *** FDR < 0.001 ** < 0.01 * < 0.05 ns ≥ 0.05 · n/a: a class with a single family` | **corrected** | colourbar + x-axis strip |

**The Figure 2 panels changed level and letter order on 2026-08-05.** They plotted the
subfamily table; the published caption is family-level, and it assigns the letters
D = observed/random OR, E = observed OR, F = random OR — the reverse of what was drawn for D
and E. Both are corrected in `nb03`, whose cell 17 now asserts the three published
significance statements rather than trusting them. Two label consequences: the point-marker
legend says **family**, not subfamily, and because Retroposon and RC have a single family each
their 9 class pairs are untestable and are annotated **`n/a`** rather than `ns` — so the pair
count in the caption strip is no longer the flat 15.
| 3A / 3B | `856:9` | two-sample Kolmogorov–Smirnov | computed in `nb03` | `Two-sample Kolmogorov–Smirnov, p = <value> (raw, single test)` | **raw** | panel title |
| 4A / 5A / 6A node colour | `859:25`, `861:28`, `861:33` | goatools GO enrichment, `p_fdr_bh` | `FDR` | `-log10(FDR-corrected GO enrichment p-value)` | **corrected** | colourbar (vector version: `Fig456A_colorbar_vector.svg`) |
| 4B / 5B / 6B / S8B stars | `859:25`, `861:28`, `861:33`, `861:32` | Fisher exact per cell, BH across all cells | computed in `nb06` | `Functional group · stars: Fisher exact, FDR-corrected across all cells · *** FDR < 0.001 ** < 0.01 * < 0.05 · GO terms counted at FDR < 0.05` | **corrected** | x-axis strip |
| 4B / 5B / 6B / S8B colourbar | as above | — (count scale) | GO term counts at FDR < 0.05 | `log10(GO terms at FDR < 0.05 + 1)` | n/a | colourbar |
| 7A, 7C, 7D | `861:34` | Pearson correlation | `output/figure7_statistics.csv` | in-panel box `Pearson R = <value>, p = <value> (raw)` · y-axis `GO terms in a group (FDR < 0.05)` | **raw** correlation of **corrected** counts | in-panel annotation + axis |
| 7B, 7E, 7F, 7G | `861:34` | Mann–Whitney U | `output/figure7_statistics.csv` | `Mann–Whitney U, raw p = <value> (single test); n = <a> vs <b>` | **raw** | **figure caption** (manuscript edit M9) — *not* a panel title |
| 7H | `861:34` | Fisher exact per family × functional group vs its own class, BH-corrected | computed in `nb06` | footer strip `Stars: Fisher exact of a family x functional group against its own class, BH-corrected` + `Ribbons with fewer than 5 GO terms are not drawn (visualisation only; bar heights are unfiltered)` | **corrected** | footer strip |
| S8C | `861:32` | same as 7H | computed in `nb06` | same footer strip **without** the ribbon-filter sentence | **corrected** | footer strip |
| 7H | `861:34` | Fisher exact per family × functional group against its own class, BH-corrected | computed in `nb06` | footer strip `Ribbon heights are log10(GO terms at FDR < 0.05 + 1); ribbon labels shown for ≥ 5 GO terms. Stars: Fisher exact of a family × functional group against its own class, BH-corrected *** FDR < 0.001 ** < 0.01 * < 0.05` | **corrected** | figure footer |
| 8B | new Figure 8 | empirical permutation p, 10,000 matched windows | `empirical_p_two_sided` | per panel: `observed <v> vs null <m> ± <sd>; z = <z>, empirical p = <p> (raw)` and `<n> random 220 kb windows` | **raw** | in-panel annotation |
| **8C** (not S8C) | new Figure 8 | Fisher exact 2×2, young vs old L1 | `empirical_p_two_sided` of `T4` | `Young vs old L1 composition: OR = <or>, Fisher p = <p> (raw)` | **raw** | panel title |
| S1 | `856:10` | Kolmogorov–Smirnov per class | computed in `nb03` | per panel `K-S p = <value>` + footer `Two-sample Kolmogorov–Smirnov per class; p-values are raw (6 independent tests, no multiplicity correction applied)` | **raw** | panel titles + figure footer |
| S2 | `856:11` | Mann–Whitney U, BH across 44 families | `revision_G3/output/S2_family_mannwhitney_fdr.csv` | legend `Mann–Whitney U, BH-corrected across the 44 families: *** FDR < 0.001, ** < 0.01, * < 0.05, ns ≥ 0.05` | **corrected** | figure legend |
| S3 | `856:12` | Mann–Whitney U, BH across 44 families | `revision_G3/output/S3_family_mannwhitney_fdr.csv` | same string as S2 | **corrected** | figure legend |
| S9 / S10 / S11 | new frames | goatools `p_fdr_bh` | `FDR` | same colourbar string as 4A | **corrected** | colourbar |
| S13A | new frame | — (descriptive OR) | `or_observed_to_random` | `Observed / random odds ratio` | n/a | axis label |
| S13B | new frame | Spearman ρ per pair | `window_sensitivity_concordance.csv` | `<pair> bias <b>, ρ = <rho>` | **raw** | panel titles |
| S13C | new frame | — (set arithmetic) | `percentile_sensitivity_summary.csv` | `GO terms at FDR < 0.05`, `Fraction of the 5 % terms still significant at 10 %` | n/a | axis labels |
| S13D | new frame | goatools `p_fdr_bh` | `FDR` of `go_grid_headline_by_condition.csv` | horizontal colourbar `-log10(FDR-corrected GO enrichment p-value)` + legend `not significant at FDR < 0.05, or absent` | **corrected** | colourbar + legend |
| S12A | new frame | — (term counts) | `GO_grid/INDEX.csv` | colourbar `GO terms relative to 10 kb / 5 %` | n/a | colourbar |
| S12B | new frame | — (set arithmetic) | `go_grid_preservation.csv` | legend `published terms preserved`, `Jaccard of the term sets` | n/a | legend |
| S12C | new frame | goatools `p_fdr_bh` | `GO_grid/*_fdr005.csv` | legend `at least one GO term at FDR < 0.05` | **corrected** | legend |
| S13E | new frame | hypergeometric / Kendall τ | `robustness_geneset_stability.csv`, `robustness_rank_stability.csv` | `Overlap coefficient of the top-5 % gene sets`, `Kendall τ of the gene ranking (95 % bootstrap CI)` | **raw** (each a single test) | axis labels |
| S14A–S14C | new frames | — (convergence diagnostics) | `permutation_convergence_checkpoints.csv` | `Running mean random OR (fraction of the N = 500 value)` etc. | n/a | axis labels |

---

## 2. Which statistic is which — one line each

| Statistic | Where it comes from | Correction |
|---|---|---|
| Fisher exact enrichment p | class/family/subfamily TE enrichment in TSS windows | BH across the groups tested at that level |
| Empirical permutation p | observed OR against the N = 500 random background | BH across the groups; floor 2/501 = 0.004 |
| goatools GO enrichment p | `GOEnrichmentStudy`, `p_fdr_bh` | BH inside each study, threshold 0.05 (D2) |
| Mann–Whitney U | class-pair and family comparisons | BH across the 15 pairs (Fig 2) or 44 families (S2/S3) |
| Kolmogorov–Smirnov | all TEs vs TEs near genes, divergence | **none** — reported raw |
| Pearson r | GO term count vs enrichment statistics (Fig 7) | **none** — reported raw |
| IFNA empirical p | 10,000 matched random 220 kb windows | **none** — one test per panel, raw |
| Hypergeometric / Kendall τ | gene-set and rank stability (S13E) | **none** — reported raw |

---

## 3. Still to edit by hand in Figma

1. **Frame `861:34` (Figure 7)** contains the literal string
   `"GO terms count in a group (FDR < 0.1)"` → must read **`FDR < 0.05`**. The regenerated
   panels already say 0.05, so leaving the frame text would make the figure contradict its own
   axes and the revised Methods (caveat C11, gap G17).
1b. **Figure 7's panels C/E/G moved and E/F/H were replaced** when the panel identities were
   corrected against the published caption on 2026-08-04 (see `PLACEMENT.md` §1). Swap all
   eight together rather than panel by panel, or the letters will not match the caption.
2. **Grep every frame** for `0.1`, `500` and `p-value` once the panels are swapped. Analysis
   parameters baked into figure text do not update when the CSVs do, and this is the most
   likely route to a figure that contradicts the Methods.
3. **Captions.** Each caption should name the test and the correction once, matching the
   artwork. Panels whose p-values are raw (3A/3B, 7A–7G, 8B, 8C, **8D**, S1, **S12B, S13A**)
   must say so in the caption as well as on the panel.
   *Renumbered 2026-08-07:* the last two were S13B and S13E before old Supplementary Figure 6
   became Figure 8A and the supplementary set became S1–S13, and 8D is the leave-one-out panel
   promoted from an inset. See `revision_G3/svg/PLACEMENT.md` §0. Item 1 above — the baked-in
   `FDR < 0.1` in frame `861:34` — is **done**, verified against the exported PDFs.

---

## 4. The Methods and Data availability sentences this requires

Add to **Methods**:

> Unless stated otherwise, all reported p-values are Benjamini–Hochberg FDR-adjusted; raw
> p-values are labelled as such in the figures and captions. Kolmogorov–Smirnov comparisons of
> divergence distributions, the Pearson correlations of GO term counts against enrichment
> statistics, and the interferon-alpha domain permutation tests are single tests per panel and
> are reported raw.

Add to **Data availability** (G3 requires raw p-values to be obtainable):

> Raw and FDR-adjusted p-values are both provided for every test:
> `enrichment_families_with_random.csv` and `enrichment_subfamilies_with_random.csv` carry
> `Enrichment_p_value` / `Enrichment_p_value_adjusted` and `p_raw_empirical` /
> `p_adjusted_empirical_bh`; the GO tables carry `P-value` and `FDR`;
> `TableS_class_enrichment_full.csv` carries the unadjusted Fisher p-values omitted from
> Table 2.
