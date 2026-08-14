# G3: Genes|Genomes|Genetics — author guidelines, distilled for manuscript G3-2026-406828

**Purpose.** A working reference for formatting the accepted manuscript
("Evolutionary arms race between transposable elements and human genes…", G3-2026-406828) to G3 house
style, with an explicit gap list against our current draft. Companion to
`G3_revision_implementation_plan_260803.md` (see WP14 there).

## 0. Provenance and a warning about it

The live page `https://academic.oup.com/g3journal/pages/author-guidelines?login=false#section-11`
**cannot be retrieved from this environment**: OUP returns HTTP 403 to every request (plain curl,
browser user-agent, Googlebot user-agent, reader proxies), and even the Internet Archive's captures of
that URL are OUP's "crawl prevention" placeholder rather than the page content. This document is
therefore assembled from three sources that *are* verifiable, each statement tagged:

| Tag | Source | Reliability |
|---|---|---|
| `[GSA-A]` | GSA's own "Preparing Manuscripts for Submission" page, retrieved from the Internet Archive (`web.archive.org/web/2021/https://www.g3journal.org/content/prep-manuscript`) | Substantive policy is still GSA's; some *print-style* details predate the 2021 move to OUP and are superseded — flagged where they conflict with `[G3-2026]` |
| `[G3-2026]` | Structure, headings, citation style and reference formatting extracted programmatically from five G3 research articles published in 2026 (PMC13365844, 13334165/167/169/170/186) | Highest reliability for *current* house style — this is what G3 actually prints today |
| `[SEARCH]` | Verbatim policy statements from the indexed current author-guidelines page returned via web search | Reliable but fragmentary |
| `[LETTER]` | The editor's decision letter of July 29, 2026 (`G3_reviewer_report_260802.md`) | Authoritative for this manuscript |
| `[VERIFY]` | Item that **Daniil must confirm on the live page** — either the archived source may be stale, or no accessible source covers it | Do not treat as settled |

> `[LETTER]` explicitly warns: *"These guidelines have been updated with new requirements and your
> careful attention is required to avoid delays."* Treat §9 (the verify list) as mandatory homework, not
> optional.

---

## 1. Article type and length

- Our paper is an **Investigation** (full-length research article). `[GSA-A]`
- **Investigations have no minimum or maximum length.** `[GSA-A]` So Reviewer 1's request to shorten the
  Discussion is an editorial/quality request, not a length-limit compliance issue — but it still must be
  honoured.
- Observed in current G3 Investigations: 49–65 references, 4–10 figures, 0–2 tables. `[G3-2026]`
  Our draft has **123 references, 9 figures (10 after the revision), 1 table (2 after the split)** —
  the reference count is roughly double the norm, which independently supports the Discussion trim
  planned in WP8.

---

## 2. Section order

**Body** (this is the current printed order): `[G3-2026]`

1. Title
2. Author(s) and affiliation(s)
3. Abstract
4. Keywords
5. Introduction
6. **Materials and methods** — note the **sentence case**; G3 prints "Materials and methods", not
   "MATERIALS AND METHODS" `[G3-2026]`
7. Results
8. Discussion (or "Results and discussion" combined; "Conclusions" is optional and used by some papers)
9. Supplementary material (a short body section pointing at the figshare deposit)

**Back matter, in this order:** `[G3-2026]`

10. Acknowledgments
11. **Data availability**
12. Funding
13. Conflicts of interest
14. (optional) Code availability — one of the sampled 2026 papers used a separate section
15. **Literature cited** — this is the heading G3 uses, *not* "References"

`[GSA-A]` places the data/reagent availability statement *at the end of Materials and Methods*; every
2026 article I sampled prints it as a **separate back-matter section** after Acknowledgments.
**Do both**: keep a one-line pointer at the end of Methods and the full statement as a back-matter
section. `[VERIFY]` which the current guidelines mandate.

`[GSA-A]` also specifies four legacy heading levels (Level 1 centred all-caps, Level 2 flush-left bold,
Level 3 paragraph-initiating bold + colon, Level 4 paragraph-initiating italic). Current G3 typesetting
does not use centred all-caps. Submit with plain Word heading styles (Heading 1/2/3) and let the
copyeditors apply house style. `[VERIFY]`

---

## 3. Manuscript file mechanics

- English, **American spelling**, correct grammar and punctuation. `[GSA-A]`
- **12-point type, double-spaced throughout** — including Literature cited and any appendices. `[GSA-A]`
- **Consecutive page numbers** beginning with the cover page; **line numbers** for review. `[GSA-A]`, `[LETTER]`
- Title: concise, informative, avoids jargon. `[GSA-A]` (For non-human-organism papers the organism name
  must be in the title; ours is human, so N/A — but the abstract must still name the organism, see §4.)
- Affiliations: department; institution; city; state/province; country; postal code. **No street
  addresses.** `[GSA-A]`
- Multiple affiliations marked with `*, †, ‡, §` (then doubled `**, ††, ‡‡, §§`). Single-author paper →
  not applicable. `[GSA-A]`
- **Short running title of ~35 characters including spaces** — we do not have one; must be added.
  Suggested: `TEs near human genes in T2T` (28 characters). `[GSA-A]`
- Corresponding author block: name, office mailing address with street, phone, email. `[GSA-A]`
  Our draft gives only an email (`danya.nikitin.orel@gmail.com`) — needs the full block, and consider
  using the institutional address at the Institute of Molecular Biology, NAS Armenia.
- ORCID for the corresponding author: standard OUP requirement; yours is 0000-0003-1029-1174. `[VERIFY]`
- Editorial style baseline: **Scientific Style and Format (Council of Science Editors, CSE)**. `[GSA-A]`

---

## 4. Abstract and keywords

- **Single paragraph, < 250 words.** `[GSA-A]` Observed range in 2026 papers: 174–264 words. `[G3-2026]`
  Ours is **205 words** — compliant.
- Must stand alone; begin with broad context, then specific background, then purpose, methods, core
  findings, conclusions; emphasise what is new; avoid jargon; **contain the full name of the organism
  studied** (i.e. "human", *Homo sapiens*). `[GSA-A]`
- Keywords: 3–10 observed. `[G3-2026]` Our draft has none — must be added. Suggested:
  *transposable elements; T2T-CHM13; transcription start site; gene ontology; interferon alpha; LINE-1;
  human genome*.

---

## 5. Statistics reporting — this is also Reviewer 1's minor comment 3

G3's own policy, verbatim in substance: `[GSA-A]`

- State **the method and model applied**, not merely the software and options.
- When many genes/phenotypes are tested, **apply a multiple-comparison correction** or justify not
  doing so; **state which correction** was used.
- **"It should also be clear whether the p-values reported are raw, or after correction."**
- **"Corrected p-values are often appropriate, but raw p-values should be available in the supporting
  materials so that others may perform their own corrections."**
- Large-scale exploration studies must describe the replication structure completely.

**Consequences for us.** Reviewer minor comment 3 is not a stylistic preference — it is journal policy.
Two obligations follow: (i) every figure axis/legend and every table column must say whether the value
is raw or FDR-adjusted (WP9 in the plan); (ii) **the supplementary files must expose the raw p-values
alongside the adjusted ones.** Our enrichment CSVs already carry both
(`Enrichment_p_value` and `Enrichment_p_value_adjusted`, `p_raw_empirical` and
`p_adjusted_empirical_bh`) and the GO tables carry both `P-value` and `FDR` — so we satisfy (ii)
already, but the Data availability statement should say so explicitly. Add one sentence to Methods:
*"Unless stated otherwise, reported p-values are Benjamini–Hochberg FDR-adjusted; raw p-values for
every test are provided in Files S2, S4, S6 and S8."*

---

## 6. Figures

- **One file per figure**, all panels of a multipanel figure on a single page in that one file. `[GSA-A]`
- **Panels labelled with a letter (A, B, C…) in the upper-left corner of each panel.** `[GSA-A]`
- **Vector formats for charts, graphs, diagrams and maps: `.eps`, `.ai`, `.pdf`** (resolution does not
  apply to vector art). `[GSA-A]` Our pipeline already emits SVG + PDF with
  `plt.rcParams["svg.fonttype"] = "none"` — export the **PDF** version for submission.
- Raster: **≥ 350 dpi** for images/photographs, **≥ 600 dpi** for line art (1200 dpi for fine line art),
  72 dpi only for on-screen. `[GSA-A]`
  **Where this actually bites for us** (verified 2026-08-03 by inspecting the Figma source, see the
  implementation plan §5b): the UCSC Genome Browser panels are **vector**, not screenshots — Figure 1C
  and Supplementary Figure 6 contain no raster fills at all. The raster exposure is (i) small pasted
  **colorbars** (23 × 58 to 30 × 175 px, `scaleMode=CROP`) inside Figures 4, 5, 6 and Supplementary
  Figure 8 — better re-emitted as vector from matplotlib than rescaled; and (ii) the **schematics**,
  Figure 8 (3 bitmaps) and Figure 9 (10 bitmaps, two carried over from the PhD thesis).
- **`.doc`/`.docx` and `.jpeg`/`.jpg` figures are not accepted.** `[GSA-A]`
- Figure **titles and legends go in the manuscript file, not in the image file.** `[GSA-A]`
  Our figures are currently exported as standalone PDF/PNG with no embedded captions — correct
  already; confirm no title text is baked into the images (our matplotlib `fig.suptitle(...)` calls
  *do* bake titles in, e.g. "Detailed TE Family Enrichment Analysis" on Figure 1D — consider removing
  the suptitles).
- Labels and legends in a **sans-serif 10-point font**. `[GSA-A]` Our `GLOBAL_FONT_SIZE = 10` matches;
  confirm the family is sans-serif (matplotlib default DejaVu Sans is).
- Number figures with Arabic numerals in citation order; cite every figure in the text. `[GSA-A]`
- Distinguish similar glyphs in labels (l vs 1, O vs 0). `[GSA-A]`
- Legends start with a **brief title**, then a self-contained description sufficient to understand the
  data. Italicise mathematical variables, genotypes and other normally-italic symbols in both the
  legend and the figure. `[GSA-A]`
- No image manipulation that misrepresents data; editors may request originals. `[GSA-A]`
- Figures are auto-converted into downloadable PowerPoint slides — another reason to keep the
  simplified network panels legible (Reviewer minor 6 / WP12). `[GSA-A]`

---

## 7. Tables

- **All tables go at the end of the main text**, in an **editable format** (Word or LaTeX) — never as
  an image. `[GSA-A]`
- Data only. No shading, colour type, line drawings or graphics inside tables. `[GSA-A]`
- Title required, concise, with the table number in Arabic numerals. **Do not number tables 1A, 1B** —
  interior parts may be labelled A, B for reference. `[GSA-A]`
  → This directly supports the WP10 decision to split the over-wide Table 1 into **Table 1 and
  Table 2**, not into "Table 1A/1B".
- Footnotes directly below the table, lowercase superscript italic letters (a, b, c…). Use
  `*`, `**`, `***` for conventional significance levels, explained below the table. `[GSA-A]`
- Top and bottom rules at 0.5 pt. Labels in all caps, sans-serif 10 pt; table body sans-serif 10 pt,
  double-spaced; totals in bold. `[GSA-A]` `[VERIFY]` (this is legacy print styling; copyeditors may
  handle it)
- Legend, if present, precedes any footnotes. `[GSA-A]`

---

## 8. References — **our largest formatting gap**

### 8.1 In-text citations: author–year, not numbers

- Two authors → both names. Three or more → first author + `et al.` `[GSA-A]`
- Example format: `(Chen et al. 1997; Scott and Rogers 1998; Isaacson 1999)` `[GSA-A]`
- Confirmed in current print: `(Meuwissen et al. 2001)`, `(Wang et al. 2021)`,
  `(Pieruschka and Schurr 2019)` `[G3-2026]`
- Same first author, same year → `1996a`, `1996b`. `[GSA-A]`
- Cite only preprints, published, or in-press works. Personal communications / unpublished results:
  list all contributors by initials and last name, never `et al.` `[GSA-A]`

**Our draft uses numbered citations throughout** — `(1)`, `(14)`, `(20–22)`, and reference-number
mentions embedded in sentences such as *"the previous landmark study by (14)"* and
*"A landmark study by (10) that showed…"*. This must be converted to author–year, and those
grammatically broken constructions rewritten (*"the landmark study by Lu et al. (2020)"*).
G3 provides **EndNote and Zotero style files for CSE style**. `[GSA-A]`

**How this is actually done in our manuscript (verified 2026-08-03).** The reference manager is
**Mendeley Cite**, the Word web add-in — not a numbered-text convention and not the legacy COM plugin.
The manuscript carries **128 in-text `<w:sdt>` citation content controls** plus one tagged
`MENDELEY_BIBLIOGRAPHY`; the payload lives in a single `MENDELEY_CITATIONS` property inside
`word/webextensions/webextension1.xml`, and the active style is
`MENDELEY_CITATIONS_STYLE` = *NLM/Vancouver: Citing Medicine 2nd edition (citation-sequence)*,
`format: numeric`, locale `en-GB`. There are **zero** legacy field codes (`ADDIN CSL_CITATION`,
`instrText`, `fldSimple`) anywhere in the archive. Consequences:

- **The style conversion is a Mendeley Cite operation, not a text edit** — switch the style to the
  G3/CSE author–year style in the Mendeley Cite pane and refresh (the bibliography is already flagged
  `MENDELEY_BIBLIOGRAPHY_IS_DIRTY = true`). **Daniil owns this step (G1).**
- **The in-sentence grammar sweep is a separate, manual job (G1b)**, because that prose sits *outside*
  the content controls and Mendeley will not fix it.
- **Any text edit that moves or deletes a citation must move or delete the whole `<w:sdt>` element.**
  Editing only the visible run leaves an orphaned control or dead plain text that will not re-render.
  Full mechanics in `G3_revision_implementation_plan_260803.md` §3.5 and caveat C10.
- Switch the locale from `en-GB` to US English at the same time — it is the source of the British
  spellings G3 style rejects (G7).

### 8.2 Literature cited entries

Current printed format `[G3-2026]`:

```
Allentoft ME, et al. 2015. Population genomics of bronze age Eurasia. Nature. 522:167–172.
    doi:10.1038/nature14507.
Ameur A, et al. 2017. SweGen: a whole-genome data resource of genetic variability in a
    cross-section of the Swedish population. Eur J Hum Genet. 25:1253–1260. doi:10.1038/ejhg.2017.130.
```

The archived GSA page gives an older comma-and-initials variant
(`Texada, M. J., R. A. Simonette, C. B. Johnson, W. J. Deery, and K. M. Beckingham, 2008 Yuri gagarin
is required for…  J. Cell Sci. 121: 1926-1936.`) `[GSA-A]`. **Follow the `[G3-2026]` form** — that is
what the journal prints now — and let copyediting reconcile any residual difference.

Rules that hold in both: `[GSA-A]`
- Order alphabetically by first author.
- More than five authors → list the first five, then `et al.`
- Same first author + same year → `2000a`, `2000b`, alphabetised by second author.
- Books: `Sturtevant AH, Beadle GW. 1939. An Introduction to Genetics. Philadelphia: W. B. Saunders.`
- Book chapters: author, year, title, `pp. X–Y in <Book>`, editors, publisher, city.
- Preprints are citable, with DOI and posting date:
  `… bioRxiv. doi:10.1101/123456 (Preprint posted March 1, 2016).`
  → Applies to several of our own citations (the 2026a/2026b preprints, refs 30, 65, 66, 84).
- Include DOIs. `[G3-2026]`

---

## 9. Data availability, code, and supplementary material

### 9.1 Policy

- G3 requires authors to **publicly release all data and software code underlying the paper as a
  condition of publication**, and to include a **data availability statement**. `[SEARCH]`
- The statement must **list accession numbers or DOIs** of anything in a public repository, **list the
  file names and descriptions of all supplemental files**, and include any applicable IRB numbers.
  `[GSA-A]`
- Accession numbers/DOIs must exist **before publication**. `[GSA-A]`
- `[LETTER]`: *"the GitHub page link provided in the data availability statement is broken; please
  ensure that the data availability statement is updated in your revision."*
- `[LETTER]`, Associate Editor: *"Please ensure that code/data links are made available in the
  revision"* and *"it would be very helpful to have a browser/track instance for the TE annotations
  made available with the GitHub repository."*

### 9.2 Supplementary material naming — second formatting gap

G3 naming convention: **`Figure S1`, `Table S1`, `File S1`** — numbered figures/tables/files, with:
`[GSA-A]`

- supplemental **figure** title and legend placed **below** the figure; title begins with
  `Figure S1` in **bold**;
- supplemental **table** title placed **above** the table, beginning `Table S1`, bold;
- **files** = anything that is not a figure or table (extra methods, datasets, code) → `File S1`…;
- tables preferred as **Excel, .csv or .txt**, not PDF;
- every supplemental item **cited at least once in the manuscript** as "Figure S1", "Table S1",
  "File S1";
- all supplemental items summarised in the Data availability statement;
- **upload to the GSA figshare portal** (`https://gsajournals.figshare.com/submit`) in a single batch,
  titled `Supplemental Material for Nikitin 2026`; individual file titles cannot be changed after
  upload. A repository of the author's choosing is also acceptable. `[GSA-A]` `[VERIFY]` — confirm
  whether OUP-era G3 still routes supplements through figshare or now takes them in the submission
  system.
- If the manuscript describes software, include a **code freeze**. `[GSA-A]` → our Zenodo snapshot
  (plan decision D9) serves this purpose.

**Our draft uses "Supplementary Figure 1–8" and "Supplementary File 1–8" throughout the text, the
captions and the file names on disk.** All of it must be renamed to `Figure S1`–`S8` and
`File S1`–`S8` (plus the new supplementary items added by the revision), in: the manuscript body, the
figure/file captions, the on-disk filenames, `CLAUDE.md`, and the response letter.

### 9.3 A model statement for our paper

```
Data availability

All code used for proximity mapping, enrichment statistics, permutation testing, and Gene Ontology
analysis is available at https://github.com/Nikit357/T2T_transposons_genes and archived at Zenodo
(doi:10.5281/zenodo.XXXXXXX). Transposable element annotations were obtained from the UCSC
RepeatMasker track for T2T-CHM13v2.0 (hs1) and gene annotations from NCBI RefSeq All (hs1) via the
UCSC Table Browser. A UCSC Genome Browser track hub containing all TE annotations used here, the
10 kb TSS neighbourhoods, and the gene sets analysed is available at
https://nikit357.github.io/T2T_transposons_genes/trackhub/hub.txt and can be loaded directly at
<one-click URL>. File S1 contains … [enumerate every File S1–Sn with a one-line description, including
which files carry raw as well as FDR-adjusted p-values]. This study analysed only publicly available
genomic data; no IRB approval was required.
```

---

## 10. Nomenclature, italics, numbers, units

- **Genetic nomenclature**: use the standard for the organism, consulting the organism database — for
  human genes, HGNC-approved symbols. `[GSA-A]` Italicise genotype names and symbols; do **not**
  italicise phenotypes; distinguish genotype from phenotype consistently. `[GSA-A]`
  → Check every human gene symbol in our text (`IFNA10`, `IFNW1`, `POLR2A`, `SSU72L1–L5`, `CIB3`,
  `TMC1/2`, `XIST`) against HGNC and apply consistent italics. `[VERIFY]` the exact human-gene italics
  convention G3 applies.
- Italicise organism names even without the species epithet. `[GSA-A]`
- Numbers: write out nine and below, except in dates, fractions/decimals, percentages and units of
  measurement; Arabic numerals above nine; avoid beginning a sentence with a numeral. `[GSA-A]`
- Units: abbreviate only after a number ("3 min", but "several minutes"). "percent" spelled out except
  with a numeral ("75%"). Temperature as `37°`. `[GSA-A]`
- **RRIDs** encouraged for reagents/organisms/tools. `[GSA-A]` For a purely computational paper the
  equivalent good practice is software versions with citations — our Methods already lists versions for
  bedtools, scipy, statsmodels, goatools, matplotlib, seaborn, plotly, networkx, pyvis, supervenn,
  pandas, numpy. Keep that.
- **Reagent Table**: required/encouraged for wet-lab reagents. `[GSA-A]`, `[SEARCH]` Not applicable to
  this study — state that explicitly in the cover letter so it does not read as an omission. `[VERIFY]`
- In-article database linking (e.g. FungiDB gene links via Word's Hyperlink style) — not applicable to
  human genes, but indicate on the submission form whether the manuscript contains linkable entities.
  `[GSA-A]`
- Mathematical characters: use MathType 5.0+ for symbols not on the keyboard; characters typed as text
  can cause processing errors. `[GSA-A]` Relevant to our `<10⁻²⁰⁰`, `9.3*10⁻¹³³` table entries — convert
  to proper superscripts/scientific notation rather than asterisk-times constructions.

---

## 11. Submission package for this revision

`[LETTER]` requires exactly three items, plus the formatted files:

1. **Clean version** of the manuscript, formatted for G3 (everything in this document).
2. **Highlighted or tracked-changes version** linking each response to the changed text.
3. **A separate response document** addressing each editor/reviewer comment.

Plus `[GSA-A]`: at revision, provide the **non-PDF source files** — manuscript in Word or LaTeX,
each figure as a separate vector/high-resolution file, tables editable at the end of the main text,
and supplementary material uploaded to figshare (or the chosen repository).

Submission link (`[LETTER]`):
`https://g3.msubmit.net/cgi-bin/main.plex?el=A7NQ1IFP2A1gOY1I4A9ftdolTB6k66m5lMLcBNIDNmgZ`

---

## 12. Gap list — what our current draft violates

Ordered by effort. Every item is also in the plan's Phase 7 checklist. Figure-side items (G11, G12, G17,
G18) are executed **by Daniil, manually, in the Figma source file** — the notebooks emit panel SVGs only
and nothing writes to Figma. See `G3_revision_implementation_plan_260803.md` §5b for the file key, the
frame-per-figure node-ID map, and the manual placement workflow, and `revision_G3/svg/PLACEMENT.md` for
which SVG goes into which frame.

**Ownership:** G1 (citation style) is **Daniil's**, done as a Mendeley Cite style switch — see §8.1.
Everything else in this list is ours, executed as Word tracked changes.

| # | Gap | Where | Fix |
|---|---|---|---|
| G1 | **Numbered citations** instead of author–year; several in-sentence uses (*"the study by (14)"*) are ungrammatical under author–year | whole manuscript, 123 refs | swap to CSE author–year via Zotero/EndNote G3 style file; manually rewrite in-sentence citations |
| G2 | **Supplementary items named "Supplementary Figure/File N"** instead of `Figure Sn` / `File Sn` | body text, captions, on-disk filenames, `CLAUDE.md` | global rename; re-cite each in text |
| G3 | **Broken repository URL** in Data availability (`T2T_genes_evolution`) | line 371 of the draft | correct URL + Zenodo DOI + track hub URL; enumerate every File Sn |
| G4 | **No keywords** | title page | add 3–10 (§4) |
| G5 | **No short running title** (~35 chars) | title page | add |
| G6 | **Incomplete corresponding-author block** (email only) | title page | full institutional address, phone, ORCID |
| G7 | **Section headings** in all-caps and non-G3 order; `REFERENCES` instead of `Literature cited`; `ETHICAL STATEMENT`/`CONFLICT OF INTEREST`/`AUTHORS CONTRIBUTION`/`ACKNOWLEDGEMENTS` ordering and spelling | whole manuscript | reorder per §2; `Acknowledgments` (US spelling), `Conflicts of interest` (plural), `Data availability`, `Funding`, `Literature cited` |
| G8 | **Table 1 is 11 columns wide** (also Reviewer minor 4) | Table 1 | split into Tables 1 and 2, not 1A/1B (§7) |
| G9 | **Scientific notation typed as `9.3*10⁻¹³³`** | Table 1 | proper notation (§10) |
| G10 | **p-value provenance unlabelled** in several figures (also Reviewer minor 3, and journal policy §5) | Figures 1D, 3, 7; captions | label raw vs FDR-adjusted everywhere; state the raw-p-value availability in Methods |
| G11 | **Figure suptitles baked into images** (e.g. "Detailed TE Family Enrichment Analysis") | Figure 1D and others | remove suptitles; titles belong in the manuscript legend |
| G12 | **Raster elements below print resolution** — pasted colorbars (23 × 58 … 30 × 175 px) and schematic bitmaps; the UCSC panels are vector and fine | Figures 4, 5, 6, 8, 9, S8 | re-emit colorbars as vector; verify schematic bitmaps ≥ 350 dpi |
| G17 | **Analysis parameters baked into figure text** — Figure 7 hard-codes `GO terms count in a group (FDR < 0.1)` | Figma frame `861:34` | edit to 0.05 after the FDR change; grep all frames for `0.1` / `500` / `p-value` |
| G18 | **Mixed fonts and axis-label sizes** — Inter everywhere except Helvetica in the vector UCSC panels; label sizes 10 / 13.3 / 14.7 / 16 px across frames | Figures 1, S6 and others | unify while the frames are open |
| G13 | **Human gene symbols not consistently italicised** | throughout | HGNC symbols + consistent italics |
| G14 | **Acknowledgments are unusually long** (4 paragraphs, ~600 words, personal) | Acknowledgments | not a rule violation, but trim toward journal norm; keep the funding statement separate under `Funding` |
| G15 | **Abstract does not name the organism explicitly** as required | Abstract | ensure "human" appears in the first sentences (it does — verify against §4 wording) |
| G16 | **Ethical statement placement** | after Methods | G3 has no `ETHICAL STATEMENT` section; fold the "no human subjects / public data only" sentence into Data availability or Methods |

---

## 13. Verify on the live page before submitting `[VERIFY]`

The editor warned the guidelines were recently updated. Open
`https://academic.oup.com/g3journal/pages/author-guidelines#section-11` in a browser and confirm:

1. Whether the Data availability statement belongs at the end of Materials and methods, in the back
   matter, or both.
2. Whether supplementary material still goes to the GSA figshare portal, or now into the submission
   system / another repository.
3. The exact current Literature cited format (initials spacing, `et al.` punctuation, author cut-off).
4. Current figure file formats and minimum resolutions (the 350/600/1200 dpi figures above come from
   the archived GSA page).
5. Whether an AI-usage disclosure section is required and where it goes — our draft has an
   `AI usage` subsection in Methods (Gemini Pro for code refinement, literature search and GO
   pre-classification; ChatGPT for grammar). OUP has a generative-AI policy; confirm the required
   wording and placement. **Do not delete this disclosure.**
6. Whether a Reagent Table is expected even for purely computational papers.
7. Any structured-abstract, graphical-abstract or plain-language-summary requirement.
8. ORCID requirement for the corresponding author.
9. Word/figure/table caps for Investigations (archived source says none).
10. Whether `Code availability` should be its own back-matter section (one 2026 paper used it) — likely
    a good idea for us given the Associate Editor's request.

---

## 14. Sources

- [G3 author guidelines (live page — HTTP 403 from this environment)](https://academic.oup.com/g3journal/pages/author-guidelines)
- [GSA "Preparing Manuscripts for Submission" (Internet Archive capture)](https://web.archive.org/web/2021/https://www.g3journal.org/content/prep-manuscript)
- [GSA "Article Types" (Internet Archive capture)](https://web.archive.org/web/2021/https://www.g3journal.org/content/article-types)
- [Overleaf/GSA LaTeX template for G3](https://www.overleaf.com/latex/templates/template-for-preparing-your-submission-to-g3-genes-genomes-genetics-using-overleaf/vffkrpmjrcgf)
- Current house style extracted from G3 articles in PubMed Central: PMC13365844, PMC13334165,
  PMC13334167, PMC13334169, PMC13334170, PMC13334186
- Decision letter G3-2026-406828, July 29 2026 (`G3_reviewer_report_260802.md`)
- Figure source file: [Figures T2T genes (Figma)](https://www.figma.com/design/WRNeTzKZObdmAQ8QG1EZlq/Figures-T2T-genes)
  — inspected read-only 2026-08-03; frame-per-figure map in the implementation plan §5b
