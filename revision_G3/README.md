# `revision_G3/` — revision package for manuscript G3-2026-406828

Everything the G3 revision changes lives in this directory. It is **authoritative** for the
revision: where a helper exists both here and in one of the parent directory's notebooks, the copy
here is the one the revised manuscript reports.

- **Manuscript:** "Telomere-to-telomere co-mapping of transposable elements and human genes
  identifies a cluster of young L1 elements in the interferon-alpha domain" (retitled per D3).
- **Journal:** G3: Genes|Genomes|Genetics — conditionally accepted, minor revisions.
- **Plan:** `../G3_revision_implementation_plan_260803.md` (work packages WP1–WP16, decisions D1–D9).
- **Reviewer report:** `../G3_reviewer_report_260802.md`.
- **House style:** `../G3_article_guidelines.md` (gap list G1–G18).

---

## 1. The frozen-notebook rule

Four files in the parent directory are **frozen for this revision — zero edits, not a cell, not a
label, not a threshold**:

```
6d59a2a735b8d0f4fcf6d9dddbb8bb39  TEs_mapped_on_TSS_analysis.ipynb
a75ceaf51c0a0d221f53357bb0040b55  Gene_ontology_analysis.ipynb
3e8aec87bd9e78fce53463a2073d968b  download_and_process_files_UCSC_genes.ipynb
cfd78a7eb38b8f5bbc76dd0fba75dc01  GO_subfamilies.py
```

These MD5s are the freeze baseline. All four are committed to git as of the first commit on branch
`g3-revision` ("Freeze baseline: commit the three analysis notebooks as-is"), so the freeze is
verifiable two ways:

```bash
md5sum -c <<'SUMS'
6d59a2a735b8d0f4fcf6d9dddbb8bb39  TEs_mapped_on_TSS_analysis.ipynb
a75ceaf51c0a0d221f53357bb0040b55  Gene_ontology_analysis.ipynb
3e8aec87bd9e78fce53463a2073d968b  download_and_process_files_UCSC_genes.ipynb
cfd78a7eb38b8f5bbc76dd0fba75dc01  GO_subfamilies.py
SUMS
git status --short -- TEs_mapped_on_TSS_analysis.ipynb Gene_ontology_analysis.ipynb \
  download_and_process_files_UCSC_genes.ipynb GO_subfamilies.py   # expect empty
```

Consequences, all deliberate:

1. Helper functions the old notebooks defined in-cell were **copied** into `revision_lib.py` and
   modified there (§3 below). Two versions of the same logic now exist at two FDR thresholds, and the
   older one still runs. For anything G3-related, `revision_lib.py` is authoritative.
2. **Never import from a frozen notebook** (no `nbimporter`, no `import_ipynb`). That would
   reintroduce the coupling the freeze removes, and would silently pick up FDR 0.1.
3. `GO_subfamilies.py` stays at `fdr_threshold=0.1` because it belongs to the companion subfamilies
   manuscript, whose refresh is deferred. The revision's subfamily-level GO at 0.05 comes from
   `06_go_rerun_fdr005.py --level subfamilies` and writes here, leaving the shared `../GO_tables/`
   untouched.

## 2. N = 500 permutations — the authoritative number

**500 permutations were run.** `../epigenomic_files/permutation_results/` held 501 files (500
per-seed BEDs plus the consolidated CSV), and Table 1's empirical p = 0.004 is exactly
2/(500+1) = 0.003992, which confirms it independently. The Results section and Figure 1D already said
500; the **Methods said 1,000 and was wrong**, and is corrected **down** to 500 in this revision.

The generator cell in the frozen `download_and_process_files_UCSC_genes.ipynb` (cell 34) still reads
`NUM_PERMUTATIONS = 1000`. **That value was never executed.** It is the origin of the manuscript's
inconsistency, it stays in the notebook because the notebook is frozen, and it is documented here, in
`../CLAUDE.md` and in `../REPRODUCE.md` rather than fixed. A reader who opens that notebook will see
`1000` while the paper says 500; the explanation is the one in this paragraph.

The permutation background exists to **correct the odds ratio for element-length bias** (random OR
vs. mean element length, **Pearson R = 0.985** across the 44 families, n = 44, p = 8.4 × 10⁻³⁴,
mean lengths 122–6,357 bp; Alu 316 bp → random OR 1.535, L1 6,357 bp → 2.664, persisted to
`output/length_bias_correlation.json`. An earlier version of this file said 0.661, which was wrong),
not to serve as the primary significance test — that is
the FDR-adjusted Fisher exact p. What N has to deliver is a stable mean and SD of the random OR, not
a small p-value floor, and `nb01_permutation_convergence.ipynb` demonstrates that stability is
reached well before 500.

## 3. `revision_lib.py` — provenance of every copied helper

| Function in `revision_lib.py` | Copied from | Changed how |
|---|---|---|
| `run_goatools_enrichment` | `Gene_ontology_analysis.ipynb` cell 6 | `fdr_threshold` 0.1 → **0.05** (D2); cached DAG/GAF |
| `run_goatools_ordered_enrichment` | `Gene_ontology_analysis.ipynb` cell 6 | same |
| `save_go_network_svg` | `Gene_ontology_analysis.ipynb` cell 6 | `+fdr_threshold`, `+min_shared_genes`, `+max_term_genes`, `+title` (suppressible, G11), enforced collision check |
| `visualize_go_class_network` | `Gene_ontology_analysis.ipynb` cell 36 | same filters; class palette switched from Tableau to the shared TE class palette |
| `save_go_network_svg_families_by_classes` | `Gene_ontology_analysis.ipynb` cell 175 | same filters; `family_to_class` is now an explicit argument instead of a notebook global |
| `CLASS_PALETTE` | `Gene_ontology_analysis.ipynb` cell 3 (`class_names_p`) | verbatim |

New in this module, with no frozen original:

| Function | Purpose |
|---|---|
| `assert_no_label_collisions(fig)` / `find_label_collisions` / `save_svg_collision_checked` | Bounding-box label-overlap check, **enforced** before every SVG write (reviewer minor comment 6, caveat C16). A colliding figure cannot be saved silently. |
| `load_permutation_counts(window=...)` | Reader for the compacted permutation store, replacing the 6.37 GB `consolidated_random_data.csv` (WP1b). `n` is a **weight** column, not a row count. |
| `permutation_totals(window, by=...)` | Per-permutation totals by class or family, streamed one seed at a time. |
| `read_counts_file` / `load_manifest` / `permutation_store_dir` | Store plumbing. |

**One change beyond FDR and the new arguments, and why it was necessary:** the frozen
`run_goatools_enrichment` parses the 31 MB GO DAG and the 190 MB GAF *inside every call*. The
revision's GO re-run makes roughly 1,200 calls (6 classes + 2 divergence groups + 44 families +
1,143 subfamilies), so the copies here load both once per process and cache them. The GO results are
unchanged — the same DAG and the same association dict reach `GOEnrichmentStudy` either way — and the
positive control (300 olfactory receptor genes → `olfactory receptor activity`, FDR < 0.05) passes.

## 4. Run order

Two ordering constraints that are not obvious:

* **`nb07` before `nb05`.** `nb05`'s panel S13D and its headline-claim table both read
  `output/go_grid_headline_by_condition.csv`, which `nb07` writes. `nb05` is the single writer of
  `robustness_headline_claims_by_condition.csv` — if both notebooks wrote it, whichever ran last
  would decide its contents and the `not evaluated` cells the GO grid was built to fill could
  silently come back.
* **`07a` before `07b`.** `07a --verify-10kb` is a gate, not a formality: it refuses to write the
  5 kb and 20 kb gene tables unless rebuilding the 10 kb one reproduces the published
  `TEs_on_genes.csv` exactly (row count, gene count, empty windows, per-class totals, mean
  divergence to 1e-9, and the per-window multiset of family / class / subfamily names).

Compute lives in `.py` / `.sh` (cacheable, re-runnable, backgroundable); every figure lives in a
notebook so each subfigure can be inspected and adjusted before it goes anywhere. **Nothing in this
directory writes to Figma.** Notebooks emit SVG panels to `svg/`, and `svg/PLACEMENT.md` maps each
one to its target Figma frame for the manual placement step.

```bash
source ~/venvs/Retroelements_3_11/bin/activate
cd /home/jovyan/Projects/Retroelements/T2T_genes_article/T2T_transposons_genes

# --- Phase 1: free the disk first, then launch the long jobs -----------------
python revision_G3/01b_compact_permutation_results.py            # 11 GB -> ~110 MB   ~1-2 h
python revision_G3/01c_expand_counts.py --seed 1 --check         # lossless proof
bash   revision_G3/05a_build_windows.sh                          # 5 kb / 20 kb windows
bash   revision_G3/01_permutations_stream.sh --window 5kb        # N=500             2-4 h
bash   revision_G3/01_permutations_stream.sh --window 20kb       # N=500             2-4 h
python revision_G3/01a_consolidate_counts.py --check-legacy      # regression vs legacy 10 kb
python revision_G3/06_go_rerun_fdr005.py --level all             # FDR 0.05          1-2 h

# --- Phase 2: analysis ------------------------------------------------------
python revision_G3/02_ifna_domain_test.py                        # 4 tests           ~20 min
python revision_G3/04a_lu2020_geneset_overlap.py                 # Lu et al. overlap
python revision_G3/04b_newly_resolved_regions.py                 # assembly bound (C6)   ~10 min
python revision_G3/05b_window_sensitivity.py                     # 6 classes + 44 families x 3 windows
python revision_G3/05c_percentile_sensitivity.py                 # GO at 10 % vs 5 %  ~10 min
python revision_G3/10_tables.py                                  # Table1.csv + Table2.csv

# --- Phase 4: figures (notebooks, SVG only) ---------------------------------
jupyter notebook revision_G3/nb01_permutation_convergence.ipynb  # -> svg/S14_*
jupyter notebook revision_G3/nb02_ifna_domain.ipynb              # -> svg/Fig8B_*, Fig8C_*
jupyter notebook revision_G3/nb03_relabelled_figures.ipynb       # -> svg/Fig1D_*, 2D-2F, 3A/3B, S1-S3
jupyter notebook revision_G3/nb06_go_networks_fdr005.ipynb       # -> svg/Fig4A_*, 5A_*, 6A_*, S8B/S8C, S9-S11
jupyter notebook revision_G3/nb07_go_grid_robustness.ipynb       # -> svg/S12A-C   (BEFORE nb05)
jupyter notebook revision_G3/nb05_sensitivity_robustness.ipynb   # -> svg/S13_*

# --- Phase 5: manuscript (needs python-docx: use the collagen venv) ---------
python revision_G3/11_results_numbers.py                         # every number the text quotes
~/venvs/collagen_3_11/bin/python revision_G3/13_manuscript_tracked_edits.py
~/venvs/collagen_3_11/bin/python revision_G3/15_house_style.py      # MUST follow 13, not replace it
~/venvs/collagen_3_11/bin/python revision_G3/14_build_extensive_discussion.py

# --- Phase 6: deliverable ---------------------------------------------------
bash revision_G3/12_build_trackhub.sh                            # bigBeds + hub, ~1 h
bash revision_G3/12b_publish_trackhub.sh --push                  # -> origin/gh-pages
bash revision_G3/12c_verify_trackhub_live.sh                     # only after Pages is enabled
```

**The hub is live only after GitHub Pages is enabled**, which is a one-time manual step on the
`gh-pages` branch that only the repository owner can take (*Settings → Pages → Deploy from a branch
→ `gh-pages`, `/ (root)`*). Until then the URL the manuscript prints returns 404 no matter how many
times the branch is pushed — that was the state of things on 2026-08-14, diagnosed in
`../trackhub_ghpages_plan_260814.md`.

**The working manuscript changed on 2026-08-04, and with it how these two scripts behave.** Daniil
edited the revision by hand and **accepted 594 of the 1,124 tracked deletions**, so
`Revised_manuscript/T2T_genes_article_G3_revision_260804_manual.docx` — not the baseline, and not the
260803 script output — is the current state of the paper. Both `13_…` and `15_…` target it.

`13_manuscript_tracked_edits.py` used to read the untouched baseline and overwrite the revision, and
that copy *was* its idempotency mechanism. Doing that now would throw away every acceptance and
manual edit, so it **edits the working file in place** and is idempotent by per-edit detection
instead: each edit reports `applied`, `already present` (skipped), or **`not found`**, and only the
last is fatal — so a re-run is no longer destructive while a silent skip is still impossible. Every
structural insertion (the two tables, the new Discussion subsections, the repository-URL paragraph)
carries a marker sentence it is checked against, because unlike a text replacement an insertion is
not self-limiting: without the marker a second run would add a second copy of the same 16 paragraphs.
`--from-baseline` still exists for a clean rebuild onto a *different* filename, and refuses to
overwrite the working file.

`15_house_style.py` edits the same file, refuses to start if the Phase 5 edits are missing, and
refuses to run twice. On the working file it correctly declines: the house style is already there.
The old rule "15 must run after 13, every time" is therefore **superseded** — 13 no longer destroys
the Phase 7 pass. What still needs solving is that the M1–M10 text has not been written yet, and when
it is, its new `Supplementary Figure`/`File` strings will need the G2 sweep re-run over the newly
inserted paragraphs only.

**`Reject All` no longer restores the baseline.** Accepted deletions are permanent by design, so the
clean-plus-tracked submission pair comes from the working file forward and the diff a reviewer sees
covers the edits made *after* the acceptance pass. If the journal wants a full tracked diff against
the submitted version, produce it as a Word **Compare** of the baseline against the final clean file
rather than from the revision marks.

All three manuscript scripts need `python-docx` + `lxml`, which are in `~/venvs/collagen_3_11` and
**not** in `~/venvs/Retroelements_3_11` (caveat C15).

**Do not start the 5 kb / 20 kb permutations before the compaction finishes** — free disk was 2.4 GB
before WP1b and the two new backgrounds need the space it frees.

## 5. Canonical input files

Use these and only these for TE ↔ TSS intersections. They supersede the broken
`../../T2T_article/T2T_repeat_masker_processed.csv` path and the missing `T2T_genes_sorted.bed`
recorded in `../CLAUDE.md`.

| File | Rows | Columns | Why it is safe to use |
|---|---|---|---|
| `../T2T_repeat_masker_processed_sorted.bed` | 3,709,429 | `chrom, start, end, score(divergence), subfamily, family, class` | Byte-size-identical to `../../epigenomic_files/repeats_all.bed`, the file the legacy permutations shuffled |
| `../T2T_genes.bed` | 38,704 | `chrom, start, end, gene` | Interval set exactly identical to the `mapped_on_TSS` bedGraph used as the legacy `-b` file |

Because both are provably interchangeable with what the legacy pipeline consumed, the retained 10 kb
N = 500 background and the new 5 kb / 20 kb runs are directly comparable — there is no hidden
geometry change to explain in the Methods.

WP4 additionally needs three **third-party** inputs, which live in `external/` and are the only files
here that this repository does not produce: Lu et al.'s supplementary Table S1, the MGI mouse–human
homology report, and the UCSC T2T → GRCh38 liftOver chains. `external/PROVENANCE.md` records the exact
URL, size and checksum of each, what it contains, and the four download routes that no longer work —
including the discovery that Lu et al.'s categories are **mouse** gene sets, which is why WP4 has to
map through orthology.

## 6. Directory layout

```
revision_G3/
  README.md                          this file
  revision_lib.py                    shared helpers (§3)
  01_permutations_stream.sh          N=500 streaming permutations, one window per invocation
  01a_consolidate_counts.py          consolidated counts + legacy regression check
  01b_compact_permutation_results.py 11 GB -> ~110 MB, verified per seed before any delete
  01c_expand_counts.py               reconstruct a legacy BED from counts on demand
  02_ifna_domain_test.py             four IFNA-domain tests
  04a_lu2020_geneset_overlap.py      gene-set overlap with Lu et al. 2020
  04b_newly_resolved_regions.py      upper bound on the assembly-attributable component (C6)
  05a_build_windows.sh               5 kb / 20 kb TSS neighbourhoods + TE mapping
  05b_window_sensitivity.py          enrichment at 5/10/20 kb
  05c_percentile_sensitivity.py      GO at top/bottom 10 % vs 5 %
  06_go_rerun_fdr005.py              every GO level at FDR 0.05
  10_tables.py                       Table1.csv + Table2.csv (reformat; values unchanged)
  11_results_numbers.py              re-derives every number the revised manuscript quotes
  12_build_trackhub.sh               bigBeds + hub.txt/genomes.txt/trackDb.txt + hubCheck
  12b_publish_trackhub.sh            publish trackhub/ to gh-pages; the push is behind --push
  12c_verify_trackhub_live.sh        live hub check: 200, 206, bigBed magic, sizes vs MANIFEST.json
  13_manuscript_tracked_edits.py     Phase 5 D-K, all as tracked changes; citation-safe
  15_house_style.py                  Phase 7 G3 house style; MUST run after 13
  14_build_extensive_discussion.py   Extensive_discussion_260803.docx by copy-and-delete
  Revised_manuscript/                the three manuscript docx files, kept together
    T2T_genes_article_G3_submitted_baseline_260418.docx   read-only; what the journal received
    T2T_genes_article_G3_revision_260803.docx             the revision, tracked changes
    Extensive_discussion_260803.docx                      WP8/D7 supplementary file
  nb01_permutation_convergence.ipynb WP1 justification figure (S14)
  nb02_ifna_domain.ipynb             WP2 figures + the four test results inline
  nb03_relabelled_figures.ipynb      WP9 label-only re-creations; no numbers change
  nb05_sensitivity_robustness.ipynb  WP5 figures + the robustness comparison (run AFTER nb07)
  nb06_go_networks_fdr005.ipynb      WP6/WP12 networks and heatmaps at FDR 0.05
  nb07_go_grid_robustness.ipynb      the 3 windows x 2 percentiles GO grid comparison (S12A-C)
  07a_build_gene_tables.py           per-window TE-TSS tables; --verify-10kb is the gate
  07b_go_grid.py                     GO across the grid; reuses the two published 10 kb cells
  08_build_supplementary.py          the five thematic Excel workbooks
  supplementary/                     the deliverable: Files S1-S5 + README + CHECKSUMS, 8.7 MB
  external/                          third-party downloads for WP4; see PROVENANCE.md
    PROVENANCE.md                    exact URLs, sizes, md5s, and the routes that fail
    lu2020/mmc2..7.xlsx              Lu et al. 2020 supplementary tables (their Table S1 = mmc2)
    HOM_MouseHumanSequence.rpt       MGI mouse-human homology (gitignored, 15 MB)
    hs1ToHg38.over.chain.gz          UCSC T2T -> GRCh38 liftOver chains
  output/                            every derived table
    legacy_fdr01_n500/               pre-revision snapshot (caveat C3)
    permutation_counts_10kb/         compacted legacy store + MANIFEST.json
    permutation_counts_5kb/          new N=500 background
    permutation_counts_20kb/         new N=500 background
  svg/                               one SVG per panel, for manual Figma placement
    PLACEMENT.md                     SVG -> Figma frame ID + panel position
  trackhub/                          the built hs1 hub, 105 MB, gitignored; published to gh-pages
```

The repository root additionally carries `.github/workflows/`: `verify-trackhub.yml` (manual
dispatch plus a weekly link-rot check of the published hub) and `deploy-trackhub.yml` (dormant;
only needed if Pages is ever switched to the "GitHub Actions" source). Neither builds the hub —
building needs the 155 MB gitignored RepeatMasker BED and UCSC's `bedToBigBed`.

## 7. Environment decisions

| Item | Decision |
|---|---|
| Analysis / plotting venv | `~/venvs/Retroelements_3_11` — has `goatools`, `adjustText` 1.3.0, `networkx`, `pyvis`, `supervenn`, `statannotations` |
| **docx / tracked-changes venv (C15)** | **`~/venvs/collagen_3_11`** — it is the only venv with `python-docx` 1.2.0 + `lxml`, which `~/.claude/skills/word_rewrite_trackchanges.py` needs. Nothing was installed into the Retroelements venv, so the two stay independent. |
| `bedtools` | 2.31.1 at `/usr/local/bin/bedtools` |
| `zstd` | 1.5.7 at `/opt/conda/bin/zstd` |
| Cores | 8 → `xargs -P 8` |

`adjustText` in this venv is **1.3.0**, which renamed the layout kwargs: `force_points` →
`force_static` and `expand_text` → `expand`. The frozen notebooks only ever passed `arrowprops`, so
they were unaffected; new code must use the 1.3.0 names.

## 8. Open items to confirm with Daniil

1. **`go-basic.obo` release date.** The Methods cite a December 31 2025 ontology snapshot, but the
   local `go-basic.obo` self-reports `rel(2025-10-10)` and both it and `goa_human.gaf` are dated
   October 17 2025 on disk. These are the files every GO number in the paper came from, so the
   Methods date needs correcting to match the files rather than the files being re-downloaded —
   re-downloading would change term memberships and invalidate the measured FDR shifts in plan §3.2
   (caveat C5).
2. **The `NUM_PERMUTATIONS = 1000` line in the frozen download notebook** (§2 above). It is a
   one-line fix and the only exception to the freeze worth arguing for; left as-is and documented
   until Daniil says otherwise.
3. **The submitted baseline docx.** `T2T_genes_article_for_plos_one.docx` was renamed to
   `T2T_genes_article_G3_submitted_baseline_260418.docx` and made read-only; the working file is
   `T2T_genes_article_G3_revision_260803.docx`. Both were originally filed with the subfamilies
   paper and now live in `Revised_manuscript/` (below). The baseline is byte-identical to the
   submitted file (md5 `1dbcbd4419987fd28ddf803129487cfd`). Nothing was deleted.
