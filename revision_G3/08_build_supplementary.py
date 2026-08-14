#!/usr/bin/env python
"""Assemble the supplementary package as five thematic Excel workbooks (decision D-d).

Why this script exists
----------------------
There was no assembled supplementary package anywhere in the repository. The only
supplementary files on disk were the April baseline in
`../manuscript_figures_supplementary_old/`, whose GO-derived members are still at
FDR 0.1 and are therefore superseded by decision D2; the revision's own new tables
were loose CSVs in `output/` with no `File Sn` number and no directory. Assembling
by script rather than by hand is what makes the package reproducible and what lets
`--verify` state, rather than assume, that every GO sheet is at 0.05.

Five workbooks, not fourteen files (decision D-c / D-d)
-----------------------------------------------------
A reader opens five things instead of fourteen, and each workbook is one subject
with one sheet per table. The gene-set sheets are **long format** — the old
Files S3 / S5 / S7 were 8, 10 and 44 single-column sheets, and reproducing 44
sheets inside one workbook would defeat the point of the regrouping. It also fixes
a real defect: Excel forbids `?` in a sheet name, so the old File S7 shipped the
`hAT?` and `hAT-Tip100?` families under the corrupted names `hAT...` and
`hAT-Tip100...`. In long format the family name is a cell value and survives intact.

What is carried over unchanged, and what is re-emitted
-----------------------------------------------------
* Gene sets are the **published** sets, read from the April workbooks. They are
  gene lists, not GO output, so the FDR change does not touch them. Note those
  sheets store the first gene in the header row, so they are read with
  `header=None`.
* Every GO sheet is **re-emitted at FDR 0.05** from `output/GO_*_fdr005.csv` and
  asserted to satisfy `FDR.max() < 0.05` before it is written. The April copies are
  at 0.1 and are not used.
* `TSS_TE_intersections` is the published `../TEs_on_genes.csv` reduced to exactly
  the columns the April `Supplementary File 1.csv` shipped — i.e. without the two
  epigenomic `Signal` / `Signals` columns, which that file already omitted.

The old File S2 mismatch, resolved by carrying both tables
---------------------------------------------------------
The manuscript describes File S2 as "enrichment statistics of TE subfamilies", but
the April `Supplementary File 2.csv` contains the **44 families**. (It does carry the
permutation columns, contrary to what the implementation plan recorded; the caption
mismatch is real, the missing-columns half of it was not.) Workbook S1 now ships
`enrichment_families` **and** `enrichment_subfamilies`, both with the full column
set including raw and adjusted p-values, so both the caption and the Data
availability claim become true without dropping anything.

Excel's limits are enforced, not hoped for
------------------------------------------
Sheet names are capped at 31 characters and may not contain ``[]:*?/\\``; a cell may
not exceed 32,767 characters. `Full Term Gene List` reaches 125,086 characters — the
full human annotation of "protein binding" — so it is dropped from the GO sheets and
the drop is declared in each workbook's README. It is a property of the GO
annotation (`go-basic.obo` + `goa_human.gaf`, both cited in the Methods), not a
result of this study, and `Overlapping Genes` — which *is* the result — is kept. Any
other over-limit cell is a hard error rather than a silent truncation.

Figures stay placeholders
-------------------------
The figures are exported from Figma by hand (Phase 4b), so the script writes
`figures/PLACEHOLDERS.md` listing the 14 expected filenames rather than inventing
content.

Usage
    python revision_G3/08_build_supplementary.py            # build all five
    python revision_G3/08_build_supplementary.py --verify   # check what is on disk
"""

from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import os
import re
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import revision_lib as rl  # noqa: E402

OUTPUT_DIR = rl.OUTPUT_DIR
GRID_DIR = os.path.join(OUTPUT_DIR, "GO_grid")
SUPP_DIR = os.path.join(rl.REVISION_DIR, "supplementary")
OLD_DIR = os.path.join(rl.REPO_DIR, "manuscript_figures_supplementary_old")

EXCEL_SHEET_NAME_LIMIT = 31
EXCEL_FORBIDDEN_IN_SHEET_NAME = re.compile(r"[\[\]:*?/\\]")
EXCEL_CELL_LIMIT = 32_767
EXCEL_ROW_LIMIT = 1_048_576

# Columns dropped because a single cell exceeds Excel's limit. Each needs a reason
# that says why losing it from the workbook loses nothing.
COLUMNS_DROPPED_FOR_EXCEL = {
    "Full Term Gene List":
        "up to 125,086 characters per cell (the full human annotation of "
        "'protein binding'), 3.8x Excel's 32,767-character cell limit. It is a "
        "property of the GO annotation (go-basic.obo + goa_human.gaf, both cited "
        "in the Methods), not a result of this study. 'Overlapping Genes', which "
        "is the result, is kept in full.",
}

# The columns the April Supplementary File 1.csv shipped: the published table
# without the two epigenomic columns, which that file already omitted.
TSS_COLUMNS = (
    ["Chromosome", "Start", "End", "Gene_name", "Divergence_scores", "TE_subfamilies",
     "TE_families", "TE_classes", "Average_Divergence_Score"]
    + [f"{c}_number" for c in rl.CLASS_NAMES]
    + ["TE_number"]
    + [f"Divergence_Avg_{c}" for c in rl.CLASS_NAMES]
)

# Old gene-set workbook sheet -> the group label it becomes in long format. The
# April sheet names are kept as the source of truth for what each set is.
GENE_SET_SOURCES = {
    "by_TE_group": ("Supplementary File 3.xlsx", None),
    "by_divergence": ("Supplementary File 5.xlsx", None),
    "by_family": ("Supplementary File 7.xlsx", None),
}

# Excel mangled these two family names in the April File S7 (it forbids `?`).
# Long format restores them.
SHEET_NAME_REPAIRS = {"hAT...": "hAT?", "hAT-Tip100...": "hAT-Tip100?"}

# Filled by sheet_gene_sets: April sheets that were EMPTY and had to be reconstructed.
# Reported in the workbook README so the repair is declared, not silent.
EMPTY_SHEETS_REPAIRED: list[str] = []

FIGURE_COUNT = 14


def sha256(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require(path: str, produced_by: str) -> str:
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"{os.path.relpath(path, rl.REPO_DIR)} is missing — produced by {produced_by}"
        )
    return path


def read_csv(path: str, **kwargs) -> pd.DataFrame:
    return pd.read_csv(path, low_memory=False, **kwargs)


# --------------------------------------------------------------------------- #
# Sheet builders
# --------------------------------------------------------------------------- #


def sheet_tss_intersections() -> pd.DataFrame:
    table = read_csv(require(os.path.join(rl.REPO_DIR, "TEs_on_genes.csv"),
                             "the published pipeline"), index_col=0)
    return table[TSS_COLUMNS]


def sheet_enrichment_classes() -> pd.DataFrame:
    return read_csv(require(os.path.join(OUTPUT_DIR, "TableS_class_enrichment_full.csv"),
                            "10_tables.py"))


def sheet_enrichment_families() -> pd.DataFrame:
    table = read_csv(require(os.path.join(rl.REPO_DIR,
                                          "enrichment_families_with_random.csv"),
                             "the published pipeline"))
    return table.drop(columns=[c for c in table.columns if c.startswith("Unnamed")])


def sheet_enrichment_subfamilies() -> pd.DataFrame:
    # The first column holds the subfamily names and is unnamed, while a second
    # `subfamily_name` column also exists — pandas resolves the clash by renaming one
    # to `subfamily_name.1`. Reading the index explicitly avoids shipping that name.
    table = read_csv(require(os.path.join(rl.REPO_DIR,
                                          "enrichment_subfamilies_with_random.csv"),
                             "the published pipeline"), index_col=0)
    table.index.name = "subfamily_name"
    table = table.reset_index()
    return table.drop(columns=[c for c in table.columns if c.endswith(".1")])


def _read_old_gene_set_workbook(filename: str) -> dict[str, list[str]]:
    """Every sheet of an April gene-set workbook as {set name: [genes]}.

    Read with `header=None`: those sheets carry no header row, so the first gene sits
    in the column name and would otherwise be silently lost from every set.
    """
    path = require(os.path.join(OLD_DIR, filename), "the April 2026 submission")
    workbook = pd.ExcelFile(path)
    sets = {}
    for sheet in workbook.sheet_names:
        frame = workbook.parse(sheet, header=None)
        name = SHEET_NAME_REPAIRS.get(sheet, sheet)
        if frame.empty:
            sets[name] = []
            continue
        # File S7's sheets DO have a header ("Unnamed: 0", "Gene name") and an index
        # column; the others are a bare single column of genes.
        column = frame.columns[-1]
        genes = [str(g) for g in frame[column].dropna().tolist()]
        if genes and genes[0] == "Gene name":
            genes = genes[1:]
        sets[name] = genes
    return sets


def _reconstruct_class_gene_sets() -> tuple[dict[str, list[str]], dict[str, dict]]:
    """The eight class-level gene sets, rebuilt from the published gene table.

    Needed because two sheets of the April `Supplementary File 3.xlsx` — `TE top` and
    `TE bottom` — are **empty**: those gene sets drive GO results the paper reports and
    the track hub ships, but they were never actually included in the submitted
    supplementary material. The construction is imported from `06_go_rerun_fdr005.py`'s
    constants rather than retyped, and it is *verified* against the six sheets that are
    not empty before it is used for the two that are.
    """
    go_rerun = importlib.import_module("06_go_rerun_fdr005")
    df = read_csv(os.path.join(rl.REPO_DIR, "TEs_on_genes.csv"), index_col=0)

    def ranked(column, ascending=False):
        return list(
            df.sort_values(column, ascending=ascending)[["Gene_name"]]
            .drop_duplicates()["Gene_name"]
        )

    def per_gene_max(column):
        """The gene's own count, which is what the sort+drop_duplicates cut ranks on."""
        return df.groupby("Gene_name")[column].max()

    def boundary(column, ascending=False):
        """The count at the cut position — the value every arbitrary tie sits at."""
        values = sorted(per_gene_max(column), reverse=not ascending)
        counts = per_gene_max(column)
        cut = values[go_rerun.CLASS_TOP_N - 1]
        record = counts.to_dict()
        record["__boundary__"] = cut
        return record

    sets, boundaries = {}, {}
    # Frozen cell 20's sheet labels: the class sets are plural, SVA/Helitron are not.
    for class_name, sheet in [("LINE", "LINEs"), ("SINE", "SINEs"), ("LTR", "LTRs"),
                              ("DNA", "DNA")]:
        column = f"{class_name}_number"
        sets[sheet] = ranked(column)[:go_rerun.CLASS_TOP_N]
        boundaries[sheet] = boundary(column)
    # SVA and Helitron are every gene carrying such an element, not a ranked cut, so they
    # have no boundary and must agree with the April sheets exactly.
    for sheet, column in [("SVA", "Retroposon_number"), ("Helitron", "RC_number")]:
        sets[sheet] = list(df[df[column] != 0]["Gene_name"].unique())
        boundaries[sheet] = {"__boundary__": None}
    sets["TE top"] = ranked("TE_number", ascending=False)[:go_rerun.CLASS_TOP_N]
    boundaries["TE top"] = boundary("TE_number", ascending=False)
    sets["TE bottom"] = ranked("TE_number", ascending=True)[:go_rerun.CLASS_TOP_N]
    boundaries["TE bottom"] = boundary("TE_number", ascending=True)
    return sets, boundaries


def _class_gene_sets_with_empty_sheets_repaired() -> tuple[dict[str, list[str]], list[str]]:
    """The April class-level sets, with the two empty sheets filled from a verified rebuild."""
    published = _read_old_gene_set_workbook("Supplementary File 3.xlsx")
    empty = [name for name, genes in published.items() if not genes]
    if not empty:
        return published, []

    rebuilt, boundary_counts = _reconstruct_class_gene_sets()
    # The rebuild is only trustworthy for the empty sheets if it reproduces the others.
    # It does NOT reproduce them exactly, and the reason is not a construction
    # difference: a top-1,436 cut lands inside a large block of genes that all carry the
    # SAME element count (for LINEs, 1,257 genes are strictly above the boundary count of
    # 9 and 1,033 genes tie at it, competing for the remaining 179 places), so which tied
    # genes get in is decided by row order and is arbitrary. The check therefore demands
    # that every disagreeing gene sit exactly ON the boundary count — if any gene above it
    # differed, the construction really would have changed.
    for name, genes in published.items():
        if not genes:
            continue
        if name not in rebuilt:
            raise KeyError(f"no reconstruction for the April sheet {name!r}")
        differing = set(genes) ^ set(rebuilt[name])
        off_boundary = {
            gene for gene in differing
            if boundary_counts[name].get(gene) != boundary_counts[name]["__boundary__"]
        }
        if off_boundary:
            raise ValueError(
                f"the reconstruction differs from the April {name!r} sheet for "
                f"{len(off_boundary)} gene(s) that are NOT at the boundary count "
                f"{boundary_counts[name]['__boundary__']} — the construction has changed, "
                f"so it cannot be trusted to fill {empty}. Examples: "
                f"{sorted(off_boundary)[:5]}"
            )
        cut = boundary_counts[name]["__boundary__"]
        print(f"    {name:10s} {len(differing):>4} of {len(genes):,} genes differ"
              + (f", all tied at the boundary count {cut} (arbitrary either way)"
                 if differing else " (not a ranked cut)" if cut is None else ""))
    print(f"  reconstruction agrees with all "
          f"{len([g for g in published.values() if g])} non-empty April sheets except at "
          f"the tie boundary; using it for the empty ones: {empty}")
    for name in empty:
        published[name] = rebuilt[name]
    return published, empty


def sheet_gene_sets(filename: str, group_label: str) -> pd.DataFrame:
    """One long-format table: group, gene, rank — replacing 8, 10 or 44 sheets."""
    if filename == "Supplementary File 3.xlsx":
        sets, repaired = _class_gene_sets_with_empty_sheets_repaired()
        EMPTY_SHEETS_REPAIRED.extend(repaired)
    else:
        sets = _read_old_gene_set_workbook(filename)
        for name, genes in sets.items():
            if not genes:
                raise ValueError(f"{filename}: sheet {name!r} is empty and there is no "
                                 f"verified reconstruction for it")
    rows = []
    for set_name, genes in sets.items():
        for rank, gene in enumerate(genes, 1):
            rows.append({group_label: set_name, "gene": gene, "rank": rank})
    return pd.DataFrame(rows)


def sheet_go(level: str, classification_file: str, sep: str = ",") -> pd.DataFrame:
    """A GO level at FDR 0.05, with the manual functional-group classification."""
    table = read_csv(require(os.path.join(OUTPUT_DIR, f"GO_{level}_fdr005.csv"),
                             "06_go_rerun_fdr005.py"))
    classification = _load_classification(classification_file, sep)
    merged = table.merge(classification, how="left", on="Term Name")
    return merged.drop(columns=[c for c in merged.columns if c.startswith("Unnamed")])


def _load_classification(filename: str, sep: str) -> pd.DataFrame:
    """The manual GO term -> functional group map, as nb06 cell 16 loads it."""
    frame = read_csv(require(os.path.join(rl.REPO_DIR, filename),
                             "manual curation"), sep=sep)
    column = ("Functional group Gemini" if "Functional group Gemini" in frame.columns
              else "Functional Group")
    name_column = "GO Term Name" if "GO Term Name" in frame.columns else "Term Name"
    frame = frame[[name_column, column]].dropna()
    frame.columns = ["Term Name", "Functional Group"]
    return frame.drop_duplicates("Term Name")


def grid_sheet(name: str, produced_by: str, optional: bool = False):
    """A GO-grid table from nb07, which must have been executed first."""
    path = os.path.join(OUTPUT_DIR, name)
    if not os.path.exists(path) and optional:
        return None
    return read_csv(require(path, produced_by))


# --------------------------------------------------------------------------- #
# Workbook definitions
# --------------------------------------------------------------------------- #


def workbook_definitions(allow_missing: bool) -> list[dict]:
    """The five workbooks, each as {filename, subject, sheets: [(name, frame, note)]}."""
    optional = allow_missing

    grid_index = grid_sheet("GO_grid_index.csv", "nb07 (or 07b's INDEX.csv)",
                            optional=True)
    if grid_index is None and os.path.exists(os.path.join(GRID_DIR, "INDEX.csv")):
        grid_index = read_csv(os.path.join(GRID_DIR, "INDEX.csv"))

    sensitivity_sheets = [
        ("window_classes", read_csv(require(
            os.path.join(OUTPUT_DIR, "window_sensitivity_classes.csv"),
            "05b_window_sensitivity.py")),
         "class-level enrichment at 5, 10 and 20 kb"),
        ("window_families", read_csv(os.path.join(
            OUTPUT_DIR, "window_sensitivity_families.csv")),
         "family-level enrichment at 5, 10 and 20 kb"),
        ("window_concordance", read_csv(os.path.join(
            OUTPUT_DIR, "window_sensitivity_concordance.csv")),
         "Spearman rho per window pair with a label-shuffling permutation test"),
        ("window_flips", read_csv(os.path.join(
            OUTPUT_DIR, "window_sensitivity_flips.csv")),
         "every group whose significance call changes between windows"),
        ("percentile_summary", read_csv(os.path.join(
            OUTPUT_DIR, "percentile_sensitivity_summary.csv")),
         "GO term counts and stability per group, top/bottom 5 % vs 10 %"),
        ("percentile_terms", read_csv(os.path.join(
            OUTPUT_DIR, "percentile_sensitivity_terms.csv")),
         "every GO term gained or lost at 10 %, named"),
        ("geneset_stability", read_csv(os.path.join(
            OUTPUT_DIR, "robustness_geneset_stability.csv")),
         "overlap of the top-5 % gene sets between windows, with a "
         "hypergeometric p-value"),
        ("rank_stability", read_csv(os.path.join(
            OUTPUT_DIR, "robustness_rank_stability.csv")),
         "Kendall tau of the per-gene ranking between windows, bootstrap CI"),
    ]
    for name, filename, note in [
        ("GO_grid_index", "GO_grid_index.csv",
         "one row per cell of the 3 windows x 2 percentiles GO grid"),
        ("GO_grid_preservation", "go_grid_preservation.csv",
         "Jaccard and fraction of published terms preserved, per cell and level"),
        ("GO_grid_terms", "go_grid_terms.csv",
         "every GO term gained or lost in each cell relative to 10 kb / 5 %"),
        ("GO_grid_concordance", "go_grid_concordance.csv",
         "Spearman of per-group term counts against the published cell, with a "
         "label-shuffling permutation test"),
        ("headline_by_condition", "go_grid_headline_by_condition.csv",
         "each headline claim under all six conditions"),
    ]:
        frame = (grid_index if name == "GO_grid_index"
                 else grid_sheet(filename, "nb07_go_grid_robustness.ipynb",
                                 optional=optional))
        if frame is None:
            print(f"  ! {name}: {filename} not on disk — sheet omitted "
                  f"(run nb07 first)")
            continue
        sensitivity_sheets.append((name, frame, note))

    return [
        {
            "filename": "File_S1_TE_TSS_map_and_enrichment.xlsx",
            "subject": "The TE-TSS map and the enrichment statistics",
            "sheets": [
                ("TSS_TE_intersections", sheet_tss_intersections(),
                 "every TSS window with the transposable elements within 10 kb of it: "
                 "counts and mean divergence per class, and the comma-joined element "
                 "lists. Counts are per TSS window, not per gene, so a gene with "
                 "several annotated TSS contributes several rows"),
                ("enrichment_classes", sheet_enrichment_classes(),
                 "the six TE classes: observed and random odds ratios, raw AND "
                 "adjusted Fisher p, and the empirical permutation p (N = 500)"),
                ("enrichment_families", sheet_enrichment_families(),
                 "the 44 TE families, same columns"),
                ("enrichment_subfamilies", sheet_enrichment_subfamilies(),
                 "all 1,143 TE subfamilies, same columns — this is the table the "
                 "File S2 caption described"),
            ],
        },
        {
            "filename": "File_S2_gene_sets.xlsx",
            "subject": "The foreground gene sets used for GO enrichment",
            "sheets": [
                ("by_TE_group", sheet_gene_sets("Supplementary File 3.xlsx",
                                                "TE_group"),
                 "top 5 % of genes by element count for each TE group, plus the "
                 "TE-top and TE-bottom sets (was File S3, 8 sheets)"),
                ("by_divergence", sheet_gene_sets("Supplementary File 5.xlsx",
                                                  "divergence_group"),
                 "highest- and lowest-divergence 5 % of genes per class "
                 "(was File S5, 10 sheets)"),
                ("by_family", sheet_gene_sets("Supplementary File 7.xlsx",
                                              "family_name"),
                 "top 5 % of genes by element count for each of the 44 families "
                 "(was File S7, 44 sheets; the hAT? and hAT-Tip100? names are "
                 "restored here — Excel had mangled them into sheet names)"),
            ],
        },
        {
            "filename": "File_S3_gene_ontology.xlsx",
            "subject": "Gene Ontology enrichment at FDR < 0.05",
            "sheets": [
                ("GO_TE_groups", sheet_go("classes_count",
                                          "GO_top_5_perc_genes_by_class_gemini_manual.csv",
                                          sep=";"),
                 "GO terms enriched in each TE group's gene set, with the manual "
                 "functional-group classification (was File S4, at FDR 0.1)"),
                ("GO_by_divergence", sheet_go(
                    "classes_divergence",
                    "GO_terms_divergence_classification_Gemini.csv"),
                 "GO terms for the highest- and lowest-divergence gene sets "
                 "(was File S6, at FDR 0.1)"),
                ("GO_by_family", sheet_go(
                    "families",
                    "GO_terms_families_classification_Gemini - families.csv"),
                 "GO terms per TE family (was File S8, at FDR 0.1)"),
            ],
        },
        {
            "filename": "File_S4_IFNA_domain_and_prior_work.xlsx",
            "subject": "The interferon-alpha domain, the assembly bound, and the "
                       "overlap with prior work",
            "sheets": [
                ("IFNA_elements", read_csv(require(
                    os.path.join(OUTPUT_DIR, "ifna_window_elements.csv"),
                    "02_ifna_domain_test.py")),
                 "all 175 transposable elements in the 220 kb interferon-alpha "
                 "domain (chr9:21,150,692-21,370,055)"),
                ("IFNA_tests", read_csv(os.path.join(
                    OUTPUT_DIR, "ifna_test_results.csv")),
                 "the four domain tests: observed value, null mean and SD, and the "
                 "empirical p-value"),
                ("IFNA_subfamily_composition", read_csv(os.path.join(
                    OUTPUT_DIR, "ifna_subfamily_composition.csv")),
                 "the L1 subfamilies represented in the domain and their divergence"),
                ("assembly_bound", _assembly_bound(),
                 "newly resolved T2T sequence and the share of elements, windows "
                 "and genes it contributes, by TE class, family and chromosome"),
                ("prior_work_overlap_matrix", read_csv(require(
                    os.path.join(OUTPUT_DIR, "lu2020_overlap_matrix.csv"),
                    "04a_lu2020_geneset_overlap.py")),
                 "Lu et al. 2020 category x TE group: overlap size, Fisher p and "
                 "Jaccard"),
                ("prior_work_categories", read_csv(os.path.join(
                    OUTPUT_DIR, "lu2020_categories_mapped.csv")),
                 "Lu et al.'s mouse categories mapped to human genes through MGI "
                 "orthology"),
                ("prior_work_shared_genes", read_csv(os.path.join(
                    OUTPUT_DIR, "lu2020_overlap_genes.csv")),
                 "the individual genes shared between each category and each group"),
            ],
        },
        {
            "filename": "File_S5_sensitivity_and_robustness.xlsx",
            "subject": "Window-size, percentile and GO-grid sensitivity analyses",
            "sheets": sensitivity_sheets,
        },
    ]


def _assembly_bound() -> pd.DataFrame:
    """The newly resolved sequence tables, stacked with a `grouping` column."""
    frames = []
    for grouping, filename in [("class", "assembly_bound_by_class.csv"),
                               ("family", "assembly_bound_by_family.csv"),
                               ("chromosome", "assembly_bound_by_chromosome.csv")]:
        frame = read_csv(require(os.path.join(OUTPUT_DIR, filename),
                                 "04b_newly_resolved_regions.py"))
        frames.append(frame.assign(grouping=grouping))
    stacked = pd.concat(frames, ignore_index=True)
    return stacked[["grouping"] + [c for c in stacked.columns if c != "grouping"]]


# --------------------------------------------------------------------------- #
# Writing, with the Excel limits enforced
# --------------------------------------------------------------------------- #


def prepare_sheet(workbook_name: str, sheet_name: str,
                  frame: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    """Validate a sheet against Excel's limits, dropping only declared columns."""
    problems = []
    if len(sheet_name) > EXCEL_SHEET_NAME_LIMIT:
        problems.append(f"sheet name is {len(sheet_name)} > "
                        f"{EXCEL_SHEET_NAME_LIMIT} characters")
    if EXCEL_FORBIDDEN_IN_SHEET_NAME.search(sheet_name):
        problems.append("sheet name contains a character Excel forbids")
    if len(frame) + 1 > EXCEL_ROW_LIMIT:
        problems.append(f"{len(frame):,} rows exceeds Excel's {EXCEL_ROW_LIMIT:,}")

    dropped = [c for c in frame.columns if c in COLUMNS_DROPPED_FOR_EXCEL]
    frame = frame.drop(columns=dropped)

    for column in frame.columns:
        if not (pd.api.types.is_string_dtype(frame[column])
                or frame[column].dtype == object):
            continue
        longest = int(frame[column].astype(str).str.len().max() or 0)
        if longest > EXCEL_CELL_LIMIT:
            problems.append(
                f"column {column!r} has a {longest:,}-character cell, over Excel's "
                f"{EXCEL_CELL_LIMIT:,} limit, and is not on the declared drop list"
            )
    if problems:
        raise ValueError(f"{workbook_name} / {sheet_name}: " + "; ".join(problems))
    return frame, dropped


def readme_frame(definition: dict, sheets: list[tuple],
                 dropped: dict[str, list[str]]) -> pd.DataFrame:
    """The README sheet: what each sheet is, its size, and which script wrote it."""
    rows = [{
        "sheet": "README",
        "rows": "",
        "columns": "",
        "description": (
            f"{definition['subject']}. Built by revision_G3/"
            f"{os.path.basename(__file__)} for manuscript G3-2026-406828. "
            f"Every Gene Ontology result in this package is filtered at "
            f"FDR < {rl.FDR_THRESHOLD} (Benjamini-Hochberg); the permutation "
            f"background is N = {rl.N_PERMUTATIONS}, so the empirical p floor is "
            f"{rl.EMPIRICAL_P_FLOOR:.4f}."
        ),
    }]
    for name, frame, note in sheets:
        description = note
        if dropped.get(name):
            reasons = "; ".join(
                f"{c} removed: {COLUMNS_DROPPED_FOR_EXCEL[c]}" for c in dropped[name]
            )
            description = f"{note}. {reasons}"
        rows.append({
            "sheet": name,
            "rows": len(frame),
            "columns": len(frame.columns),
            "description": description,
        })
    return pd.DataFrame(rows)


def build() -> int:
    os.makedirs(SUPP_DIR, exist_ok=True)
    os.makedirs(os.path.join(SUPP_DIR, "figures"), exist_ok=True)

    inventory = {}
    for definition in workbook_definitions(allow_missing=True):
        print(f"\n=== {definition['filename']}")
        prepared, dropped = [], {}
        for name, frame, note in definition["sheets"]:
            frame, dropped_columns = prepare_sheet(definition["filename"], name, frame)
            if dropped_columns:
                dropped[name] = dropped_columns
            prepared.append((name, frame, note))
            print(f"  {name:28s} {len(frame):>8,} rows x {len(frame.columns):>3} cols"
                  + (f"  (dropped {', '.join(dropped_columns)})"
                     if dropped_columns else ""))

        path = os.path.join(SUPP_DIR, definition["filename"])
        with pd.ExcelWriter(path, engine="openpyxl") as writer:
            readme_frame(definition, prepared, dropped).to_excel(
                writer, sheet_name="README", index=False)
            for name, frame, _note in prepared:
                frame.to_excel(writer, sheet_name=name, index=False)
        size = os.path.getsize(path)
        print(f"  wrote {definition['filename']} ({size / 1e6:.1f} MB, "
              f"{len(prepared) + 1} sheets)")
        inventory[definition["filename"]] = {
            "subject": definition["subject"],
            "size_bytes": size,
            "sheets": {name: {"rows": len(frame), "columns": list(frame.columns),
                              "description": note}
                       for name, frame, note in prepared},
            "columns_dropped_for_excel": dropped,
        }

    write_readme(inventory)
    write_placeholders()
    with open(os.path.join(SUPP_DIR, "INVENTORY.json"), "w") as handle:
        json.dump(inventory, handle, indent=2)
    # Last, so it covers every file this run wrote — including INVENTORY.json itself.
    write_checksums()
    print(f"\nwrote INVENTORY.json, README.md, CHECKSUMS.sha256 and "
          f"figures/PLACEHOLDERS.md in {os.path.relpath(SUPP_DIR, rl.REPO_DIR)}/")
    return 0


def write_readme(inventory: dict) -> None:
    lines = [
        "# Supplementary material — G3-2026-406828",
        "",
        "*Telomere-to-telomere co-mapping of transposable elements and human genes "
        "identifies a cluster of young L1 elements in the interferon-alpha domain*",
        "",
        f"Assembled by `revision_G3/{os.path.basename(__file__)}`. Five thematic "
        "workbooks, one sheet per table, each workbook opening with a README sheet "
        "that describes its sheets.",
        "",
        "## What changed relative to the April 2026 submission",
        "",
        "1. **Every Gene Ontology result is now at FDR < 0.05** (Benjamini-Hochberg), "
        "with no \"suggestive\" band. The April supplementary files were at FDR 0.1. "
        "This is the change reviewer minor comment 2 asked for.",
        "2. **Fourteen candidate files became five workbooks**, so every `File Sn` "
        "citation in the manuscript is renumbered and now names a sheet as well as a "
        "file. The mapping from the old numbering is below.",
        "3. **The gene-set tables are long format** rather than one sheet per group. "
        "This also restores the `hAT?` and `hAT-Tip100?` family names, which Excel "
        "had mangled into sheet names in the April File S7.",
        "4. **Workbook S1 carries both the family and the subfamily enrichment "
        "tables.** The April File S2 caption described subfamilies but the file "
        "contained the 44 families; both are now present with the full column set.",
        "",
        "## Old numbering -> new",
        "",
        "| April 2026 | Now |",
        "|---|---|",
        "| File S1 (TSS/TE coordinates) | File S1, sheet `TSS_TE_intersections` |",
        "| File S2 (enrichment) | File S1, sheets `enrichment_families` / "
        "`enrichment_subfamilies` |",
        "| File S3 (gene sets by TE group) | File S2, sheet `by_TE_group` |",
        "| File S4 (GO by TE group) | File S3, sheet `GO_TE_groups` |",
        "| File S5 (gene sets by divergence) | File S2, sheet `by_divergence` |",
        "| File S6 (GO by divergence) | File S3, sheet `GO_by_divergence` |",
        "| File S7 (gene sets by family) | File S2, sheet `by_family` |",
        "| File S8 (GO by family) | File S3, sheet `GO_by_family` |",
        "| — (Lu et al. overlap) | File S4, sheet `prior_work_overlap_matrix` |",
        "| — (\"the accompanying tables\") | File S5 |",
        "",
        "## The workbooks",
        "",
    ]
    for filename, record in inventory.items():
        lines += [f"### `{filename}`", "", record["subject"] + ".", "",
                  "| Sheet | Rows | Columns | Contents |", "|---|---|---|---|"]
        for name, sheet in record["sheets"].items():
            lines.append(f"| `{name}` | {sheet['rows']:,} | "
                         f"{len(sheet['columns'])} | {sheet['description']} |")
        lines.append("")
    lines += [
        "## Conventions that apply throughout",
        "",
        f"* **Element counts are per TSS window, not per gene.** A gene with several "
        f"annotated TSS contributes several windows, so an element within 10 kb of "
        f"two TSS of the same gene is counted twice. This is a property of the "
        f"published design and is stated in the manuscript's Limitations.",
        f"* **GO FDR is {rl.FDR_THRESHOLD}** everywhere, Benjamini-Hochberg corrected.",
        f"* **The permutation background is N = {rl.N_PERMUTATIONS}**, so the "
        f"empirical p-value floor is 2/{rl.N_PERMUTATIONS + 1} = "
        f"{rl.EMPIRICAL_P_FLOOR:.4f}.",
        "* **Raw and adjusted p-values are both reported** wherever a correction was "
        "applied, per G3's statistics policy.",
        "",
        "## One column is not in the workbooks",
        "",
    ]
    for column, reason in COLUMNS_DROPPED_FOR_EXCEL.items():
        lines.append(f"`{column}` — {reason}")
    if EMPTY_SHEETS_REPAIRED:
        lines += [
            "## Two gene sets were empty in the April submission",
            "",
            f"The April `Supplementary File 3.xlsx` shipped "
            f"**{len(EMPTY_SHEETS_REPAIRED)} empty sheets**: "
            + ", ".join(f"`{name}`" for name in EMPTY_SHEETS_REPAIRED)
            + ". Those gene sets drive GO results the paper reports and gene tracks the "
            "UCSC hub ships, so the sets themselves were never in doubt — only their "
            "inclusion in the supplementary file. They are reconstructed here by the "
            "same construction `revision_G3/06_go_rerun_fdr005.py` uses (top and bottom "
            "1,436 genes by total element count), which was checked against the six "
            "non-empty April sheets first.",
            "",
            "That check is not an exact match, and the reason matters: a top-1,436 cut "
            "lands inside a large block of genes carrying the *same* element count — for "
            "the LINEs set, 1,257 genes are strictly above the boundary count of 9 while "
            "1,033 tie at it, competing for the remaining 179 places. Which tied genes "
            "make the cut is decided by row order and is arbitrary. The verification "
            "therefore requires every disagreeing gene to sit exactly **on** the boundary "
            "count, and it does; no gene above the boundary differs, in any set. The same "
            "arbitrariness applies to the two reconstructed sets, and to the published "
            "ones.",
            "",
        ]
    lines += ["", "## Figures", "",
              "The figures are exported from Figma by hand and are listed in "
              "`figures/PLACEHOLDERS.md`.", ""]
    with open(os.path.join(SUPP_DIR, "README.md"), "w") as handle:
        handle.write("\n".join(lines))


def write_placeholders() -> None:
    lines = [
        "# Expected supplementary figure exports",
        "",
        "These are exported from Figma by hand (revision Phase 4b) and are not "
        "produced by any script here — `revision_G3/svg/PLACEMENT.md` maps every SVG "
        "panel to its Figma frame. Drop the PDFs into this directory under exactly "
        "these names.",
        "",
    ]
    titles = {
        1: "Divergence distributions per TE class, all TEs vs TEs near genes",
        2: "Per-family divergence distributions, 44 families",
        3: "Per-family element length distributions, 44 families",
        4: "Observed vs random odds ratios per subfamily",
        5: "Subfamily enrichment against element abundance",
        6: "GO term counts against enrichment statistics (promoted to main Figure 8A)",
        7: "Functional-group clustermap, subfamily level",
        8: "Combined class and family functional-group clustermap, and the "
           "unfiltered class-group-family Sankey (panel C)",
        9: "Full class-level GO network by element count",
        10: "Full class-level GO network by divergence",
        11: "Full family-level GO network",
        12: "GO robustness across the 3 windows x 2 percentiles grid",
        13: "Window and percentile concordance",
        14: "Permutation background convergence at N = 500",
    }
    for number in range(1, FIGURE_COUNT + 1):
        lines.append(f"- [ ] `Figure_S{number}.pdf` — {titles[number]}")
    lines += ["", f"{FIGURE_COUNT} figures, S1-S{FIGURE_COUNT}, no gaps.", ""]
    with open(os.path.join(SUPP_DIR, "figures", "PLACEHOLDERS.md"), "w") as handle:
        handle.write("\n".join(lines))


def write_checksums() -> None:
    lines = []
    for name in sorted(os.listdir(SUPP_DIR)):
        path = os.path.join(SUPP_DIR, name)
        if os.path.isfile(path) and name != "CHECKSUMS.sha256":
            lines.append(f"{sha256(path)}  {name}")
    with open(os.path.join(SUPP_DIR, "CHECKSUMS.sha256"), "w") as handle:
        handle.write("\n".join(lines) + "\n")


def verify() -> int:
    """Check the package on disk: five workbooks, sheet limits, GO at 0.05."""
    print("=" * 78)
    print("08 verification — the supplementary package on disk")
    print("=" * 78)
    failures = []

    workbooks = sorted(f for f in os.listdir(SUPP_DIR)
                       if f.startswith("File_S") and f.endswith(".xlsx")) \
        if os.path.isdir(SUPP_DIR) else []
    print(f"\n{len(workbooks)} workbooks: {workbooks}")
    if len(workbooks) != 5:
        failures.append(f"expected 5 workbooks, found {len(workbooks)}")

    for filename in workbooks:
        path = os.path.join(SUPP_DIR, filename)
        workbook = pd.ExcelFile(path)
        print(f"\n{filename} ({os.path.getsize(path) / 1e6:.1f} MB, "
              f"{len(workbook.sheet_names)} sheets)")
        if "README" not in workbook.sheet_names:
            failures.append(f"{filename}: no README sheet")
        for sheet in workbook.sheet_names:
            if len(sheet) > EXCEL_SHEET_NAME_LIMIT:
                failures.append(f"{filename}/{sheet}: sheet name too long")
            frame = workbook.parse(sheet)
            note = ""
            if "FDR" in frame.columns and len(frame):
                worst = float(frame["FDR"].max())
                ok = worst < rl.FDR_THRESHOLD
                note = f"  max FDR {worst:.4g} {'OK' if ok else 'ABOVE THRESHOLD'}"
                if not ok:
                    failures.append(f"{filename}/{sheet}: FDR up to {worst:.4g}")
            print(f"  {sheet:28s} {len(frame):>8,} rows x "
                  f"{len(frame.columns):>3} cols{note}")

    checksums = os.path.join(SUPP_DIR, "CHECKSUMS.sha256")
    if os.path.exists(checksums):
        bad = []
        for line in open(checksums):
            digest, name = line.strip().split("  ", 1)
            path = os.path.join(SUPP_DIR, name)
            if not os.path.exists(path) or sha256(path) != digest:
                bad.append(name)
        print(f"\nchecksums: {'all match' if not bad else f'MISMATCH {bad}'}")
        if bad:
            failures.append(f"checksum mismatch: {bad}")
    else:
        failures.append("CHECKSUMS.sha256 missing")

    print("-" * 78)
    if failures:
        for failure in failures:
            print(f"FAIL {failure}")
        return 1
    print("VERIFIED: five workbooks, every sheet within Excel's limits, every GO "
          f"sheet at FDR < {rl.FDR_THRESHOLD}, checksums match.")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--verify", action="store_true",
                        help="check the package on disk instead of rebuilding it")
    args = parser.parse_args()
    return verify() if args.verify else build()


if __name__ == "__main__":
    sys.exit(main())
