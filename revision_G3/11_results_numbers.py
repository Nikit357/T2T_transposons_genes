#!/usr/bin/env python
"""Collect every number the revised manuscript text quotes, in one place.

Phase 5 of the G3 revision rewrites Methods, Results, Discussion and Data availability.
Each of those edits quotes a number that came out of Phase 2-4. This script re-derives all
of them from the persisted CSV/JSON outputs so that (a) nothing in the manuscript is typed
from memory, and (b) if any upstream output is regenerated the discrepancy shows up here
rather than in proof.

It computes nothing new about the biology: every value either comes straight from an output
file or is a count/aggregation over one. The two exceptions are the functional-group
contingency counts and the Figure 7 correlations, which nb06 renders into SVG without
writing the underlying numbers to disk -- they are recomputed here with the same code path.

Outputs
    output/results_numbers.json  machine-readable, one key per manuscript claim
    output/results_numbers.txt   the same, formatted for reading next to the docx

Usage
    python 11_results_numbers.py
"""

import json
import os

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
OUT = os.path.join(HERE, "output")

FDR = 0.05

# The two levels the manuscript quotes GO counts for, plus their FDR 0.1 references so that
# every "lost at 0.05" statement in the revised text can be justified.
LEVELS = {
    "classes_count": ("GO_classes_count_fdr005.csv",
                      "GO_classes_count_fdr01_reference.csv", "class_name"),
    "classes_divergence": ("GO_classes_divergence_fdr005.csv",
                           "GO_classes_divergence_fdr01_reference.csv", "class_name"),
    "families": ("GO_families_fdr005.csv",
                 "GO_families_fdr01_reference.csv", "family_name"),
    "subfamilies": ("GO_subfamilies_fdr005.csv",
                    "GO_subfamilies_fdr01_reference.csv", "subfamily_name"),
}


def read_output(name, **kwargs):
    return pd.read_csv(os.path.join(OUT, name), **kwargs)


def load_classification(path, sep=","):
    """Manual GO term -> functional group map, exactly as nb06 loads it."""
    frame = pd.read_csv(os.path.join(REPO, path), sep=sep, low_memory=False)
    column = ("Functional group Gemini" if "Functional group Gemini" in frame.columns
              else "Functional Group")
    frame = frame[["GO Term Name" if "GO Term Name" in frame.columns else "Term Name",
                   column]].dropna()
    frame.columns = ["Term Name", "Functional Group"]
    return frame.drop_duplicates("Term Name")


def functional_group_counts(go_table, classification, group_column):
    """GO terms per (group, functional group) plus the BH-corrected Fisher stars.

    Same contingency construction as frozen cell 70 and nb06 cell 17: each cell is that
    group's terms in that functional category against every other cell, one-sided.
    """
    merged = go_table.merge(classification, how="left", on="Term Name")
    unclassified = int(merged["Functional Group"].isna().sum())
    counts = (merged.dropna(subset=["Functional Group"])[[group_column, "Functional Group"]]
              .value_counts().reset_index(name="count"))
    pivot = counts.pivot(index=group_column, columns="Functional Group",
                         values="count").fillna(0)

    total = pivot.values.sum()
    rows, cols = pivot.sum(axis=1), pivot.sum(axis=0)
    p_raw = np.ones(pivot.shape)
    for i, r in enumerate(pivot.index):
        for j, c in enumerate(pivot.columns):
            a = pivot.loc[r, c]
            b, d = rows[r] - a, cols[c] - a
            p_raw[i, j] = fisher_exact([[a, b], [d, total - (a + b + d)]],
                                       alternative="greater")[1]
    p_adj = pd.DataFrame(
        multipletests(p_raw.flatten(), method="fdr_bh")[1].reshape(p_raw.shape),
        index=pivot.index, columns=pivot.columns)
    return pivot, p_adj, unclassified


def significant_cells(pivot, p_adj, threshold=FDR):
    """Every (group, functional group) pair that is significantly overrepresented."""
    out = []
    for r in pivot.index:
        for c in pivot.columns:
            if p_adj.loc[r, c] < threshold and pivot.loc[r, c] > 0:
                out.append({"group": r, "functional_group": c,
                            "n_terms": int(pivot.loc[r, c]),
                            "fdr": float(p_adj.loc[r, c])})
    return sorted(out, key=lambda d: d["fdr"])


def terms_of(table, group_column, group, extra=()):
    sub = table[table[group_column] == group]
    keep = ["Term ID", "Term Name", "FDR", "Overlap Count", *extra]
    return sub[keep].sort_values("FDR").to_dict("records")


def lost_terms(table005, table01, group_column, group):
    kept = set(table005.loc[table005[group_column] == group, "Term Name"])
    sub = table01[(table01[group_column] == group) & (~table01["Term Name"].isin(kept))]
    return [{"Term Name": r["Term Name"], "FDR": float(r["FDR"])}
            for _, r in sub.sort_values("FDR").iterrows()]


results = {}

# ---------------------------------------------------------------------------------------
# 1. GO term totals at both thresholds, every level
# ---------------------------------------------------------------------------------------
tables = {}
totals = {}
for level, (f005, f01, group_column) in LEVELS.items():
    t005, t01 = read_output(f005), read_output(f01)
    tables[level] = (t005, t01, group_column)
    if level == "classes_divergence":
        for frame in (t005, t01):
            frame["group_label"] = (frame["class_name"].astype(str) + " / "
                                    + frame["divergence_group"].astype(str))
        group_column = "group_label"
        tables[level] = (t005, t01, group_column)
    totals[level] = {
        "n_terms_fdr01": int(len(t01)),
        "n_terms_fdr005": int(len(t005)),
        "n_groups_fdr01": int(t01[group_column].nunique()),
        "n_groups_fdr005": int(t005[group_column].nunique()),
        "groups_fdr005": sorted(t005[group_column].unique().tolist()),
        "groups_lost_at_005": sorted(set(t01[group_column]) - set(t005[group_column])),
        "terms_per_group_fdr005": {
            k: int(v) for k, v in t005[group_column].value_counts().items()},
    }
results["go_totals"] = totals

# ---------------------------------------------------------------------------------------
# 2. Class level by count -- Figures 4A / 4B, manuscript lines 93 and 97
# ---------------------------------------------------------------------------------------
t005, t01, gc = tables["classes_count"]
class_classification = (
    pd.read_csv(os.path.join(REPO, "GO_top_5_perc_genes_by_class_gemini_manual.csv"),
                sep=";", low_memory=False)[["Term Name", "Functional group Gemini"]]
    .dropna().drop_duplicates("Term Name")
    .rename(columns={"Functional group Gemini": "Functional Group"}))
pivot_4b, padj_4b, unclassified_4b = functional_group_counts(
    t005, class_classification, "class_name")
results["classes_count"] = {
    "terms_per_class": totals["classes_count"]["terms_per_group_fdr005"],
    "classes_with_no_terms": sorted(
        {"LINE", "LTR", "SINE", "DNA", "SVA", "RC", "TE_top", "TE_bottom"}
        - set(t005["class_name"])),
    "LINE_terms": terms_of(t005, "class_name", "LINE"),
    "LINE_terms_lost_at_005": lost_terms(t005, t01, "class_name", "LINE"),
    "SVA_terms": terms_of(t005, "class_name", "SVA"),
    "LTR_terms": terms_of(t005, "class_name", "LTR"),
    "TE_bottom_top3_lowest_fdr": terms_of(t005, "class_name", "TE_bottom")[:3],
    "TE_top_terms": terms_of(t005, "class_name", "TE_top"),
    "n_functional_groups": int(pivot_4b.shape[1]),
    "unclassified_terms": unclassified_4b,
    "significant_functional_cells": significant_cells(pivot_4b, padj_4b),
    "functional_group_counts": {
        r: {c: int(pivot_4b.loc[r, c]) for c in pivot_4b.columns if pivot_4b.loc[r, c] > 0}
        for r in pivot_4b.index},
}

# ---------------------------------------------------------------------------------------
# 3. Divergence level -- Figures 5A / 5B, manuscript lines 103, 105, 107
# ---------------------------------------------------------------------------------------
t005, t01, gc = tables["classes_divergence"]
divergence_classification = load_classification(
    "GO_terms_divergence_classification_Gemini.csv")
t005_grouped = t005.copy()
t005_grouped["class_name"] = (t005_grouped["class_name"].astype(str) + "_"
                              + t005_grouped["divergence_group"].astype(str))
pivot_5b, padj_5b, unclassified_5b = functional_group_counts(
    t005_grouped, divergence_classification, "class_name")
results["classes_divergence"] = {
    "terms_per_group": totals["classes_divergence"]["terms_per_group_fdr005"],
    "groups_lost_at_005": totals["classes_divergence"]["groups_lost_at_005"],
    "LINE_lowest_terms": terms_of(t005, gc, "LINE / lowest"),
    "LINE_lowest_lost": lost_terms(t005, t01, gc, "LINE / lowest"),
    "LINE_highest_terms": terms_of(t005, gc, "LINE / highest"),
    "LTR_lowest_terms": terms_of(t005, gc, "LTR / lowest"),
    "LTR_highest_terms": terms_of(t005, gc, "LTR / highest"),
    "SINE_lowest_terms": terms_of(t005, gc, "SINE / lowest"),
    "TE_all_lowest_terms": terms_of(t005, gc, "TE_all / lowest"),
    "top3_lowest_fdr_overall": (
        t005.sort_values("FDR")[[gc, "Term Name", "FDR"]].head(3).to_dict("records")),
    "n_functional_groups": int(pivot_5b.shape[1]),
    "unclassified_terms": unclassified_5b,
    "significant_functional_cells": significant_cells(pivot_5b, padj_5b),
}

# ---------------------------------------------------------------------------------------
# 4. Family level -- Figures 6A / 6B, manuscript lines 119-135
# ---------------------------------------------------------------------------------------
t005, t01, gc = tables["families"]
family_classification = load_classification(
    "GO_terms_families_classification_Gemini - families.csv")
pivot_6b, padj_6b, unclassified_6b = functional_group_counts(
    t005, family_classification, "family_name")
family_to_class = (t005.drop_duplicates("family_name")
                   .set_index("family_name")["class_name"].to_dict())
counts_by_family = t005["family_name"].value_counts()
enrichment_families = pd.read_csv(
    os.path.join(REPO, "enrichment_families_with_random.csv"))

group_frequency = pivot_6b.sum(axis=0).sort_values(ascending=False)
results["families"] = {
    "n_families_with_terms": int(t005["family_name"].nunique()),
    "n_families_total": int(len(enrichment_families)),
    "families_lost_at_005": totals["families"]["groups_lost_at_005"],
    "families_by_class": {
        cls: sorted(f for f, c in family_to_class.items() if c == cls)
        for cls in sorted(set(family_to_class.values()))},
    "terms_per_family": {k: int(v) for k, v in counts_by_family.items()},
    "families_with_single_term": sorted(counts_by_family[counts_by_family == 1].index),
    "terms": {f: terms_of(t005, "family_name", f) for f in sorted(counts_by_family.index)},
    "terms_lost_at_005": {f: lost_terms(t005, t01, "family_name", f)
                          for f in sorted(set(t01["family_name"]))},
    "n_ltr_family_terms": int(
        t005.loc[t005["class_name"] == "LTR", "Term Name"].nunique()),
    "n_functional_groups": int(pivot_6b.shape[1]),
    "unclassified_terms": unclassified_6b,
    "functional_group_frequency": {k: int(v) for k, v in group_frequency.items()},
    "functional_groups_found_once": sorted(
        group_frequency[group_frequency == 1].index.tolist()),
    "significant_functional_cells": significant_cells(pivot_6b, padj_6b),
}
# The LTR class-level comparison the manuscript makes at line 123.
divergence_005 = tables["classes_divergence"][0]
class_count_005 = tables["classes_count"][0]
results["families"]["n_ltr_class_terms"] = int(
    class_count_005.loc[class_count_005["class_name"] == "LTR", "Term Name"].nunique())

# ---------------------------------------------------------------------------------------
# 5. Supplementary Figure 8B -- the combined class + family map, manuscript line 137
#
# The published caption defines S8B as one clustermap whose rows are the class-level groups
# (TE top, TE bottom, the classes with terms) *and* the family-level groups. Building it that
# way at FDR 0.1 reproduces every claim in line 137, which is what identifies the panel.
# ---------------------------------------------------------------------------------------
combined_classification = pd.concat(
    [family_classification, class_classification]).drop_duplicates("Term Name")


def combined_table(threshold):
    classes = read_output(f"GO_classes_count_{threshold}.csv")[["Term Name", "class_name"]]
    families = read_output(f"GO_families_{threshold}.csv")[["Term Name", "family_name"]]
    return pd.concat([classes.rename(columns={"class_name": "te_group"}),
                      families.rename(columns={"family_name": "te_group"})],
                     ignore_index=True)


s8b = {}
for label, threshold in [("fdr005", "fdr005"), ("fdr01", "fdr01_reference")]:
    pivot, padj, unclassified = functional_group_counts(
        combined_table(threshold), combined_classification, "te_group")
    s8b[label] = {
        "n_te_groups": int(pivot.shape[0]),
        "n_functional_groups": int(pivot.shape[1]),
        "unclassified_terms": unclassified,
        "significant_functional_cells": significant_cells(pivot, padj),
    }
    if label == "fdr005":
        pivot_s8b, padj_s8b, unclassified_s8b = pivot, padj, unclassified
results["supplementary_8b"] = s8b

# The subfamily-level clustermap is a companion-paper panel, kept here for the record only.
t005, t01, gc = tables["subfamilies"]
subfamily_classification = load_classification(
    "Subfamilies_classified - Classification.csv")
pivot_sub, padj_sub, unclassified_sub = functional_group_counts(
    t005, subfamily_classification, "subfamily_name")
results["subfamilies_companion_panel"] = {
    "n_subfamilies_with_terms": int(t005["subfamily_name"].nunique()),
    "n_terms": int(len(t005)),
    "n_functional_groups": int(pivot_sub.shape[1]),
    "unclassified_terms": unclassified_sub,
}

# ---------------------------------------------------------------------------------------
# 6. Figure 7 -- what determines a family's GO term count, at FDR 0.05
# ---------------------------------------------------------------------------------------
families_005 = tables["families"][0]
enrichment_families = enrichment_families.set_index("family_name", drop=False)
enrichment_families["Status"] = np.where(
    enrichment_families["p_adjusted_empirical_bh"] < FDR, "Significant",
    "Non-Significant")
enrichment_families["GO_terms_count"] = families_005.value_counts("family_name")
enrichment_families["GO_status"] = np.where(
    enrichment_families["GO_terms_count"].isna(), "No", "Yes")

# Panel letters follow the published Figure 7 caption, which is entirely family-level:
# A obs/random OR, B GO count by enrichment significance, C total copy number, D average
# divergence, E divergence by GO status, F obs/random OR by GO status, G copy number by GO
# status, H the class -> functional group -> family Sankey.
figure7 = {}
for panel, column, label in [
        ("7A", "OR_Observed_to_Random", "observed / random OR"),
        ("7C", "Total_TEs_number_log", "log10(total elements in a family)"),
        ("7D", "Average_divergence_all", "average divergence in a family")]:
    frame = enrichment_families.dropna(subset=[column, "GO_terms_count"])
    result = stats.pearsonr(frame[column], frame["GO_terms_count"])
    figure7[panel] = {"x": label, "n": int(len(frame)),
                      "pearson_r": float(result.statistic),
                      "pearson_p_raw": float(result.pvalue)}
for panel, split, value, label in [
        ("7B", "Status", "GO_terms_count", "GO terms by enrichment significance"),
        ("7E", "GO_status", "Average_divergence_all", "divergence by GO status"),
        ("7F", "GO_status", "OR_Observed_to_Random", "obs/random OR by GO status"),
        ("7G", "GO_status", "Total_TEs_number_log", "copy number by GO status")]:
    frame = enrichment_families.dropna(subset=[value])
    groups = [g for g in ["Significant", "Non-Significant", "Yes", "No"]
              if g in set(frame[split])]
    a = frame.loc[frame[split] == groups[0], value]
    b = frame.loc[frame[split] == groups[1], value]
    figure7[panel] = {
        "comparison": label, "groups": groups, "n": [int(len(a)), int(len(b))],
        "median": [float(a.median()), float(b.median())],
        "mannwhitney_p_raw": float(stats.mannwhitneyu(a, b).pvalue)}
results["figure7"] = figure7
results["figure7"]["families_with_zero_terms"] = int(
    (enrichment_families["GO_status"] == "No").sum())
results["figure7"]["zero_term_families_by_class"] = {
    cls: int(n) for cls, n in
    enrichment_families.loc[enrichment_families["GO_status"] == "No", "class_name"]
    .value_counts().items()}

# ---------------------------------------------------------------------------------------
# 7. IFNA domain -- Figure 8, Results new paragraph
# ---------------------------------------------------------------------------------------
ifna_tests = read_output("ifna_test_results.csv")
with open(os.path.join(OUT, "ifna_qc.json")) as handle:
    ifna_qc = json.load(handle)
results["ifna"] = {
    "qc": ifna_qc,
    "tests": [{k: (None if pd.isna(v) else v) for k, v in row.items()}
              for row in ifna_tests.to_dict("records")],
}

# ---------------------------------------------------------------------------------------
# 8. Window and percentile sensitivity -- Results robustness paragraph, Methods
# ---------------------------------------------------------------------------------------
concordance = read_output("window_sensitivity_concordance.csv")
flips = read_output("window_sensitivity_flips.csv")
percentile = read_output("percentile_sensitivity_summary.csv")
headline = read_output("percentile_sensitivity_headline.csv")
geneset = read_output("robustness_geneset_stability.csv")
ranks = read_output("robustness_rank_stability.csv")
with open(os.path.join(OUT, "window_sensitivity_ntss.json")) as handle:
    ntss = json.load(handle)

preserved = percentile["fraction_of_5pct_preserved"].dropna()
results["sensitivity"] = {
    "n_tss_bp_by_window": ntss,
    "concordance": concordance.to_dict("records"),
    "spearman_rho_range": [float(concordance["spearman_rho"].min()),
                           float(concordance["spearman_rho"].max())],
    "max_concordance_permutation_p": float(concordance["concordance_permutation_p"].max()),
    "n_significance_flips": int(len(flips)),
    "flip_groups": flips.iloc[:, 1].tolist(),
    "percentile_fraction_preserved_median": float(preserved.median()),
    "percentile_fraction_preserved_range": [float(preserved.min()),
                                            float(preserved.max())],
    "percentile_terms_lost_total": int(percentile["n_lost_at_10pct"].sum()),
    "percentile_terms_gained_total": int(percentile["n_gained_at_10pct"].sum()),
    "headline_claims": headline.to_dict("records"),
    "geneset_overlap_coefficient_range": [
        float(geneset["overlap_coefficient"].min()),
        float(geneset["overlap_coefficient"].max())],
    "geneset_max_hypergeometric_p": float(geneset["hypergeometric_p"].max()),
    "kendall_tau_range": [float(ranks["kendall_tau"].min()),
                          float(ranks["kendall_tau"].max())],
}

# ---------------------------------------------------------------------------------------
# 9. Permutation convergence -- Methods N = 500 justification
# ---------------------------------------------------------------------------------------
checkpoints = read_output("permutation_convergence_checkpoints.csv")
results["permutation_convergence"] = {
    "columns": list(checkpoints.columns),
    "checkpoints": checkpoints.to_dict("records"),
}

# ---------------------------------------------------------------------------------------
# 10. Lu et al. 2020 gene-set overlap -- Results comparison paragraph
# ---------------------------------------------------------------------------------------
lu_matrix = read_output("lu2020_overlap_matrix.csv")
lu_mapping = read_output("lu2020_mapping_summary.csv")
headline_columns = ["our_set", "their_category", "n_shared", "fold_over_expected",
                    "fisher_odds_ratio", "jaccard", "fisher_p_adjusted", "significant",
                    "direction", "n_ours_in_universe", "n_theirs_in_universe"]
results["lu2020"] = {
    "columns": list(lu_matrix.columns),
    "mapping_summary": lu_mapping.to_dict("records"),
    "counterpart_pairs": lu_matrix.loc[
        lu_matrix["our_set"].isin(["L1 (family)", "Alu (family)", "MIR (family)", "LINE", "SINE", "LTR", "DNA", "SVA"]),
        headline_columns].sort_values("fisher_p_adjusted").to_dict("records"),
    "all_overlaps": lu_matrix.to_dict("records"),
}

# ---------------------------------------------------------------------------------------
# 11. Newly resolved regions -- the bounded assembly argument in the Discussion
# ---------------------------------------------------------------------------------------
with open(os.path.join(OUT, "assembly_bound_summary.json")) as handle:
    results["assembly_bound"] = json.load(handle)

# ---------------------------------------------------------------------------------------
# 12. The GO grid -- every number manuscript edit M4 introduces
# ---------------------------------------------------------------------------------------
grid_dir = os.path.join(OUT, "GO_grid")
if os.path.isdir(grid_dir) and os.path.exists(os.path.join(grid_dir, "INDEX.csv")):
    grid_index = pd.read_csv(os.path.join(grid_dir, "INDEX.csv"))
    grid_preservation = read_output("go_grid_preservation.csv")
    grid_concordance = read_output("go_grid_concordance.csv")
    grid_headline = read_output("go_grid_headline_by_condition.csv")
    grid_survival = read_output("go_grid_group_survival.csv")
    with open(os.path.join(OUT, "go_grid_summary.json")) as handle:
        grid_summary = json.load(handle)

    published = grid_preservation["is_published"]
    alternatives = grid_preservation[~published]
    results["go_grid"] = {
        "n_cells": int(len(grid_index)),
        "windows": sorted(grid_index["window"].unique().tolist()),
        "percentiles": sorted(int(p) for p in grid_index["percentile"].unique()),
        # INDEX counts rows, i.e. one per (group, term); go_grid_preservation counts
        # DISTINCT terms. Both are reported because the manuscript could quote either and
        # they are not the same number.
        "term_rows_per_cell": {
            f"{r.level}/{r.window}/p{r.percentile}": int(r.n_terms_005)
            for r in grid_index.itertuples()
        },
        "distinct_terms_per_cell": {
            f"{r.level}/{r.window}/p{r.percentile}": int(r.n_terms)
            for r in grid_preservation.itertuples()
        },
        "window_trend_at_5pct": grid_summary["window_trend_at_5pct"],
        "percentile_always_increases_term_count":
            bool(grid_summary["percentile_always_increases"]),
        "fraction_of_published_preserved": {
            f"{r.level}/{r.window}/p{r.percentile}":
                round(float(r.fraction_of_published_preserved), 3)
            for r in alternatives.itertuples()
        },
        "fraction_of_published_preserved_min":
            round(float(alternatives["fraction_of_published_preserved"].min()), 3),
        "fraction_of_published_preserved_median":
            round(float(alternatives["fraction_of_published_preserved"].median()), 3),
        "jaccard_min": round(float(alternatives["jaccard"].min()), 3),
        "spearman_rho_min": round(float(grid_concordance["spearman_rho"].min()), 3),
        "spearman_rho_max": round(float(grid_concordance["spearman_rho"].max()), 3),
        "concordance_permutation_p_max":
            float(grid_concordance["concordance_permutation_p"].max()),
        "headline_claims_total": int(len(grid_headline)),
        "headline_claims_robust_all_six":
            int((grid_headline["n_conditions_survived"] == 6).sum()),
        "headline_claims_surviving_published":
            int(grid_headline["survives_published"].sum()),
        "headline_claims_by_conditions_survived": {
            f"{r.group}: {r.term_name}": int(r.n_conditions_survived)
            for r in grid_headline.itertuples()
        },
        "groups_surviving_all_six": grid_summary["groups_surviving_all_six"],
        "groups_total": grid_summary["groups_total"],
        # Caveat S8: the grid is GO only. A difference between cells is a gene-set effect,
        # not a background effect -- no permutations were re-run.
        "permutations_rerun_for_the_grid": False,
    }
else:
    results["go_grid"] = {"status": "not run — 07b_go_grid.py and nb07 have not been run"}

# ---------------------------------------------------------------------------------------
# 13. The supplementary package -- so the Supplementary material section can be checked
#     against what the five workbooks actually contain
# ---------------------------------------------------------------------------------------
inventory_path = os.path.join(HERE, "supplementary", "INVENTORY.json")
if os.path.exists(inventory_path):
    with open(inventory_path) as handle:
        inventory = json.load(handle)
    results["supplementary"] = {
        "n_workbooks": len(inventory),
        "workbooks": {
            filename: {
                "subject": record["subject"],
                "size_mb": round(record["size_bytes"] / 1e6, 2),
                "n_sheets": len(record["sheets"]),
                "sheets": {name: sheet["rows"] for name, sheet in record["sheets"].items()},
                "columns_dropped_for_excel": record["columns_dropped_for_excel"],
            }
            for filename, record in inventory.items()
        },
        "total_size_mb": round(
            sum(r["size_bytes"] for r in inventory.values()) / 1e6, 2),
    }
else:
    results["supplementary"] = {"status": "not built — run 08_build_supplementary.py"}

# ---------------------------------------------------------------------------------------
# Write
# ---------------------------------------------------------------------------------------
json_path = os.path.join(OUT, "results_numbers.json")
with open(json_path, "w") as handle:
    json.dump(results, handle, indent=2, default=str)


def render(value, indent=0):
    pad = "  " * indent
    if isinstance(value, dict):
        for key, item in value.items():
            if isinstance(item, (dict, list)):
                lines.append(f"{pad}{key}:")
                render(item, indent + 1)
            else:
                lines.append(f"{pad}{key}: {item}")
    elif isinstance(value, list):
        for item in value:
            if isinstance(item, dict):
                lines.append(f"{pad}- " + "; ".join(f"{k}={v}" for k, v in item.items()))
            else:
                lines.append(f"{pad}- {item}")
    else:
        lines.append(f"{pad}{value}")


lines = ["Numbers quoted by the revised manuscript, re-derived from revision_G3/output/",
         "=" * 78, ""]
for section, payload in results.items():
    lines += [f"### {section}", ""]
    render(payload, 1)
    lines.append("")
text_path = os.path.join(OUT, "results_numbers.txt")
with open(text_path, "w") as handle:
    handle.write("\n".join(lines) + "\n")

print(f"wrote {json_path}")
print(f"wrote {text_path}")
print()
print(f"GO terms at FDR {FDR}: "
      + ", ".join(f"{k} {v['n_terms_fdr005']} (was {v['n_terms_fdr01']})"
                  for k, v in totals.items()))
print(f"families with >= 1 term: {results['families']['n_families_with_terms']} of "
      f"{results['families']['n_families_total']}; lost: "
      f"{results['families']['families_lost_at_005']}")
print(f"unclassified terms: 4B {unclassified_4b}, 5B {unclassified_5b}, "
      f"6B {unclassified_6b}, S8B {unclassified_s8b}")
