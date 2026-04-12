import warnings
import io
import itertools
import json
import os
import re
import time  # Added for performance comparison
import urllib.parse
from collections import defaultdict

import matplotlib
import matplotlib.colors as mcolors
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import requests
import seaborn as sns
from adjustText import adjust_text
from goatools.anno.gaf_reader import GafReader
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.obo_parser import GODag
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.path import Path
from matplotlib.ticker import LogLocator, NullFormatter, ScalarFormatter
from pyvis.network import Network
from scipy import stats
from scipy.stats import fisher_exact, kstest
from statannotations.Annotator import Annotator
from statsmodels.stats.multitest import multipletests
from supervenn import supervenn
from tqdm import tqdm

plt.rcParams["svg.fonttype"] = "none"
plt.rcParams.update({"font.size": 10})

def run_goatools_enrichment(
    gene_list,
    background_list,
    obo_path="go-basic.obo",
    gaf_path="goa_human.gaf",
    fdr_threshold=0.1,
):
    # 1. Load the Ontology DAG
    print("Loading GO DAG...")
    godag = GODag(obo_path)

    # 2. Load Human Annotations and Map Symbols manually
    print("Loading Human Annotations and mapping symbols...")
    ogaf = GafReader(gaf_path)

    # Extract associations: Symbol -> Set of GO IDs
    full_assoc = {}
    for ntf in ogaf.associations:
        symbol = ntf.DB_Symbol
        if symbol not in full_assoc:
            full_assoc[symbol] = set()
        full_assoc[symbol].add(ntf.GO_ID)

    # 3. Initialize Enrichment Study
    print("Initializing Enrichment Study...")
    goeaobj = GOEnrichmentStudy(
        background_list,
        full_assoc,
        godag,
        propagate_counts=False,
        alpha=fdr_threshold,
        methods=["fdr_bh"],
    )

    # 4. Run Study
    print(f"Running Analysis on {len(gene_list)} genes...")
    results_all = goeaobj.run_study(gene_list)

    # Filter for significant results (FDR )
    results_sig = [r for r in results_all if r.p_fdr_bh < fdr_threshold]

    # Use the precomputed population total from the study object keys you provided
    pop_total = goeaobj.pop_n

    parsed_results = []
    for res in results_sig:
        # Use the precomputed counts from the record keys you provided
        pop_count = res.pop_count
        study_count = res.study_count

        fold_enrichment = (study_count / res.study_n) / (pop_count / res.pop_n)

        parsed_results.append(
            {
                "Term ID": res.GO,
                "Term Name": res.name,
                "Term Database": res.NS,
                "P-value": res.p_uncorrected,
                "FDR": res.p_fdr_bh,  # FDR correctly distinct from P-value
                "Fold Enrichment": fold_enrichment,
                "Overlap Count": study_count,
                "Total Term Genes (Human)": pop_count,
                # Use precomputed lists from the record keys you provided
                "Overlapping Genes": sorted(res.study_items),
                "Full Term Gene List": sorted(res.pop_items),
            }
        )

    df = pd.DataFrame(parsed_results)
    if not df.empty:
        # Sort by adjusted p-value as requested
        return df.sort_values("FDR")
    else:
        print("No significant terms passed the filters.")
        return pd.DataFrame()

TEs_on_genes_subfamilies = pd.read_csv('TEs_on_genes_counts_subfamilies.csv')

# gather top  5% genes by subfamilies and analyze them on GO
sets_genes = []
labels_genes = []

background_list = list(TEs_on_genes_subfamilies["Gene_name"].unique())

significant_terms_dfs = []

genes_to_take_number = len(TEs_on_genes_subfamilies["Gene_name"].unique()) // 20

subfamilies_by_classes = (
    pd.read_csv("individuals_by_classes_TE.csv")#.dropna()#.set_index("subfamily_name")
)[['class_name', 'individual_name']]
subfamilies_by_classes.columns = ['class_name', 'subfamily_name']
subfamilies_by_classes.index = subfamilies_by_classes['subfamily_name']
subfamily_to_class = subfamilies_by_classes["class_name"].to_dict()

for subfamily_name in tqdm(subfamilies_by_classes['subfamily_name'][::-1]):
    print(
        f"""
    {subfamily_name}
    
    """
    )
    genes_by_TE_number = list(
        TEs_on_genes_subfamilies[
            TEs_on_genes_subfamilies[f"{subfamily_name}_number"] != 0
        ]
        .sort_values(f"{subfamily_name}_number", ascending=False)[["Gene_name"]]
        .drop_duplicates()["Gene_name"]
    )
    if len(genes_by_TE_number) >= genes_to_take_number:
        genes_to_take = genes_by_TE_number[:genes_to_take_number]
    elif len(genes_by_TE_number) > 0:
        genes_to_take = genes_by_TE_number
    else:
        continue

#    sets_genes.append(set(genes_to_take))

#    labels_genes.append(subfamily_name)

    results_GO = run_goatools_enrichment(genes_to_take, background_list=background_list)
    results_GO["subfamily_name"] = subfamily_name
    results_GO["class_name"] = subfamily_to_class[subfamily_name]
    if len(results_GO) > 0:
        results_GO.to_csv(f'GO_tables/GO_terms_by_subfamilies_{subfamily_name}.csv')
#    significant_terms_dfs.append(results_GO)


#top_5_perc_genes_by_subfamilies = pd.concat(significant_terms_dfs, axis=0)

#top_5_perc_genes_by_subfamilies.to_csv('top_5_perc_genes_by_subfamilies.csv')
