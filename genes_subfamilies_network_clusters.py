import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import matplotlib.colors as mcolors
from adjustText import adjust_text
from matplotlib.lines import Line2D
from tqdm import tqdm
from sklearn.metrics import silhouette_score
from scipy.sparse import issparse

# Set font type for clean SVG export
plt.rcParams["svg.fonttype"] = "none"

# --- 1. Global Configurations (Tune these to build optimal structure) ---

# Text and Data limits
GLOBAL_FONT_SIZE = 10
SUBFAMILY_LIMIT = 1143  # Limited for performance during clustering & layout
JACCARD_THRESHOLD = 0.01

# MCL Clustering Parameters
# Inflation controls cluster granularity. Higher = more/smaller clusters.
MCL_INFLATION = 1.14 
MCL_MAX_ITER = 100
MCL_BENCHMARK_RANGE = np.linspace(1, 1.5, 25) # For the performance plot

# Layout Dynamics (Force-directed parameters)
REPULSION_K_FACTOR = 4.0   #300  # Higher = stronger push between nodes in clusters
LAYOUT_ITERATIONS = 500      # More iterations = more stable cluster separation
CLUSTER_SPACING = 30.0      # Physical distance between unrelated clusters
ISOLATE_RADIUS_SCALE = 1.5   # Factor to push zero-connected nodes into the outer ring
JITTER_STRENGTH = 0.15       # Intensity of random randomness added to positions

# Collision & Overlap Management
COLLISION_MIN_DIST = 210     # Overlap penalty (multiplier for node radii sum)
COLLISION_ITERATIONS = 300   # Intensity of the physics relaxation step
NODE_RADIUS_DIVISOR = 450    # Controls the coordinate-space "size" of nodes

BASE_WIDTH = 0.3
WIDTH_MULTIPLIER = 15


class_names_p = {
    "LINE": "#cc660b",
    "LTR": "#70453c",
    "SINE": "#ab1f20",
    "DNA": "#195f90",
    "Retroposon": "#765297",
    "RC": "#238023",
}

COLUMN_NAME = "Average_divergence_all"
TEs_on_genes_subfamilies = pd.read_csv("TEs_on_genes_counts_subfamilies.csv")

# gather top  5% genes by subfamilies and analyze them on GO
sets_genes = []
labels_genes = []

background_list = list(TEs_on_genes_subfamilies["Gene_name"].unique())

significant_terms_dfs = []

genes_to_take_number = len(TEs_on_genes_subfamilies["Gene_name"].unique()) // 20

subfamilies_by_classes = (
    pd.read_csv(
        "individuals_by_classes_TE.csv"
    )  # .dropna()#.set_index("subfamily_name")
)[["class_name", "individual_name"]]
subfamilies_by_classes.columns = ["class_name", "subfamily_name"]
subfamilies_by_classes.index = subfamilies_by_classes["subfamily_name"]
subfamily_to_class = subfamilies_by_classes["class_name"].to_dict()

for subfamily_name in tqdm(subfamilies_by_classes["subfamily_name"]):

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

    sets_genes.append(set(genes_to_take))

    labels_genes.append(subfamily_name)


enrichment_df_subfamilies = pd.read_csv(
    "enrichment_subfamilies_with_random.csv", index_col=0
).reset_index()
enrichment_df_subfamilies["Status"] = "Non-Significant"
enrichment_df_subfamilies["Status"][
    enrichment_df_subfamilies["p_adjusted_empirical_bh"] < 0.05
] = "Significant"
enrichment_df_subfamilies["OR_Observed_to_Random"] = (
    enrichment_df_subfamilies["Observed_Odds_Ratio"]
    / enrichment_df_subfamilies["Random_Odds_Ratio_Mean"]
)
enrichment_df_subfamilies["Enrichment_p_value_adjusted_log"] = np.log10(
    enrichment_df_subfamilies.Enrichment_p_value_adjusted
)

enrichment_df_subfamilies.index = enrichment_df_subfamilies["subfamily_name"]

# --- 2. Clustering & Calculation Functions ---

def run_mcl(matrix, inflation=2.0, max_iter=100):
    """NumPy-based Markov Clustering (MCL) implementation."""
    A = matrix + np.eye(matrix.shape[0])
    M = A / A.sum(axis=0)
    
    for i in range(max_iter):
        M_prev = M.copy()
        M = np.linalg.matrix_power(M, 2)
        M = np.power(M, inflation)
        M = M / M.sum(axis=0)
        if np.allclose(M, M_prev, atol=1e-5):
            break
            
    clusters = {}
    cluster_id = 0
    visited = set()
    for col in range(M.shape[1]):
        if col in visited: continue
        attractor = np.argmax(M[:, col])
        indices = np.where(np.isclose(M[attractor, :], M[attractor, col], atol=1e-3))[0]
        for idx in indices:
            clusters[idx] = cluster_id
            visited.add(idx)
        cluster_id += 1
    return clusters

def jaccard_index(set1, set2):
    """Calculates the Jaccard similarity index between two sets."""
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    return intersection / union if union > 0 else 0

def prevent_collisions(pos, node_sizes, iterations, min_dist_factor, jitter):
    """Relaxes node positions to prevent overlaps in coordinate space."""
    new_pos = {n: np.array(p) for n, p in pos.items()}
    nodes = list(new_pos.keys())
    radii = {n: np.sqrt(node_sizes[i]) / NODE_RADIUS_DIVISOR for i, n in enumerate(nodes)}
    
    for _ in range(iterations):
        if jitter > 0:
            for n in nodes:
                new_pos[n] += np.random.normal(0, jitter * 0.04, 2)
        for i in range(len(nodes)):
            for j in range(i + 1, len(nodes)):
                n1, n2 = nodes[i], nodes[j]
                vec = new_pos[n1] - new_pos[n2]
                dist = np.linalg.norm(vec)
                target_dist = (radii[n1] + radii[n2]) * min_dist_factor
                if dist < target_dist and dist > 0:
                    push = (target_dist - dist) / dist * 0.5 * vec
                    new_pos[n1] += push
                    new_pos[n2] -= push
    return new_pos

# --- 3. Graph Construction ---

G = nx.Graph()
subset_labels = labels_genes[:SUBFAMILY_LIMIT]
subset_sets = sets_genes[:SUBFAMILY_LIMIT]

adj_matrix = np.zeros((len(subset_labels), len(subset_labels)))
for i in range(len(subset_labels)):
    for j in range(i + 1, len(subset_labels)):
        j_idx = jaccard_index(subset_sets[i], subset_sets[j])
        if j_idx > JACCARD_THRESHOLD:
            G.add_edge(subset_labels[i], subset_labels[j], weight=j_idx)
            adj_matrix[i, j] = j_idx
            adj_matrix[j, i] = j_idx

for i, (label, gene_set) in enumerate(zip(subset_labels, subset_sets)):
    size = len(gene_set)
    log_size = np.log10(size + 1) * 450 
    if label not in G:
        G.add_node(label)
    G.nodes[label]['size'] = log_size

# --- 4. MCL Performance Benchmarking ---

print("Benchmarking MCL cluster performance...")
perf_stats = []
distance_matrix = 1.0 - adj_matrix

for infl in tqdm(MCL_BENCHMARK_RANGE, desc="Testing Inflation Values"):
    labels_dict = run_mcl(adj_matrix, inflation=infl, max_iter=MCL_MAX_ITER)
    labels_list = [labels_dict[i] for i in range(len(subset_labels))]
    unique_ids = set(labels_list)
    n_clusters = len(unique_ids)
    
    # Calculate unsupervised metrics
    try:
        sil = silhouette_score(distance_matrix, labels_list, metric='precomputed') if n_clusters > 1 else 0
    except:
        sil = 0
    
    comms = [set([subset_labels[idx] for idx, cid in labels_dict.items() if cid == target]) for target in unique_ids]
    mod = nx.community.modularity(G, comms) if n_clusters > 1 else 0
    
    # Added Coverage: Fraction of edges that are internal to clusters
    coverage, _ = nx.community.partition_quality(G, comms)
    
    perf_stats.append({
        'inflation': infl, 
        'clusters': n_clusters, 
        'silhouette': sil, 
        'modularity': mod,
        'coverage': coverage
    })

perf_df = pd.DataFrame(perf_stats)

# --- 5. Final MCL & Cluster-Aware Layout ---

print(f"Applying final MCL with inflation {MCL_INFLATION}...")
final_labels_dict = run_mcl(adj_matrix, inflation=MCL_INFLATION)
for idx, c_id in final_labels_dict.items():
    G.nodes[subset_labels[idx]]['cluster'] = c_id

# Group by Cluster ID
cluster_groups = {}
for idx, c_id in final_labels_dict.items():
    if c_id not in cluster_groups:
        cluster_groups[c_id] = []
    cluster_groups[c_id].append(subset_labels[idx])

# Cluster Color Setup
cmap = plt.get_cmap('tab20')
cluster_colors = {c: cmap(i % 20) for i, c in enumerate(cluster_groups.keys())}

# Calculate Layout: Cluster-by-Cluster
print("Calculating Cluster-Aware Layout...")
pos = {}
num_clusters = len(cluster_groups)
grid_cols = int(np.ceil(np.sqrt(num_clusters)))

# Sort clusters by size (largest first)
sorted_clusters = sorted(cluster_groups.items(), key=lambda x: len(x[1]), reverse=True)

for i, (c_id, node_list) in enumerate(sorted_clusters):
    # Cluster center on a global grid
    center_x = (i % grid_cols) * CLUSTER_SPACING
    center_y = (i // grid_cols) * CLUSTER_SPACING
    
    # Isolate clusters that are single nodes (isolates)
    if len(node_list) == 1:
        continue # Handled later in periphery ring
    
    S = G.subgraph(node_list)
    # Strong local layout for the cluster
    local_pos = nx.spring_layout(S, k=REPULSION_K_FACTOR/np.sqrt(len(node_list)), 
                                 iterations=LAYOUT_ITERATIONS, weight='weight', seed=42)
    
    for node, p in local_pos.items():
        pos[node] = p + np.array([center_x, center_y]) + np.random.normal(0, JITTER_STRENGTH, 2)

# Periphery layout for single-node clusters
isolates = [nodes[0] for c_id, nodes in cluster_groups.items() if len(nodes) == 1]
if isolates:
    extent = CLUSTER_SPACING * grid_cols
    radius = extent * ISOLATE_RADIUS_SCALE
    center = (extent/2, extent/2)
    angles = np.linspace(0, 2*np.pi, len(isolates), endpoint=False)
    for node, angle in zip(isolates, angles):
        pos[node] = center + np.array([radius * np.cos(angle), radius * np.sin(angle)]) + np.random.normal(0, JITTER_STRENGTH, 2)

node_sizes = [G.nodes[n]['size'] for n in subset_labels]
pos = prevent_collisions(pos, node_sizes, COLLISION_ITERATIONS, COLLISION_MIN_DIST, JITTER_STRENGTH)

# --- 6. Visualization ---

# Panel 1: Multi-Metric Performance Plot
fig_p, ax_p = plt.subplots(figsize=(10, 6))
ax_p.plot(perf_df['inflation'], perf_df['silhouette'], label='Silhouette (Separation)', marker='o', lw=2)
ax_p.plot(perf_df['inflation'], perf_df['modularity'], label='Modularity (Topology)', marker='s', lw=2)
ax_p.plot(perf_df['inflation'], perf_df['coverage'], label='Coverage (Edge Density)', marker='^', lw=2)

# Fix non-overlapping cluster bars
ax_p2 = ax_p.twinx()
bar_width = (MCL_BENCHMARK_RANGE[1] - MCL_BENCHMARK_RANGE[0]) * 0.4
ax_p2.bar(perf_df['inflation'], perf_df['clusters'], alpha=0.15, color='gray', width=bar_width, label='Cluster Count')

ax_p.set_xlabel('MCL Inflation Parameter', fontsize=GLOBAL_FONT_SIZE)
ax_p.set_ylabel('Performance Score', fontsize=GLOBAL_FONT_SIZE)
ax_p2.set_ylabel('Number of Clusters', fontsize=GLOBAL_FONT_SIZE)
ax_p.set_title('Hierarchical MCL Performance Benchmarking', fontsize=GLOBAL_FONT_SIZE+2, fontweight='bold')
ax_p.legend(loc='upper left', frameon=False)
ax_p2.legend(loc='upper right', frameon=False)
plt.tight_layout()
plt.savefig("plots/MCL_Hierarchical_Metrics.svg", format='svg')

# Panel 2: The Network Map
fig, ax = plt.subplots(figsize=(26, 26))
edges = G.edges(data=True)
edge_widths = [BASE_WIDTH + (d['weight'] * WIDTH_MULTIPLIER) for u, v, d in edges]
nx.draw_networkx_edges(G, pos, width=edge_widths, alpha=0.1, edge_color='gray', ax=ax)

# Node Colors
node_colors = [cluster_colors[G.nodes[n]['cluster']] for n in subset_labels]
nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors, 
                       edgecolors='white', linewidths=0.5, ax=ax)

# Labels (Significant or High Degree)
texts = []
for node, (x, y) in pos.items():
    if G.degree(node) >= 2:
        texts.append(ax.text(x, y, node, fontsize=GLOBAL_FONT_SIZE, fontweight='bold', color='black'))

print("Polishing label layout...")
adjust_text(texts, ax=ax, expand_points=(1.8, 1.8), 
            arrowprops=dict(arrowstyle='-', color='gray', lw=0.4, alpha=0.5), force_points=0.1)

plt.title(f"Markov Clustered TE Subfamily Map (Inflation={MCL_INFLATION}, N={num_clusters} Clusters)", 
          fontsize=GLOBAL_FONT_SIZE + 8, pad=45, weight='bold')

# Legend: TOP 20 CLUSTERS
top_20 = sorted(cluster_groups.keys(), key=lambda x: len(cluster_groups[x]), reverse=True)[:20]
legend_els = [Line2D([0], [0], marker='o', color='w', label=f'Cluster {c} (n={len(cluster_groups[c])})', 
                     markerfacecolor=cluster_colors[c], markersize=12) for c in top_20]
ax.legend(handles=legend_els, title="Top 20 MCL Clusters", loc='upper left', bbox_to_anchor=(1.02, 1), 
          frameon=False, title_fontsize=GLOBAL_FONT_SIZE+2, fontsize=GLOBAL_FONT_SIZE)

ax.axis('off')
plt.tight_layout(rect=[0, 0, 0.85, 1])

# --- 7. Saving ---
os.makedirs('plots', exist_ok=True)
save_path_net = "plots/TE_subfamilies_MCL_clustered_map.svg"
plt.savefig(save_path_net, format='svg', transparent=True, bbox_inches='tight')

plt.show()

print(f"Success: Cluster-aware network saved to {save_path_net}")
print(f"Success: Benchmark plot saved to plots/MCL_Hierarchical_Metrics.svg")
