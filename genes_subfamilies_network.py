import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import matplotlib.colors as mcolors
from adjustText import adjust_text
from matplotlib.lines import Line2D
from tqdm import tqdm

plt.rcParams["svg.fonttype"] = "none"
plt.rcParams.update({"font.size": 10})
# --- 1. Global Configurations (Tune these to build optimal structure) ---
save_path = "plots/TE_subfamilies_clustered_connectivity_network_less_repulsion1_jaccard_0025_divergence.svg"
COLUMN_NAME = "Average_divergence_all"

# Text and Data limits
GLOBAL_FONT_SIZE = 10
SUBFAMILY_LIMIT = 1143 #len(sets_genes) # Process all subfamilies
JACCARD_THRESHOLD = 0.025

# Edge visual parameters
BASE_WIDTH = 0.3
WIDTH_MULTIPLIER = 15

# Layout Dynamics (Force-directed parameters)
REPULSION_K_FACTOR = 1.0   #300  # Higher = stronger push between nodes in clusters
LAYOUT_ITERATIONS = 500      # More iterations = more stable cluster separation
COMPONENT_SPACING = 7.0      # Physical distance between unrelated clusters
ISOLATE_RADIUS_SCALE = 1.3   # Factor to push zero-connected nodes into the outer ring
JITTER_STRENGTH = 0.15       # Intensity of random randomness added to positions

# Collision & Overlap Management
COLLISION_MIN_DIST = 110     # Overlap penalty (multiplier for node radii sum)
COLLISION_ITERATIONS = 300   # Intensity of the physics relaxation step
NODE_RADIUS_DIVISOR = 450    # Controls the coordinate-space "size" of nodes

class_names_p = {
    "LINE": "#cc660b",
    "LTR": "#70453c",
    "SINE": "#ab1f20",
    "DNA": "#195f90",
    "Retroposon": "#765297",
    "RC": "#238023",
}


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

for subfamily_name in tqdm(subfamilies_by_classes['subfamily_name']):

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

enrichment_df_subfamilies.index = enrichment_df_subfamilies['subfamily_name']


# --- 2. Calculation Functions ---

def jaccard_index(set1, set2):
    """Calculates the Jaccard similarity index between two sets."""
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    return intersection / union if union > 0 else 0

def prevent_collisions(pos, node_sizes, iterations, min_dist_factor, jitter):
    """
    Manually relaxes node positions to prevent 3+ nodes from overlapping
    and to ensure white space borders between clusters.
    """
    new_pos = {n: np.array(p) for n, p in pos.items()}
    nodes = list(new_pos.keys())
    
    # Map node size (area) to approximate radius in coordinate space
    radii = {n: np.sqrt(node_sizes[i]) / NODE_RADIUS_DIVISOR for i, n in enumerate(nodes)}
    
    for _ in range(iterations):
        # Apply slight jitter at the start of each relaxation step
        if jitter > 0:
            for n in nodes:
                new_pos[n] += np.random.normal(0, jitter * 0.1, 2)

        for i in range(len(nodes)):
            for j in range(i + 1, len(nodes)):
                n1, n2 = nodes[i], nodes[j]
                vec = new_pos[n1] - new_pos[n2]
                dist = np.linalg.norm(vec)
                
                target_dist = (radii[n1] + radii[n2]) * min_dist_factor
                
                if dist < target_dist and dist > 0:
                    # Push nodes apart with force proportional to overlap
                    push = (target_dist - dist) / dist * 0.5 * vec
                    new_pos[n1] += push
                    new_pos[n2] -= push
                elif dist == 0:
                    # Resolve identical positions with random nudge
                    new_pos[n1] += np.random.normal(0, 0.05, 2)
    return new_pos

# --- 3. Graph Construction ---

G = nx.Graph()

# Enrichment bounds calculation with epsilon for log-safety
eps = 1e-3
log_or_series = np.log10(enrichment_df_subfamilies[COLUMN_NAME] + eps)
min_log_or = log_or_series.min()
max_log_or = log_or_series.max()

# Keep raw bounds for legend labels
min_or = enrichment_df_subfamilies[COLUMN_NAME].min()
max_or = enrichment_df_subfamilies[COLUMN_NAME].max()

subset_labels = labels_genes[:SUBFAMILY_LIMIT]
subset_sets = sets_genes[:SUBFAMILY_LIMIT]

print(f"Building graph for {len(subset_labels)} subfamilies...")

# Add nodes with attributes
for i, (label, gene_set) in enumerate(zip(subset_labels, subset_sets)):
    cls = subfamily_to_class.get(label, 'Other')
    base_color = class_names_p.get(cls, '#808080')
    
    # Robust lookup for enrichment data
    if label in enrichment_df_subfamilies.index:
        target_row = enrichment_df_subfamilies.loc[label]
    elif 'subfamily_name' in enrichment_df_subfamilies.columns and (enrichment_df_subfamilies['subfamily_name'] == label).any():
        target_row = enrichment_df_subfamilies[enrichment_df_subfamilies['subfamily_name'] == label].iloc[0]
    else:
        target_row = None

    if target_row is not None:
        or_val = target_row[COLUMN_NAME]
        status = target_row['Status']
    else:
        or_val = min_or
        status = 'Non-Significant'
    
    # Log-scaled alpha mapping
    current_log_or = np.log10(or_val + eps)
    if max_log_or != min_log_or:
        alpha = 0.15 + 0.85 * (current_log_or - min_log_or) / (max_log_or - min_log_or)
    else:
        alpha = 1.0
    alpha = np.clip(alpha, 0.15, 1.0)
        
    rgb = mcolors.to_rgb(base_color)
    rgba = (*rgb, alpha)
    
    size = len(gene_set)
    log_size = np.log10(size + 1) * 450 
    
    G.add_node(label, 
               size=log_size, 
               color=rgba, 
               class_name=cls, 
               status=status, 
               or_val=or_val)

# Edge calculation
for i in range(len(subset_labels)):
    for j in range(i + 1, len(subset_labels)):
        j_idx = jaccard_index(subset_sets[i], subset_sets[j])
        if j_idx > JACCARD_THRESHOLD:
            G.add_edge(subset_labels[i], subset_labels[j], weight=j_idx)

# --- 4. Refactored Connectivity-based Layout ---

print("Calculating Cluster and Connectivity-aware Layout...")

# Find connected components to isolate clusters
components = list(nx.connected_components(G))
components.sort(key=len, reverse=True)

pos = {}
num_cols = int(np.ceil(np.sqrt(len(components))))

# Separate isolates from connected clusters
isolates = [list(c)[0] for c in components if len(c) == 1]
connected_clusters = [c for c in components if len(c) > 1]

# Layout connected clusters
for i, component in enumerate(connected_clusters):
    # Cluster center based on spacing parameter
    center_x = (i % num_cols) * COMPONENT_SPACING
    center_y = (i // num_cols) * COMPONENT_SPACING
    
    S = G.subgraph(component)
    # k parameter is multiplied by the REPULSION_K_FACTOR to push nodes apart
    local_pos = nx.spring_layout(S, 
                                 k=REPULSION_K_FACTOR/np.sqrt(len(component)), 
                                 iterations=LAYOUT_ITERATIONS, 
                                 weight='weight', 
                                 seed=42)
    
    # Translate and apply jitter
    for node, p in local_pos.items():
        rand_nudge = np.random.normal(0, JITTER_STRENGTH, 2)
        pos[node] = p + np.array([center_x, center_y]) + rand_nudge

# Layout isolates on a large periphery ring
if isolates:
    grid_extent = COMPONENT_SPACING * num_cols
    periphery_radius = grid_extent * ISOLATE_RADIUS_SCALE
    center_offset = (grid_extent / 2, grid_extent / 2)
    
    # Randomly shuffle angles to prevent "ordered" look
    iso_angles = np.linspace(0, 2 * np.pi, len(isolates), endpoint=False)
    np.random.shuffle(iso_angles)
    
    for node, angle in zip(isolates, iso_angles):
        # Vary the radius slightly for each isolate to introduce depth/jitter
        var_radius = periphery_radius * (1 + np.random.normal(0, JITTER_STRENGTH))
        pos[node] = center_offset + np.array([var_radius * np.cos(angle), var_radius * np.sin(angle)])

# Safely extract attributes for plotting
node_colors = [G.nodes[n]['color'] for n in G.nodes]
node_sizes = [G.nodes[n]['size'] for n in G.nodes]

print("Applying physics relaxation to enforce white space and resolve overlaps...")
pos = prevent_collisions(pos, 
                         node_sizes, 
                         iterations=COLLISION_ITERATIONS, 
                         min_dist_factor=COLLISION_MIN_DIST,
                         jitter=JITTER_STRENGTH)

# --- 5. Visualization ---

fig, ax = plt.subplots(figsize=(28, 28))

# Draw Edges
edges = G.edges(data=True)
edge_widths = [BASE_WIDTH + (d['weight'] * WIDTH_MULTIPLIER) for u, v, d in edges]
nx.draw_networkx_edges(G, pos, width=edge_widths, alpha=0.15, edge_color='gray', ax=ax)

# Draw Nodes
nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors, 
                       edgecolors='white', linewidths=0.4, ax=ax)

# Draw Labels (Significant only)
texts = []
for node, (x, y) in pos.items():
    status = G.nodes[node]['status']
    if status == 'Significant':
        texts.append(ax.text(x, y, node, fontsize=GLOBAL_FONT_SIZE, fontweight='bold', color='black'))

print("Adjusting label positions...")
adjust_text(
    texts, 
    ax=ax, 
    expand_points=(2.0, 2.0), 
    arrowprops=dict(arrowstyle='-', color='gray', lw=0.4, alpha=0.6),
    force_points=0.1
)

plt.title(f"Tunable Cluster Map of TE Subfamilies\n(Connectivity hubs central, Jaccard > {JACCARD_THRESHOLD}, N={len(G.nodes)})", 
          fontsize=GLOBAL_FONT_SIZE + 4, pad=40, weight='bold')

# --- 6. Legends ---

# A) Node Size Legend
legend_sizes = [10, 100, 1000]
legend_nodes = []
for size in legend_sizes:
    log_s = np.log10(size + 1) * 450
    legend_nodes.append(Line2D([0], [0], marker='o', color='w', label=f'{size} genes',
                               markerfacecolor='gray', markersize=np.sqrt(log_s)))

size_legend = ax.legend(handles=legend_nodes, title="Gene Set Size", 
                         loc='upper left', bbox_to_anchor=(1, 1), frameon=False, 
                         labelspacing=1.5, title_fontsize=GLOBAL_FONT_SIZE, fontsize=GLOBAL_FONT_SIZE)
ax.add_artist(size_legend)

# B) Edge Width Legend
legend_jaccard_vals = [0.05, 0.1, 0.2]
legend_edges = []
for j in legend_jaccard_vals:
    lw = BASE_WIDTH + (j * WIDTH_MULTIPLIER)
    legend_edges.append(Line2D([0], [0], color='gray', lw=lw, alpha=0.4, label=f'J={j:.2f}'))

jaccard_legend = ax.legend(handles=legend_edges, title="Jaccard Index", 
                           loc='lower left', bbox_to_anchor=(1, 0), frameon=False, 
                           labelspacing=1.5, title_fontsize=GLOBAL_FONT_SIZE, fontsize=GLOBAL_FONT_SIZE)
ax.add_artist(jaccard_legend)

# C) Class Color Legend
legend_classes = [Line2D([0], [0], marker='o', color='w', label=cls,
                         markerfacecolor=color, markersize=12) 
                  for cls, color in class_names_p.items()]

class_legend = ax.legend(handles=legend_classes, title="Genomic Class", 
                         loc='center left', bbox_to_anchor=(1, 0.7), frameon=False, 
                         labelspacing=1.2, title_fontsize=GLOBAL_FONT_SIZE, fontsize=GLOBAL_FONT_SIZE)
ax.add_artist(class_legend)

# D) Transparency Legend (Log-scaled)
legend_or_vals = [min_or, (min_or + max_or)/4, (min_or + max_or)/2, max_or]
legend_or = []
for val in legend_or_vals:
    current_log = np.log10(val + eps)
    if max_log_or != min_log_or:
        l_alpha = 0.15 + 0.85 * (current_log - min_log_or) / (max_log_or - min_log_or)
    else:
        l_alpha = 1.0
    l_alpha = np.clip(l_alpha, 0.15, 1.0)
    legend_or.append(Line2D([0], [0], marker='o', color='w', label=f'{val:.2f}',
                            markerfacecolor=(0.5, 0.5, 0.5, l_alpha), markersize=12))

or_legend = ax.legend(handles=legend_or, title="Enrichment (OR)", 
                      loc='center left', bbox_to_anchor=(1, 0.35), frameon=False, 
                      labelspacing=1.2, title_fontsize=GLOBAL_FONT_SIZE, fontsize=GLOBAL_FONT_SIZE)

plt.tight_layout(rect=[0, 0, 0.85, 1])
ax.axis('off')

# --- 7. Saving ---

os.makedirs('plots', exist_ok=True)

plt.savefig(
    save_path, 
    format='svg', 
    transparent=True, 
    bbox_inches='tight'
)

plt.show()

print(f"Success: Refactored network with global tuning parameters saved to {save_path}.")
