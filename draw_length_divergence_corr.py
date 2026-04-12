import itertools
import os
import time  # Added for performance comparison
from collections import defaultdict

import matplotlib
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.patches import Patch
from scipy import stats
from scipy.stats import fisher_exact, kstest
from statannotations.Annotator import Annotator
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm

class_names = ["LINE", "LTR", "SINE", "DNA", "Retroposon", "RC"]

class_names_p = {
    "LINE": "#cc660b",
    "LTR": "#70453c",
    "SINE": "#ab1f20",
    "DNA": "#195f90",
    "Retroposon": "#765297",
    "RC": "#238023",
}
repeat_masker = pd.read_csv("../T2T_article/T2T_repeat_masker_processed.csv")
repeat_masker.columns = [
    "Unnamed: 0",
    '#"chrom"',
    "chromStart",
    "chromEnd",
    "score",
    "subfamily_name",
    "family_name",
    "class_name",
]
repeat_masker['length'] = repeat_masker['chromEnd'] - repeat_masker['chromStart']

column1 = "score"
column2 = "length"

figure, axs = plt.subplots(
    nrows=2,
    ncols=3,
    sharex=True,
    sharey=False,
    squeeze=True,
    # width_ratios=None,
    # height_ratios=[1, 1, 1, 1, 1, 1, 0.3, 0.3],
    figsize=(14, 8),
)


for class_name, ax in zip(class_names, axs.flatten()):
    
    data_to_plot_local = repeat_masker.query(f"class_name == '{class_name}'")
    corr_data = data_to_plot_local.dropna(subset=[column1, column2])
    corr_result = stats.pearsonr(corr_data[column1], corr_data[column2])
    corr_r = corr_result.statistic

    corr_p = corr_result.pvalue
    print(class_name, corr_r, corr_p)
    # Format the statistics string
    if corr_p < 0.001:
        p_val_str = "p < 0.001"
    else:
        p_val_str = f"p = {corr_p:.3f}"

    stats_text = f"Pearson R: {corr_r:.3f}\n{p_val_str}"

    # Use ax.set_title() to set the title for the subplot
    ax.set_title(class_name)

    # Draw the scatterplot onto the axis
    sns.scatterplot(
        data=data_to_plot_local,
        hue="class_name",
        x=column1,
        y=column2,
        palette=class_names_p,
        ax=ax,
        alpha=0.02,
        legend=False
    )
    ax.text(
        0.05,
        0.95,  # Position: 5% from left, 95% from bottom
        stats_text,
        transform=ax.transAxes,
        fontsize=9,
        verticalalignment="top",
        bbox=dict(boxstyle="round,pad=0.5", facecolor="white", alpha=0.8),
    )
    # Perform linear regression: y = slope * x + intercept
    slope, intercept = np.polyfit(
        corr_data[column1], corr_data[column2], 1  # Degree 1 for linear fit
    )
    # Create the line coordinates
    x_fit = np.linspace(corr_data[column1].min(), corr_data[column1].max(), 100)
    y_fit = slope * x_fit + intercept

    # Plot the line
    ax.plot(
        x_fit,
        y_fit,
        color="black",
        linestyle="-",
        linewidth=1.5,
        label=f"Linear Trend (R={corr_r:.2f})"
    )
#    ax.set( yscale="log")

plt.tight_layout(rect=[0, 0, 0.85, 1])  # Make space for the legend

# Save to SVG
save_path = f"./plots/{column2}_vs_{column1}_repeats_by_class.png"
figure.savefig(save_path, format="png", transparent=True)
print(f"Plot saved to {save_path}")

plt.show()
