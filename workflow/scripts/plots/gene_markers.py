import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from scipy.stats import zscore
import logging
import argparse

from utils import set_logging, get_marker_gene_values

_log = logging.getLogger("gene_markers")

def plot_gene_markers(markers_values, celltype, ax):
    
    sns.heatmap(markers_values.T, cmap='coolwarm', center=0, ax=ax, cbar_kws={'label': 'z-score'})

    ax.set_title(celltype, fontsize=16)
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=12)
    ax.set_xlabel('Treatment', fontsize=14)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=12)
    ax.set_ylabel('Gene', fontsize=14)

    return ax

def select_top_markers(gene_markers, top_n=10):
    if 'coef' in gene_markers.columns:
        # Sort markers by absolute value of coefficient in descending order
        sorted_markers = gene_markers.reindex(gene_markers['coef'].abs().sort_values(ascending=False).index)
        # Select top N markers
        top_markers = sorted_markers.iloc[:top_n]
        return top_markers
    else:
        # If 'coef' column is not present, return all markers
        _log.info("No 'coef' column found in gene markers. Returning all markers.")
        return gene_markers

def main():

    parser = argparse.ArgumentParser(description='Plot gene markers')
    parser.add_argument('--input-h5ad', type=str, help='Path to the AnnData object', required=True)
    parser.add_argument('--gene-markers-csv', type=str, help='Path to the gene markers', required=True)
    parser.add_argument('--output-plot', type=str, help='Path to save plot', required=True)
    parser.add_argument('--log', type=str, help='Path to log file', required=True)
    parser.add_argument('--log-level', type=str, help='Log level', default='INFO')
    args = parser.parse_args()

    set_logging(_log, args.log, args.log_level)
    _log.debug(f"Command line arguments: {args}")

    adata = sc.read(args.input_h5ad)
    _log.info(f"Read AnnData object with shape {adata.X.shape}")
    gene_markers = pd.read_csv(args.gene_markers_csv, index_col=0)
    _log.info(f"Read gene markers with shape {gene_markers.shape}")
    _log.debug(f"Gene markers: {gene_markers}")

    # Select top 10 markers
    top_gene_markers = select_top_markers(gene_markers)
    _log.info(f"Selected top {len(top_gene_markers)} gene markers")
    _log.debug(f"Top gene markers: {top_gene_markers}")

    fig, axs = plt.subplots(7, 2, figsize=(15, 20))

    for i, celltype in enumerate(adata.obs["celltype"].unique()):

        markers_values = get_marker_gene_values(adata, celltype, top_gene_markers)
        markers_values = markers_values.apply(zscore, axis=0)
        ax = axs[i//2, i%2]
        ax = plot_gene_markers(markers_values, celltype, ax)

    plt.tight_layout()
    fig.savefig(args.output_plot, dpi=300)

if __name__ == '__main__':
    main()