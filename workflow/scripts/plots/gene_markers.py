import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from scipy.stats import zscore
import logging
import argparse

_log = logging.getLogger("gene_markers")

def set_logging(log_file, log_level):
    _log.setLevel(log_level)
    # create file handler that logs debug and higher level messages
    fh = logging.FileHandler(log_file)
    # create console handler with a higher log level
    ch = logging.StreamHandler()
    # create formatter and add it to the handlers
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    fh.setFormatter(formatter)
    # add the handlers to logger
    _log.addHandler(ch)
    _log.addHandler(fh)

def get_gene_markers(adata, celltype, marker_genes):

    gene_markers = adata[adata.obs["celltype"] == celltype, marker_genes['gene_symbol']].to_df()
    gene_markers.columns = marker_genes['gene_name']
    gene_markers.index = adata[adata.obs["celltype"] == celltype].obs['treatment']
    
    gene_markers = gene_markers.apply(zscore, axis=0)

    return gene_markers

def plot_gene_markers(markers_values, celltype, ax):
    
    sns.heatmap(markers_values.T, cmap='coolwarm', center=0, ax=ax, cbar_kws={'label': 'z-score'})

    ax.set_title(celltype, fontsize=16)
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=12)
    ax.set_xlabel('Gene', fontsize=14)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=12)
    ax.set_ylabel('Treatment', fontsize=14)

    return ax

def main():

    parser = argparse.ArgumentParser(description='Plot gene markers')
    parser.add_argument('--input-h5ad', type=str, help='Path to the AnnData object', required=True)
    parser.add_argument('--gene-markers-csv', type=str, help='Path to the gene markers', required=True)
    parser.add_argument('--output-plot', type=str, help='Path to save plot', required=True)
    parser.add_argument('--log', type=str, help='Path to log file', required=True)
    parser.add_argument('--log-level', type=str, help='Log level', default='INFO')
    args = parser.parse_args()

    set_logging(args.log, args.log_level)
    _log.debug(f"Command line arguments: {args}")

    adata = sc.read(args.input_h5ad)
    _log.info(f"Read AnnData object with shape {adata.X.shape}")
    gene_markers = pd.read_csv(args.gene_markers_csv, index_col=0)
    _log.info(f"Read gene markers with shape {gene_markers.shape}")
    _log.debug(f"Gene markers: {gene_markers}")


    fig, axs = plt.subplots(7, 2, figsize=(15, 20))

    for i, celltype in enumerate(adata.obs["celltype"].unique()):
        markers_values = get_gene_markers(adata, celltype, gene_markers)
        ax = axs[i//2, i%2]
        ax = plot_gene_markers(markers_values, celltype, ax)

    plt.tight_layout()
    fig.savefig(args.output_plot, dpi=300)

if __name__ == '__main__':
    main()