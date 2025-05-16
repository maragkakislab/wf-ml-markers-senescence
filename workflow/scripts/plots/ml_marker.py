import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import logging
import argparse

from utils import set_logging, get_cmap

_log = logging.getLogger("ml_marker_plotter")

def main():

    parser = argparse.ArgumentParser(description='Plot ML marker')
    parser.add_argument('--input-h5ad', type=str, help='Path to the AnnData object', required=True)
    parser.add_argument('--results-csv', type=str, help='Path to the classification results.', required=True)
    parser.add_argument('--output-plot', type=str, help='Path to save plot', required=True)
    parser.add_argument('--log', type=str, help='Path to log file', required=True)
    parser.add_argument('--log-level', type=str, help='Log level', default='INFO')
    args = parser.parse_args()

    set_logging(_log, args.log, args.log_level)
    _log.debug(f"Command line arguments: {args}")

    adata = sc.read(args.input_h5ad)
    _log.info(f"Read AnnData object with shape {adata.X.shape}")
    # drop samples where index contains Ribo
    adata = adata[~adata.obs_names.str.contains('Ribo')]

    if args.results_csv.endswith(('.tsv', '.txt')):
        results = pd.read_csv(args.results_csv, sep='\t', index_col=0)
    elif args.results_csv.endswith('.csv'):
        results = pd.read_csv(args.results_csv, index_col=0)
    else:
        _log.error('Input file must be a tsv or csv file')

    _log.info(f"Read classification results with a shape {results.shape}")
    _log.debug(f"Classification results: {results.head()}")
    
    # Check if 'sample' column exists, if not, create it
    if 'sample' not in adata.obs.columns:
        if 'treatment' in adata.obs.columns and 'replicate' in adata.obs.columns:
            adata.obs['sample'] = adata.obs.apply(lambda x: f"{x['treatment']}_{x['replicate']}", axis=1)
        else:
            adata.obs['sample'] = adata.obs_names
        _log.info("Created 'sample' column in adata.obs")

    # Check if 'celltype' column exists, if not, use a default value
    if 'celltype' not in adata.obs.columns:
        adata.obs['celltype'] = 'default_celltype'
        _log.warning("No 'celltype' column found in adata.obs. Using 'default_celltype'.")

    adata.obs = adata.obs.merge(results['score'], left_index=True, right_index=True)
    ml_marker = adata.obs.pivot_table(index='celltype', columns='sample', values='score').reindex(columns=adata.obs.index)

    fig, ax = plt.subplots(1, 1, figsize=(15, 5))
    _log.info(f"Loaded ML marker with shape {ml_marker.shape}")
    _log.debug(f"ML marker: {ml_marker}")

    cmap = get_cmap()
    sns.heatmap(ml_marker, yticklabels=ml_marker.index, xticklabels=ml_marker.columns, cmap=cmap, cbar_kws={'label': 'ML marker'}, ax=ax)
    ax.set_xlabel('Sample', fontsize=14)
    ax.set_ylabel('Cell type', fontsize=14)
    ax.set_title('ML marker', fontsize=16)
    plt.tight_layout()
    fig.savefig(args.output_plot, dpi=300)
    _log.info(f"Saved plot to {args.output_plot}")

if __name__ == '__main__':
    main()