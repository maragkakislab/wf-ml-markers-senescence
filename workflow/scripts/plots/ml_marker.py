import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import logging
import argparse

from utils import set_logging

_log = logging.getLogger("ml_marker_plotter")

def compute_ml_marker(adata, ml_marker):
    sc.pp.log1p(adata)
    
    # Check if 'sample' column exists, if not, create it
    if 'sample' not in adata.obs.columns:
        if 'treatment' in adata.obs.columns and 'replicate' in adata.obs.columns:
            adata.obs['sample'] = adata.obs.apply(lambda x: f"{x['treatment']}_{x['replicate']}", axis=1)
        else:
            adata.obs['sample'] = adata.obs_names
        _log.info("Created 'sample' column in adata.obs")

    # Check which markers are available in the data
    available_markers = list(set(ml_marker.index) & set(adata.var_names))
    missing_markers = list(set(ml_marker.index) - set(adata.var_names))

    _log.info(f"Number of markers provided: {len(ml_marker)}")
    _log.info(f"Number of markers available in the data: {len(available_markers)}")
    _log.info(f"Number of markers missing from the data: {len(missing_markers)}")
    if missing_markers:
        _log.warning(f"Missing markers: {', '.join(missing_markers)}")

    # Use only available markers
    ml_marker = ml_marker.loc[available_markers]

    # Check if 'celltype' column exists, if not, use a default value
    if 'celltype' not in adata.obs.columns:
        adata.obs['celltype'] = 'default_celltype'
        _log.warning("No 'celltype' column found in adata.obs. Using 'default_celltype'.")

    ml_markers = pd.DataFrame(index=adata.obs['sample'].unique(), columns=adata.obs["celltype"].unique())
    for i, celltype in enumerate(adata.obs["celltype"].unique()):

        coef_markers = adata[adata.obs["celltype"] == celltype, ml_marker.index.to_list()].to_df()
        coef_markers.index = adata[adata.obs["celltype"] == celltype].obs['sample']

        ml = coef_markers.apply(lambda x: np.sum(np.multiply(x, ml_marker['coef'].values)), axis=1)

        for m in ml.index:
            ml_markers.loc[m, celltype] = ml[m]
            

    ml_markers = ml_markers.astype(float)
    return ml_markers

def main():

    parser = argparse.ArgumentParser(description='Plot ML marker')
    parser.add_argument('--input-h5ad', type=str, help='Path to the AnnData object', required=True)
    parser.add_argument('--ml-marker-csv', type=str, help='Path to the machine learning marker', required=True)
    parser.add_argument('--output-plot', type=str, help='Path to save plot', required=True)
    parser.add_argument('--log', type=str, help='Path to log file', required=True)
    parser.add_argument('--log-level', type=str, help='Log level', default='INFO')
    args = parser.parse_args()

    set_logging(_log, args.log, args.log_level)
    _log.debug(f"Command line arguments: {args}")

    try:
        adata = sc.read(args.input_h5ad)
        _log.info(f"Read AnnData object with shape {adata.X.shape}")
    except Exception as e:
        _log.error(f"Failed to read AnnData object: {e}")
        raise

    try:
        ml_markers = pd.read_csv(args.ml_marker_csv, index_col=0)
        _log.info(f"Read machine learning marker with shape {ml_markers.shape}")
        _log.debug(f"Machine learning marker: {ml_markers}")
    except Exception as e:
        _log.error(f"Failed to read ML marker CSV: {e}")
        raise

    fig, ax = plt.subplots(1, 1, figsize=(15, 5))
    try:
        computed_ml_marker = compute_ml_marker(adata, ml_markers)
        _log.info(f"Computed ML marker with shape {computed_ml_marker.shape}")
        _log.debug(f"ML marker: {computed_ml_marker}")
    except Exception as e:
        _log.error(f"Failed to compute ML marker: {e}")
        raise

    try:
        sns.heatmap(computed_ml_marker.T, yticklabels=computed_ml_marker.columns, xticklabels=computed_ml_marker.index, cmap='coolwarm', center=0, cbar_kws={'label': 'ML marker'}, ax=ax)
        ax.set_xlabel('Sample', fontsize=14)
        ax.set_ylabel('Cell type', fontsize=14)
        ax.set_title('ML marker', fontsize=16)
        plt.tight_layout()
        fig.savefig(args.output_plot, dpi=300)
        _log.info(f"Saved plot to {args.output_plot}")
    except Exception as e:
        _log.error(f"Failed to plot ML marker: {e}")
        raise

if __name__ == '__main__':
    main()