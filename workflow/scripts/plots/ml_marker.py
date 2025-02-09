import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import logging
import argparse

from utils import set_logging, get_marker_gene_values

_log = logging.getLogger("gene_markers_with_ml_marker")

def compute_ml_marker(adata, ml_marker):
    sc.pp.log1p(adata)
    adata.obs['sample'] = adata.obs.apply(lambda x: f"{x['treatment']}_{x['replicate']}", axis=1)
    ml_markers = pd.DataFrame(index=adata.obs['sample'].unique(), columns=adata.obs["celltype"].unique())
    for i, celltype in enumerate(adata.obs["celltype"].unique()):

        coef_markers =  adata[adata.obs["celltype"] == celltype, ml_marker.index.to_list()].to_df()
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

    adata = sc.read(args.input_h5ad)
    _log.info(f"Read AnnData object with shape {adata.X.shape}")
    ml_markers = pd.read_csv(args.ml_marker_csv, index_col=0)
    _log.info(f"Read machine learning marker with shape {ml_markers.shape}")
    _log.debug(f"Machine learning marker: {ml_markers}")

    fig, ax = plt.subplots(1, 1, figsize=(15, 5))
    computed_ml_marker = compute_ml_marker(adata, ml_markers)
    _log.info(f"Computed ML marker with shape {computed_ml_marker.shape}")
    _log.debug(f"ML marker: {computed_ml_marker}")
    sns.heatmap(computed_ml_marker.T,  yticklabels=computed_ml_marker.columns, xticklabels=computed_ml_marker.index, cmap='coolwarm', center=0, cbar_kws={'label': 'ML marker'}, ax=ax)
    ax.set_xlabel('Sample', fontsize=14)
    ax.set_ylabel('Cell type', fontsize=14)
    ax.set_title('ML marker', fontsize=16)
    plt.tight_layout()
    fig.savefig(args.output_plot, dpi=300)

if __name__ == '__main__':
    main()