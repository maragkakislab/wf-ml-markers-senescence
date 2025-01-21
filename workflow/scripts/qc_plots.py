import numpy as np
import scanpy as sc
import seaborn as sns
from scipy.stats import median_abs_deviation
import argparse
import logging
import matplotlib.pyplot as plt
import os

from utils import set_logging, get_default_parser

def calculate_qc_metrics(adata, gene_name_col='Gene_name', species='mouse'):

    logging.info(f"Using species {species}")

    if species == 'mouse':
        # mitochondrial genes
        adata.var["mt"] = adata.var[gene_name_col].str.startswith("mt-")
        # ribosomal genes
        adata.var["ribo"] = adata.var[gene_name_col].str.startswith(("Rps", "Rpl")) | adata.var[gene_name_col].str.contains("rRNA")
        # hemoglobin genes
        adata.var["hb"] = adata.var[gene_name_col].str.contains("^Hba|Hbb")

    elif species == 'human':
        # mitochondrial genes
        adata.var["mt"] = adata.var[gene_name_col].str.startswith("MT-")
        # ribosomal genes
        adata.var["ribo"] = adata.var[gene_name_col].str.startswith(("RPS", "RPL")) | adata.var[gene_name_col].str.contains("rRNA")
        # hemoglobin genes
        adata.var["hb"] = adata.var[gene_name_col].str.contains("^HB[^(P)]")

    else:
        raise ValueError(f"Species {species} not supported")
    
    # convert NaN to False
    adata.var["mt"] = adata.var["mt"].fillna(False)
    adata.var["ribo"] = adata.var["ribo"].fillna(False)
    adata.var["hb"] = adata.var["hb"].fillna(False)

    # convert metrics columns to bool
    adata.var["mt"] = adata.var["mt"].astype(bool)
    adata.var["ribo"] = adata.var["ribo"].astype(bool)
    adata.var["hb"] = adata.var["hb"].astype(bool)

    # calculate QC metrics
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
    )

    logging.debug(f"mitochondrial genes: {adata.var['mt'].sum()}")
    logging.debug(f"ribosomal genes: {adata.var['ribo'].sum()}")
    logging.debug(f"hemoglobin genes: {adata.var['hb'].sum()}")

    return adata

def plot_qc_metrics(adata, library_size_path, mito_path, ribo_path, hb_path):
    
    # library size plot
    fig, ax = plt.subplots(1, 1, figsize=(5, 3), dpi=300)
    sc.pl.violin(
        adata,
        ["total_counts"],
        groupby="celltype",
        size=3,
        rotation=90,
        show=False,
        ax=ax
    )
    ax.set_ylabel("Library size", fontsize=16)
    ax.set_xlabel("Cell type", fontsize=16)
    plt.savefig(library_size_path, bbox_inches="tight", dpi=300, transparent=True, format="pdf")
    logging.info(f"Saving library size plot to {library_size_path}")

    # mitochondrial gene plot
    if mito_path is not None:
        fig, ax = plt.subplots(1, 1, figsize=(5, 3), dpi=300)
        sc.pl.violin(adata, ["pct_counts_mt"], groupby="celltype", size=3, rotation=90, show=False, ax=ax)
        ax.set_ylabel("% mitochondrial genes", fontsize=16)
        ax.set_xlabel("Cell type", fontsize=16)
        plt.savefig(mito_path, bbox_inches="tight", dpi=300, transparent=True, format="pdf")
        logging.info(f"Saving mitochondrial gene plot to {mito_path}")


    # ribosomal gene plot
    if ribo_path is not None:
        fig, ax = plt.subplots(1, 1, figsize=(5, 3), dpi=300)
        sc.pl.violin(adata, ["pct_counts_ribo"], groupby="celltype", size=3, rotation=90, show=False, ax=ax)
        ax.set_ylabel("% ribosomal genes", fontsize=16)
        ax.set_xlabel("Cell type", fontsize=16)
        plt.savefig(ribo_path, bbox_inches="tight", dpi=300, transparent=True, format="pdf")
        logging.info(f"Saving ribosomal gene plot to {ribo_path}")

    # hemoglobin gene plot
    if hb_path is not None:
        fig, ax = plt.subplots(1, 1, figsize=(5, 3), dpi=300)
        sc.pl.violin(adata, ["pct_counts_hb"], groupby="celltype", size=3, rotation=90, show=False, ax=ax)
        ax.set_ylabel("% hemoglobin genes", fontsize=16)
        ax.set_xlabel("Cell type", fontsize=16)
        plt.savefig(hb_path, bbox_inches="tight", dpi=300, transparent=True, format="pdf")
        logging.info(f"Saving ribosomal gene plot to {ribo_path}")

def main():

    parser = get_default_parser(__doc__)
    parser.add_argument('--library_size_path', type=str, help='Path to save library size plot')
    parser.add_argument('--mito_path', type=str, help='Path to save mitochondrial gene plot', default=None)
    parser.add_argument('--ribo_path', type=str, help='Path to save ribosomal gene plot', default=None)
    parser.add_argument('--hb_path', type=str, help='Path to save hemoglobin gene plot', default=None)
    parser.add_argument('input', type=str, help='Input anndata file')
    args = parser.parse_args()

    set_logging(args.log, args.log_level)
    logging.debug(f"Command line arguments: {args}")

    # read anndata object
    adata = sc.read(args.input)
    logging.info(f"Read AnnData object with shape {adata.X.shape}")

    # calculate QC metrics
    adata = calculate_qc_metrics(adata, species="human")

    # plot QC metrics
    plot_qc_metrics(adata, args.library_size_path, args.mito_path, args.ribo_path, args.hb_path)

if __name__ == '__main__':
    main()