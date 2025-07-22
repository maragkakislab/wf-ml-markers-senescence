#!/usr/bin/env python3

"""
Assembles AnnData object from a counts files from the SenCat project. 

The counts files are expected to have the following format:
- Each file contains counts for one cell type, including all treatments and replicates.
- Files are tab-separated.
- First column is the gene symbol.
- Columns containing counts have names ending with 'count' or 'counts'.
- Column names contain information about the cell type, treatment and replicate in the following format:
    *_{celltype}_{treatment}_{replicate}.count[s]
    regex: .*_(\w+)_(\w+)_(\d+).count[s]*
- Files additionally contain columns 'Gene_name' and 'Gene_Length'.

The script does the following:
- Reads all the input files.
- Concatenates the counts from all files and fills missing values with 0.
- Creates following observations for AnnData:
    * index: sample name in the format {celltype}_{treatment}_{replicate}.
    * `celltype`: cell type of the sample.
    * `treatment`: treatment of the sample.
    * `replicate`: replicate number.
    * `label`: 0 if the sample is control, 1 if the sample is treatment.
- Treatments `P` and `EV` are considered as control, all other treatments are considered as treatment.
- Creates following variables for AnnData:
    * index: gene symbol.
    * `Gene_name`: gene name.
    * `Gene_Length`: gene length.
- Writes the AnnData object to the output h5ad file.

Additionally, if `--output_txt` is provided, the script writes the counts to a tab-separated file with gene names as index.
Only genes with non-NaN gene names are kept. This output is used for a downstream analysis with SenCID tool.
"""

import pandas as pd
import anndata as an
import logging
import os
from scipy.sparse import csr_matrix
import argparse

from utils import set_logging

_log = logging.getLogger("assemble_anndata")

def assemble_anndata(input_files, output_path, output_txt=None):

    # read all files into a dictionary
    data = {}
    var = []
    for file in input_files:

        _log.info(f"Reading {file}")

        if not os.path.isfile(file):
            raise ValueError(f"File {file} not found")
        if not file.endswith('.txt'):
            raise ValueError(f"File {file} is not a txt file")

        df = pd.read_csv(file, sep='\t', index_col=0, low_memory=False)
        _log.debug(f"Read {df.shape[0]} genes and {df.shape[1]} columns")

        # remove white spaces from column names
        df.columns = df.columns.str.strip()
        # remove white spaces from gene names
        df.index = df.index.str.strip()

        # get columns with raw counts
        counts_cols = [col for col in df.columns if (col.split('.')[-1] in ['count', 'counts'])]
        if len(counts_cols) == 0:
            raise ValueError(f"No count columns found in {file}")
        _log.debug(f"Found {len(counts_cols)} count columns: {counts_cols}")
        data[file] = df[counts_cols]

        var.append(df[['Gene_name', 'Gene_Length']])

    # get union of all indices
    all_indices = set()
    for k, v in data.items():
        all_indices = all_indices.union(set(v.index))
    _log.info(f"Found {len(all_indices)} unique genes in all files combined.")

    # fill missing counts with 0
    for k, v in data.items():
        data[k] = v.reindex(all_indices).fillna(0)

    # concatenate all dataframes with counts
    df = pd.concat(data.values(), axis=1)
    _log.info(f"Concatenated counts from {len(data)} files. Final shape: {df.shape}.")
    
    def obs_from_col_name(col_name):
        parts = col_name.split('_')
        # indexing from the back, because some samples have additional prefixes
        celltype = parts[-3]
        treatment = parts[-2]
        replicate = parts[-1].split('.')[0]

        return {
            'celltype': celltype,
            'treatment': treatment,
            'replicate': replicate
        }

    # create observations from column names
    obs = pd.DataFrame([obs_from_col_name(col) for col in df.columns])
    # create label column based on treatment - P and EV are control, all others are treatment
    obs['is_sen'] = obs.apply(lambda x: 0 if x['treatment'] in ['P', 'EV'] else 1, axis=1)
    _log.debug(f"Created observations with {obs.shape[0]} samples.")
    _log.debug(f"Control samples: {obs[obs['is_sen'] == 0].shape[0]}. Treatment samples: {obs[obs['is_sen'] == 1].shape[0]}.")
    _log.debug(f"Cell types: {obs['celltype'].unique()}. Treatments: {obs['treatment'].unique()}. Replicates: {obs['replicate'].unique()}.")

    # create sample name from cell type, treatment and replicate
    obs.index = obs.apply(lambda x: f"{x['celltype']}_{x['treatment']}_{x['replicate']}", axis=1)

    var = pd.concat(var)
    var = var[~var.index.duplicated(keep='first')]
    _log.debug(f"Created variables with {var.shape[0]} genes.")

    # reorder rows in df to match var
    df = df.reindex(var.index)

    # assemble AnnData object from counts, metadata about samples and genes
    adata = an.AnnData(X=df.values.T, obs=obs, var=var)
    _log.debug(f"Created AnnData object with shape {adata.X.shape}.")

    # make adata.X sparse
    adata.X = csr_matrix(adata.X)
    _log.debug(f"Converted counts to sparse matrix.")

    adata.write(output_path)
    _log.info(f"Saved AnnData object to {output_path}")

    # write txt output if requested
    if output_txt is not None:
        # reindex to gene names using mapping from var
        df.index = var.loc[df.index, 'Gene_name']
        # drop nan indices
        df = df[~df.index.isna()]
        _log.debug(f"Kept {df.shape[0]} genes with non-NaN gene names.")
        df.columns = obs.index
        df.to_csv(output_txt, sep='\t')
        _log.info(f"Saved gene counts to {output_txt}")

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input_files', type=str, nargs='+', help='Input counts files')
    parser.add_argument('--output_h5ad', type=str, required=True, help='Output h5ad file')
    parser.add_argument('--output_txt', type=str, default=None, help='Output txt file')
    parser.add_argument('--log', type=str, default="log.log", help='Log file')
    parser.add_argument('--log-level', type=str, default="INFO", help='Log level')
    args = parser.parse_args()

    set_logging(_log, args.log, args.log_level)
    _log.debug(f"Command line arguments: {args}")

    assemble_anndata(args.input_files, args.output_h5ad, args.output_txt)

if __name__ == '__main__':
    main()
