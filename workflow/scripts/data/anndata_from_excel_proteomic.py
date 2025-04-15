import pandas as pd
import scanpy as sc
import argparse
import logging 

from utils import set_logging

_log = logging.getLogger("anndata_from_excel_proteomic")

def main():
    parser = argparse.ArgumentParser(description='Convert an Excel file to an AnnData object.')
    parser.add_argument('--input-excel', type=str, required=True, help='Input Excel file with counts.')
    parser.add_argument('--sheet', type=str, required=True, help='Sheet name in the Excel file.')
    parser.add_argument('--output-h5ad', type=str, required=True, help='Output h5ad file.')
    parser.add_argument('--id_column', type=str, required=True, help='Column name with protein ID.')
    parser.add_argument('--counts_columns', type=str, nargs='+', required=True, help='Column names with counts.')
    parser.add_argument('--log', type=str, required=True, help='Log file.')
    parser.add_argument('--log-level', type=str, required=True, help='Log level.')

    args = parser.parse_args()
    set_logging(_log, args.log, args.log_level)

    df = pd.read_excel(args.input_excel, sheet_name=args.sheet, index_col=args.id_column)
    _log.info(f"Read {df.shape[1]} rows and {df.shape[0]} columns from {args.input_excel}")

    counts = df[args.counts_columns]
    _log.info(f"Extracted counts for {counts.shape[1]} samples")

    # convert float to int by multiplying by 10e14
    counts = (counts * 10e14).astype(int)
    _log.info(f"Converted counts to integers")

    obs = pd.DataFrame(index=counts.columns)
    obs['sample'] = obs.index
    obs['is_sen'] = obs['sample'].apply(lambda x: 1 if x.startswith('Etop') else 0)
    _log.info(f"Extracted sample and is_sen columns from obs")

    adata = sc.AnnData(X=counts.T, obs=obs)
    _log.info(f"Created AnnData object with shape {adata.X.shape}")

    adata.var_names_make_unique()
    _log.info(f"Made variable names unique")

    _log.info(f"Saving AnnData object to {args.output_h5ad}")
    adata.write(args.output_h5ad)

if __name__ == '__main__':
    main()
