import pandas as pd
import scanpy as sc
import argparse
import logging 

from utils import set_logging

_log = logging.getLogger("anndata_from_excel")

def main():
    parser = argparse.ArgumentParser(description='Convert an Excel file to an AnnData object.')
    parser.add_argument('--input-excel', type=str, required=True, help='Input Excel file with counts.')
    parser.add_argument('--output-h5ad', type=str, required=True, help='Output h5ad file.')
    parser.add_argument('--output-txt', type=str, default=None, help='Output txt file with counts.')
    parser.add_argument('--var-columns', type=int, required=True, help='Number of columns with variables.')
    parser.add_argument('--index-column', type=str, required=True, help='Column name with count IDs.')
    parser.add_argument('--gene-name-column', type=str, required=True, help='Column name with gene names.')
    parser.add_argument('--log', type=str, required=True, help='Log file.')
    parser.add_argument('--log-level', type=str, required=True, help='Log level.')

    args = parser.parse_args()
    set_logging(_log, args.log, args.log_level)

    df = pd.read_excel(args.input_excel)
    _log.info(f"Read {df.shape[1]} rows and {df.shape[0]} columns from {args.input_excel}")
    df.fillna(0, inplace=True)
    # drop rows where index columns is NaN
    df = df.dropna(subset=[args.index_column])
    _log.info(f"Kept {df.shape[0]} rows with non-NaN {args.index_column} values")

    var = df.iloc[:, :args.var_columns]
    _log.info(f"Extracted {var.shape[1]} variables: {var.columns}")
    var.set_index(args.index_column, inplace=True)
    _log.info(f"Set {args.index_column} as index for variables")

    counts = df.iloc[:, args.var_columns:]
    _log.info(f"Extracted counts for {counts.shape[1]} samples")
    counts.index = df[args.index_column]
    _log.info(f"Set {args.index_column} as index for counts")


    obs = pd.DataFrame(df.columns.values[args.var_columns:], columns=['obs'])
    _log.info(f"Extracted {obs.shape[0]} observations")
    obs['celltype'] = obs['obs'].str.split('_').str[0]
    obs['treatment'] = obs['obs'].str.split('_').str[1]
    obs['replicate'] = obs['obs'].str.split('_').str[2]
    obs['is_sen'] = obs.apply(lambda x: 0 if x['treatment'] in ['P', 'EV'] else 1, axis=1)
    _log.info(f"Extracted celltype, treatment and replicate from obs")
    obs.set_index('obs', inplace=True)
    _log.info(f"Set obs as index for observations")

    var.rename(columns={args.gene_name_column: 'Gene_name'}, inplace=True)
    _log.info(f"Renamed {args.gene_name_column} to Gene_name")

    adata = sc.AnnData(X=counts.T, obs=obs, var=var)
    _log.info(f"Created AnnData object with shape {adata.X.shape}")

    adata.var_names_make_unique()
    _log.info(f"Made variable names unique")

    # convert each column in var to string
    for col in adata.var.columns:
        adata.var[col] = adata.var[col].astype(str)
        _log.info(f"Converted {col} to string")

    # check if some column in var contains list 
    for col in adata.var.columns:
        if adata.var[col].apply(lambda x: isinstance(x, list)).any():
            adata.var[col] = adata.var[col].tolist()
            _log.info(f"Converted {col} to list")

    _log.info(f"Saving AnnData object to {args.output_h5ad}")
    adata.write(args.output_h5ad)

    if args.output_txt is not None:
        adata.var.set_index('Gene_name', inplace=True)
        adata.var.index = adata.var.index.astype(str)
        adata.var = adata.var[~adata.var.index.isna()]
        adata.var_names_make_unique()

        adata_txt = adata.to_df().T
        adata_txt.columns = adata.obs.index
        adata_txt.to_csv(args.output_txt, sep='\t')
        _log.info(f"Saved counts to {args.output_txt}")

if __name__ == '__main__':
    main()