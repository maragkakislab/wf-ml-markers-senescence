import pandas as pd
import scanpy as sc
import argparse
import logging 

from utils import set_logging

_log = logging.getLogger("anndata_from_txt")

class dotdict(dict):
    """dot.notation access to dictionary attributes"""
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

def assemble_anndata(input_file, output_h5ad, mapping_tsv, metadata_path,
                     input_gene_col='ENSG', mapping_gene_col='Gene stable ID', 
                     mapping_name_col='Gene name', output_txt=None):
    
    counts = pd.read_csv(input_file, sep = '\t', index_col=0).T
    _log.info(f"Read {counts.shape[1]} genes and {counts.shape[0]} samples")
    adata = sc.AnnData(counts)

    adata.obs['sample'] = adata.obs_names 

    metadata = pd.read_excel(metadata_path)
    # add metadata to observations
    adata.obs = adata.obs.reset_index().merge(metadata, how="left").set_index('index')
    _log.info(f"Added metadata to observations. {adata.obs.columns}")
    adata.obs['is_sen'] = adata.obs.apply(lambda x: 1 if str(x['condition']).startswith('sen') else 0, axis=1)
    _log.info(f"Added is_sen column to observations. {adata.obs['is_sen'].sum()} senescent samples.")

    # Load the mapping file and create a dictionary for fast lookup
    mapping_df = pd.read_csv(mapping_tsv, sep='\t')
    mapping_dict = dict(zip(mapping_df[mapping_gene_col], mapping_df[mapping_name_col]))

    adata.var_names_make_unique()
    adata.var[input_gene_col] = adata.var.index
    # add gene names to variables
    adata.var['Gene_name'] = adata.var[input_gene_col].map(mapping_dict)
    _log.info(f"Added gene names to variables. NaNs: {adata.var['Gene_name'].isna().sum()}")
    # fill missing gene names with gene IDs
    adata.var['Gene_name'] = adata.var['Gene_name'].fillna(adata.var[input_gene_col])
    _log.info(f"Filled missing gene names with gene IDs. NaNs: {adata.var['Gene_name'].isna().sum()}")
    # set gene symbol as index
    adata.var.set_index(input_gene_col, inplace=True)

    _log.info(f"Saving AnnData object to {output_h5ad}")
    adata.write(output_h5ad)

    if output_txt is not None:
        counts = counts.T
        counts['Gene_name'] = counts.index.map(mapping_dict)
        counts.index = counts['Gene_name']
        counts = counts.drop(columns=['Gene_name'])
        # drop rows with NaN gene names
        counts = counts[~counts.index.isna()]
        _log.info(f"Kept {counts.shape[0]} genes with non-NaN gene names.")
        counts.to_csv(output_txt, sep='\t')
        _log.info(f"Saved gene counts to {output_txt}")

    return adata

def main():

    if 'snakemake' not in globals():
        # Parse arguments if the script is run from the command line
        parser = argparse.ArgumentParser(description='Create anndata object from a count file.')
        parser.add_argument('--input-tsv', required=True, help='Path to the input TSV file with counts.')
        parser.add_argument('--mapping-tsv', required=True, help='Path to the TSV file with columns for mapping gene IDs to gene names.')
        parser.add_argument('--metadata-excel', required=True, help='Path to the metadata Excel file.')
        parser.add_argument('--output-h5ad', required=True, help='Path to the output H5AD file.')
        parser.add_argument('--output-txt', default=None, help='Path to the output TSV file with gene names. Used as input for SenCID.')
        parser.add_argument('--input-gene-col', default='ENSG', help='Column name in the input TSV containing gene IDs.')
        parser.add_argument('--mapping-gene-col', default='Gene stable ID', help='Column name in the mapping TSV containing gene IDs.')
        parser.add_argument('--mapping-name-col', default='Gene name', help='Column name in the mapping TSV containing gene names.')
        parser.add_argument('--log', default="log.log", help='Path to log file.')
        parser.add_argument('--log-level', default="INFO", help='Log level.')
        args = parser.parse_args()
    else:
        # Arguments will be passed from the Snakefile
        args = {}
        args['input_tsv'] = snakemake.input['counts']
        args['mapping_tsv'] = snakemake.input['mapping']
        args['metadata_excel'] = snakemake.input['metadata']
        args['output_h5ad'] = snakemake.output['h5ad']
        args['output_txt'] = snakemake.output.get('txt')
        args['input_gene_col'] = snakemake.params.get('input_gene_col', 'ENSG')
        args['mapping_gene_col'] = snakemake.params.get('mapping_gene_col', 'Gene stable ID')
        args['mapping_name_col'] = snakemake.params.get('mapping_name_col', 'Gene name')
        args = dotdict(args)

    set_logging(_log, args.log, args.log_level)
    _log.debug(f"Command line arguments: {args}")

    assemble_anndata(
        input_file=args.input_tsv,
        mapping_tsv=args.mapping_tsv,
        metadata_path=args.metadata_excel,
        output_h5ad=args.output_h5ad,
        output_txt=args.output_txt,
        input_gene_col=args.input_gene_col,
        mapping_gene_col=args.mapping_gene_col,
        mapping_name_col=args.mapping_name_col,
    )

if __name__ == "__main__":
    main()