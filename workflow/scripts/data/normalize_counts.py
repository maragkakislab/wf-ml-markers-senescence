from pydeseq2.dds import DeseqDataSet
import scanpy as sc
import logging
import argparse
from scipy.sparse import csr_matrix

_log = logging.getLogger("normalize_counts")

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

def main():
    parser = argparse.ArgumentParser(description='Normalize counts using DESeq2')
    parser.add_argument('--input-h5ad', type=str, help='Path to the AnnData object')
    parser.add_argument('--design', type=str, help='Design factors for DESeq2.')
    parser.add_argument('--output-h5ad', type=str, help='Path to save normnalized AnnData object', required=True)
    parser.add_argument('--log', type=str, help='Path to log file', required=True)
    parser.add_argument('--log-level', type=str, help='Log level', default='INFO')
    args = parser.parse_args()

    set_logging(args.log, args.log_level)
    logging.debug(f"Command line arguments: {args}")

    adata = sc.read(args.input_h5ad)
    logging.info(f"Read AnnData object with shape {adata.X.shape}")

    # if sparse matrix, convert to dense matrix
    if isinstance(adata.X, csr_matrix):
        adata.X = adata.X.todense()
        logging.info(f"Converted sparse matrix to dense matrix")

    logging.info(f"Converted sparse matrix to dense matrix")

    logging.info(f"Removing genes with counts for less than 5% of samples")
    adata = adata[:, adata.X.astype(bool).sum(0) > 0.05 * adata.shape[0]]
    logging.debug(f"Filtered genes. Shape: {adata.X.shape}")

    logging.info(f"Removing genes with total counts less than 380")
    adata = adata[:, adata.X.sum(0) > 380]
    logging.debug(f"Filtered genes. Shape: {adata.X.shape}")

    nonzero_genes = adata[:, adata.X.astype(bool).sum(0) == adata.shape[0]].shape[1]
    logging.info(f"Number of genes with counts in all samples: {nonzero_genes}")
    if nonzero_genes < 100:
        logging.warning(f"Number of genes with counts in all samples is less than 100. Adding count of 1 to all genes.")
        logging.info(f"Storing counts in a layer 'counts'")
        adata.layers['counts'] = adata.X
        adata.X = adata.X + 1

    logging.info(f"Creating DeseqDataSet object")
    dds = DeseqDataSet(adata=adata, design=args.design)
    logging.debug(f"Created DeseqDataSet object. Shape: {dds.X.shape}")

    logging.info(f"Estimating size factors")
    dds.fit_size_factors()
    logging.debug(f"{dds}")

    logging.info(f"Loading normalized counts to adata.X")
    dds.X = dds.layers['normed_counts']

    del dds.obsm['design_matrix']

    logging.info(f"Saving normalized counts to {args.output_h5ad}")
    dds.write_h5ad(args.output_h5ad)

if __name__ == '__main__':
    main()