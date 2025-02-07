import argparse
import logging
import scanpy as sc
import pandas as pd

_log = logging.getLogger("marker_classifier")

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

def classify_samples(adata, markers):

    features = markers.index
    coefs = markers['coef'].values

    X = adata[:, features].X
    y = adata.obs['is_sen']
    preds = X @ coefs
    return pd.DataFrame({'label': y, 'score': preds})

def main():
    parser = argparse.ArgumentParser(description='Classify samples using marker genes')
    parser.add_argument('--input-h5ad', type=str, help='Path to the AnnData object', required=True)
    parser.add_argument('--markers', type=str, help='Path to the marker genes', required=True)
    parser.add_argument('--output-results-csv', type=str, help='Path to save results', required=True)
    parser.add_argument('--log', type=str, help='Path to log file', required=True)
    parser.add_argument('--log-level', type=str, help='Log level', default='INFO')
    args = parser.parse_args()

    set_logging(args.log, args.log_level)
    logging.debug(f"Command line arguments: {args}")

    adata = sc.read(args.input_h5ad)
    logging.info(f"Read AnnData object with shape {adata.X.shape}")
    
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    markers = pd.read_csv(args.markers, index_col=0)

    result = classify_samples(adata, markers)
    result.to_csv(args.output_results_csv)

if __name__ == '__main__':
    main()


