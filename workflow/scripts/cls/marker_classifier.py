import argparse
import logging

import numpy as np
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

    # Check which markers are available in the data
    data_features = set(adata.var_names)
    available_markers = [
        marker for marker in markers.index if marker in data_features
    ]
    missing_markers = [
        marker for marker in markers.index if marker not in data_features
    ]

    _log.info(f"Number of markers provided: {len(markers)}")
    _log.info(f"Number of markers available in the data: {len(available_markers)}")
    _log.info(f"Number of markers missing from the data: {len(missing_markers)}")

    if available_markers == []:
        raise ValueError(
            "No markers are available in the data. Unable to classify samples."
        )

    # Use only available markers and keep marker/coefficient order aligned.
    coefs = markers.loc[available_markers, 'coef'].values
    X = adata[:, available_markers].X
    preds = np.asarray(X @ coefs).reshape(-1)

    result = pd.DataFrame({'score': preds}, index=adata.obs_names.copy())
    if 'is_sen' in adata.obs:
        result.insert(0, 'label', adata.obs['is_sen'].to_numpy(copy=True))
    return result

def main(argv=None):
    parser = argparse.ArgumentParser(description='Classify samples using marker genes')
    parser.add_argument('--input-h5ad', type=str, help='Path to the AnnData object', required=True)
    parser.add_argument('--markers', type=str, help='Path to the marker genes', required=True)
    parser.add_argument('--output-results-csv', type=str, help='Path to save results', required=True)
    parser.add_argument('--log', type=str, help='Path to log file', required=True)
    parser.add_argument('--log-level', type=str, help='Log level', default='INFO')
    args = parser.parse_args(argv)

    set_logging(args.log, args.log_level)
    logging.debug(f"Command line arguments: {args}")

    adata = sc.read(args.input_h5ad)
    logging.info(f"Read AnnData object with shape {adata.X.shape}")
    
    _log.info(f"Do log1p transformation")
    sc.pp.log1p(adata)

    markers = pd.read_csv(args.markers, index_col=0)

    try:
        result = classify_samples(adata, markers)
    except ValueError as error:
        _log.error(str(error))
        parser.exit(status=1, message=f"error: {error}\n")

    result.to_csv(args.output_results_csv)

if __name__ == '__main__':
    main()
