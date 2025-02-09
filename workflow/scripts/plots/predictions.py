import argparse
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import precision_recall_curve, auc
import logging

from utils import set_logging

_log = logging.getLogger("predictions")

def plot_results(results, save_pdf, label_col='is_sen', score_col='score'):
    fig, axs = plt.subplots(1, 2, figsize=(7, 3), dpi=600)

    precision, recall, _ = precision_recall_curve(results[label_col], results[score_col])
    axs[0].plot(recall, precision, label=f'AUC: {round(auc(recall, precision), 3)}')
    axs[0].set_xlabel('Recall', fontsize=12)
    axs[0].set_ylabel('Precision', fontsize=12)
    # add random prediction score
    axs[0].plot([0, 1], [results[label_col].mean(), results[label_col].mean()], linestyle='--', color='gray', label='Random')
    axs[0].legend()

    axs[0].set_xlim([0.0, 1.1])
    axs[0].set_ylim([0.0, 1.1])

    sns.violinplot(x=label_col, y=score_col, data=results, cut=0, ax=axs[1])
    sns.stripplot(x=label_col, y=score_col, data=results, jitter=True, color='black', alpha=0.5, ax=axs[1])

    axs[1].set_xlabel('Senescent', fontsize=12)
    axs[1].set_ylabel('Prediction', fontsize=12)

    plt.tight_layout()
    plt.savefig(save_pdf, dpi=600)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('results_path', type=str)
    parser.add_argument('output_pdf', type=str)
    parser.add_argument('h5ad_path', type=str)
    parser.add_argument('--log', type=str, help='Path to log file', required=True)
    parser.add_argument('--log-level', type=str, help='Log level', default='INFO')
    args = parser.parse_args()
    
    set_logging(_log, args.log, args.log_level)
    _log.debug(f"Command line arguments: {args}")

    if args.results_path.endswith(('.tsv', '.txt')):
        results = pd.read_csv(args.results_path, sep='\t', index_col=0)
    elif args.results_path.endswith('.csv'):
        results = pd.read_csv(args.results_path, index_col=0)
    else:
        _log.error('Input file must be a tsv or csv file')
    

    adata = sc.read(args.h5ad_path)

    results = results.merge(adata.obs["is_sen"], left_index=True, right_index=True)
    # make sure the column is 0/1
    results["is_sen"] = results["is_sen"].astype(int)
    print(results.head())

    # drop rows where index contains Ribo
    results = results[~results.index.str.contains('Ribo')]

    plot_results(results, args.output_pdf)

if __name__ == '__main__':
    main()