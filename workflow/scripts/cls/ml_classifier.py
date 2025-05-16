from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import precision_recall_curve, auc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from collections import Counter
import logging
import argparse
import random

_log = logging.getLogger("ml_classifier")

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

def univariate_feature_selection(adata, n_features=500, label_col='is_sen'):
    _log.info(f"Selecting {n_features} features")
    selector = SelectKBest(f_classif, k=n_features)
    selector.fit(adata.X, adata.obs[label_col])
    features = selector.get_support(indices=True)
    return features

def train_test_split(adata, celltype, features, label_col='is_sen', celltype_col='celltype'):

    train = adata[adata.obs[celltype_col] != celltype, features]
    test = adata[adata.obs[celltype_col] == celltype, features]
    
    return train.X, train.obs[label_col], test.X, test.obs[label_col]

def train_models_for_celltypes(adata, features, label_col='is_sen', celltype_col='celltype', penalty=None):
    models = {}
    for celltype in adata.obs[celltype_col].unique():

        X_train, y_train, _, _ = train_test_split(adata, celltype, features, label_col=label_col, celltype_col=celltype_col)
        
        if penalty is None:
            model = LogisticRegression(solver='liblinear')
        else:
            model = LogisticRegression(penalty=penalty, solver='liblinear')
        model.fit(X_train, y_train)
        
        models[celltype] = model

    return models

def evaluate_models_for_celltypes(adata, models, features, celltype_col='celltype'):
    result_test = []
    result_score = []
    for celltype in adata.obs[celltype_col].unique():

        _, _, X_test, y_test = train_test_split(adata, celltype, features)
        result_test.extend(y_test)
        model = models[celltype]
        
        y_score = model.predict_proba(X_test)[:, 1]
        result_score.extend(y_score)
    
    return pd.DataFrame({'label': result_test, 'score': result_score, 'celltype': adata.obs[celltype_col]})

def select_optimal_num_features(
        adata, 
        label_col='is_sen', 
        celltype_col='celltype', 
        feat_select_plot=None,
        min_features=1500,
        max_features=2500,
        step=50
    ):
    
    feat_to_auc = []
    for n_features in range(min_features, max_features, step):
        features = univariate_feature_selection(adata, n_features=n_features, label_col=label_col)
        models = train_models_for_celltypes(adata, features, label_col=label_col, celltype_col=celltype_col)
        result = evaluate_models_for_celltypes(adata, models, features, celltype_col=celltype_col)
        precision, recall, _ = precision_recall_curve(result['label'], result['score'])
        auc_score = auc(recall, precision)
        feat_to_auc.append((n_features, auc_score))

    if feat_select_plot is not None:
        _log.info(f"Saving feature selection plot to {feat_select_plot}")
        df = pd.DataFrame(feat_to_auc, columns=['n_features', 'auc'])
        fig, ax = plt.subplots(1, 1, figsize=(5, 3), dpi=300)
        sns.lineplot(x='n_features', y='auc', data=df, ax=ax)
        ax.set_xlabel('Number of features', fontsize=16)
        ax.set_ylabel('AUC', fontsize=16)
        plt.tight_layout()
        plt.savefig(feat_select_plot, dpi=300)

    return max(feat_to_auc, key=lambda x: x[1])[0]

def get_common_features(adata, models, features, importance_threshold, n_models, gene_col='Gene_name'):

    def get_top_features_for_celltype(models, features, importance_threshold):
        top_features = {}
        for celltype, model in models.items():
            coef = model.coef_[0]
            top_features[celltype] = [features[np.abs(coef) > importance_threshold], coef[np.abs(coef) > importance_threshold]]
        
        return top_features

    def get_overlap(top_genes, n_models):
        gene_counts = Counter()
        for genes in top_genes.values():
            gene_counts.update(genes[0])

        common_genes = [gene for gene, count in gene_counts.items() if count >= n_models]
        return np.array(common_genes)
    
    def get_feature_importance(top_features, common_features):
        feature_importance = {}
        # if some of the cell types do not have common features, we need to use zero importance for them
        for celltype, genes in top_features.items():
            for feature in common_features:
                if feature in genes[0]:
                    feature_importance.setdefault(celltype, []).append(genes[1][genes[0] == feature][0])
                else:
                    feature_importance.setdefault(celltype, []).append(0)
        
        return feature_importance

    top_features = get_top_features_for_celltype(models, features, importance_threshold)
    _log.debug(f"Top features: {top_features}")
    common_features = get_overlap(top_features, n_models)
    _log.debug(f"Common features: {common_features}")
    feature_importance = get_feature_importance(top_features, common_features)
    _log.debug(f"Feature importance: {feature_importance}")
    mean_coefs = np.array(list(feature_importance.values())).mean(axis=0)
    _log.debug(f"Mean coefs: {mean_coefs}")
    gene_names = adata.var[gene_col][common_features]

    return pd.DataFrame({'gene': gene_names, 'coef': mean_coefs})


def main():

    # fix all seeds
    random.seed(42)
    np.random.seed(42)

    parser = argparse.ArgumentParser(description='Train a logistic regression model for each cell type')
    parser.add_argument('adata_path', type=str)
    parser.add_argument('--results_csv', type=str, default=None)
    parser.add_argument('--common_features_csv', type=str, default=None)
    parser.add_argument('--importance_threshold', type=float, default=0.01)
    parser.add_argument('--tuned_common_features_csv', type=str, default=None)
    parser.add_argument('--tuned_results_csv', type=str, default=None)
    parser.add_argument('--allowed_missing_samples', type=int, default=0)
    parser.add_argument('--blacklist_genes_csv', type=str, default=None, help='Path to blacklist genes file')
    parser.add_argument('--penalty', type=str, default=None)
    parser.add_argument('--quantile', type=float, default=0.8, help='Quantile for filtering genes')
    parser.add_argument('--num_features', type=int, default=1500)
    parser.add_argument('--allowed_missing_models', type=int, default=0)
    parser.add_argument('--pca_plot', type=str, default=None, help='Path to save PCA plot of selected features')
    parser.add_argument('--log', type=str, help='Path to log file', required=True)
    parser.add_argument('--log-level', type=str, help='Log level', default='INFO')
    args = parser.parse_args()

    set_logging(args.log, args.log_level)
    _log.debug(f"Command line arguments: {args}")

    adata = sc.read_h5ad(args.adata_path)
    _log.debug(f"Read AnnData object with shape {adata.X.shape}")

    if args.pca_plot is not None:
        _log.info(f"Saving PCA plot to {args.pca_plot}")
        adata_pca = adata.copy()

        sc.pp.normalize_total(adata_pca, target_sum=1e4)
        sc.pp.log1p(adata_pca)
        sc.pp.highly_variable_genes(adata_pca, n_top_genes=2000, subset=True)

        sc.pp.pca(adata_pca, n_comps=2, use_highly_variable=False, svd_solver='arpack')
        ax = sc.pl.pca(adata_pca, color=['celltype', 'is_sen'], show=False, size=50)
        data_type = args.pca_plot.split('.')[0]
        plt.title(f'PCA of selected features, {data_type} data', fontsize=16)
        plt.savefig(args.pca_plot, dpi=300)
        plt.close()
        del adata_pca

    allowed_missing_samples = min(adata.shape[0], args.allowed_missing_samples)
    _log.info(f"Feature pre-selection. Keep only genes with counts in at least {adata.shape[0] - allowed_missing_samples} samples")   
    adata = adata[:, adata.X.astype(bool).sum(0) >= adata.shape[0] - allowed_missing_samples]
    _log.debug(f"Filtered genes. Shape: {adata.X.shape}")

    if args.blacklist_genes_csv is None or args.blacklist_genes_csv == "None":
        _log.info(f"No cell cycle genes file provided. Skipping cell cycle genes filtering")
    else:
        _log.info(f"Loading cell cycle genes from {args.blacklist_genes_csv}")
        blacklist_genes = pd.read_csv(args.blacklist_genes_csv)
        _log.debug(f"Cell cycle genes: {blacklist_genes.shape}")
        _log.debug(f"Cell cycle genes: {blacklist_genes.head()}")

        # if intersection of cell_cycle_genes["Gene"] and adata.var_names is not empty, then remove them from adata
        blacklist_genes = blacklist_genes[blacklist_genes["Gene"].isin(adata.var["Gene_name"])]
        if blacklist_genes.shape[0] > 0:
            _log.info(f"Removing cell cycle genes from adata")
            _log.debug(f"Cell cycle genes: {blacklist_genes.shape}")
            _log.debug(f"Cell cycle genes: {blacklist_genes.head()}")
            adata = adata[:, ~adata.var["Gene_name"].isin(blacklist_genes["Gene"])]
            _log.debug(f"Filtered genes. Shape: {adata.X.shape}")
        else:
            _log.info(f"No cell cycle genes found in adata")

    _log.info(f"Feature pre-selection. Keep only genes with counts above {args.quantile} quantile")
    adata.var['total_counts'] = adata.X.sum(0)
    adata = adata[:, adata.var['total_counts'] > adata.var['total_counts'].quantile(args.quantile)]
    _log.debug(f"Filtered genes. Shape: {adata.X.shape}")

    _log.info(f"Do log1p transformation")
    sc.pp.log1p(adata)

    features = univariate_feature_selection(adata, n_features=args.num_features)
    _log.info(f"Selected {args.num_features} features")
    _log.debug(f"Selected features: {features}")
    if args.penalty == "None":
        args.penalty = None
    models = train_models_for_celltypes(adata, features, penalty=args.penalty)
    _log.info(f"Trained models for cell types")
    result = evaluate_models_for_celltypes(adata, models, features)
    _log.info(f"Evaluated models for cell types")
    allowed_missing_models = min(len(models), args.allowed_missing_models)
    common_features = get_common_features(adata, models, features, n_models = len(models) - allowed_missing_models, importance_threshold = args.importance_threshold)
    _log.info(f"Selected common features")

    if args.results_csv is not None:
        _log.info(f"Saving results to {args.results_csv}")
        result.to_csv(args.results_csv)

    if args.common_features_csv is not None:
        _log.info(f"Saving common features to {args.common_features_csv}")
        common_features.to_csv(args.common_features_csv)

    if args.tuned_common_features_csv is not None:

        coefs_models = train_models_for_celltypes(adata, common_features.index)
        coefs_results = evaluate_models_for_celltypes(adata, coefs_models, common_features.index)
        coefs_common_features = get_common_features(adata, coefs_models, common_features.index, n_models = 0, importance_threshold = 0)
        coefs_common_features.to_csv(args.tuned_common_features_csv)

        if args.tuned_results_csv is not None:
            _log.info(f"Saving tuned results to {args.tuned_results_csv}")
            coefs_results.to_csv(args.tuned_results_csv)


if __name__ == '__main__':
    main()
