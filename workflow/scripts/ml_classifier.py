from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import precision_recall_curve, auc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import anndata as an
from collections import Counter
import logging

from utils import set_logging, get_default_parser

def univariate_feature_selection(adata, n_features=500, label_col='label'):
    logging.info(f"Selecting {n_features} features")
    selector = SelectKBest(f_classif, k=n_features)
    selector.fit(adata.X, adata.obs[label_col])
    features = selector.get_support(indices=True)
    return features

def train_test_split(adata, celltype, features, label_col='label', celltype_col='celltype'):

    logging.debug(f"Splitting data for celltype {celltype}")

    train = adata[adata.obs[celltype_col] != celltype, features]
    logging.debug(f"Training data shape: {train.X.shape}")
    test = adata[adata.obs[celltype_col] == celltype, features]
    logging.debug(f"Testing data shape: {test.X.shape}")
    
    return train.X, train.obs[label_col], test.X, test.obs[label_col]

def train_models_for_celltypes(adata, features, label_col='label', celltype_col='celltype'):
    models = {}
    for celltype in adata.obs[celltype_col].unique():
        logging.debug(f"Training model for celltype {celltype}")

        X_train, y_train, _, _ = train_test_split(adata, celltype, features, label_col=label_col, celltype_col=celltype_col)
        
        model = LogisticRegression(penalty='l1', solver='liblinear')
        model.fit(X_train, y_train)
        
        models[celltype] = model

    return models

def evaluate_models_for_celltypes(adata, models, features, celltype_col='celltype'):
    result_test = []
    result_score = []
    for celltype in adata.obs[celltype_col].unique():
        logging.debug(f"Evaluating model for celltype {celltype}")

        _, _, X_test, y_test = train_test_split(adata, celltype, features)
        result_test.extend(y_test)
        logging.debug(f"Testing data shape: {X_test.shape}")
        logging.debug(f"Testing data labels: {y_test}")
        model = models[celltype]
        
        y_score = model.predict_proba(X_test)[:, 1]
        logging.debug(f"Testing data scores: {y_score}")
        result_score.extend(y_score)
    
    return pd.DataFrame({'label': result_test, 'score': result_score, 'celltype': adata.obs[celltype_col]})

def select_optimal_num_features(
        adata, 
        label_col='label', 
        celltype_col='celltype', 
        feat_select_plot=None,
        min_features=50,
        max_features=2000,
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
        logging.info(f"Saving feature selection plot to {feat_select_plot}")
        df = pd.DataFrame(feat_to_auc, columns=['n_features', 'auc'])
        fig, ax = plt.subplots(1, 1, figsize=(5, 3), dpi=300)
        sns.lineplot(x='n_features', y='auc', data=df, ax=ax)
        ax.set_xlabel('Number of features', fontsize=16)
        ax.set_ylabel('AUC', fontsize=16)
        plt.tight_layout()
        plt.savefig(feat_select_plot, dpi=300)

    return max(feat_to_auc, key=lambda x: x[1])[0]

def get_common_features(adata, models, features, importance_threshold, gene_col='Gene_name'):

    def get_top_features_for_celltype(models, features, importance_threshold):
        top_features = {}
        for celltype, model in models.items():
            coef = model.coef_[0]
            top_features[celltype] = [features[np.abs(coef) > importance_threshold], coef[np.abs(coef) > importance_threshold]]
        
        return top_features

    def get_overlap(top_genes):
        gene_counts = Counter()
        for genes in top_genes.values():
            gene_counts.update(genes[0])

        common_genes = [gene for gene, count in gene_counts.items() if count == 14]
        return np.array(common_genes)
    
    def get_feature_importance(top_features, common_features):
        feature_importance = {}
        for celltype, genes in top_features.items():
            feature_importance[celltype] = genes[1][np.isin(genes[0], common_features)]
        
        return feature_importance

    top_features = get_top_features_for_celltype(models, features, importance_threshold)
    common_features = get_overlap(top_features)
    feature_importance = get_feature_importance(top_features, common_features)
    mean_coefs = np.array(list(feature_importance.values())).mean(axis=0)
    gene_names = adata.var[gene_col][common_features]

    return pd.DataFrame({'gene': gene_names, 'mean_coef': mean_coefs})


def main():
    parser = get_default_parser(__doc__)
    parser.add_argument('adata_path', type=str)
    parser.add_argument('--feat_select_plot', type=str, default=None)
    parser.add_argument('--results_csv', type=str, default=None)
    parser.add_argument('--common_features_csv', type=str, default=None)
    parser.add_argument('--importance_threshold', type=float, default=0.01)
    args = parser.parse_args()

    set_logging(args.log, args.log_level)
    logging.debug(f"Command line arguments: {args}")

    adata = an.read_h5ad(args.adata_path)
    logging.debug(f"Read AnnData object with shape {adata.X.shape}")
    n_features = select_optimal_num_features(adata, feat_select_plot=args.feat_select_plot)
    features = univariate_feature_selection(adata, n_features=n_features)
    logging.info(f"Selected {n_features} features")
    logging.debug(f"Selected features: {features}")
    models = train_models_for_celltypes(adata, features)
    logging.info(f"Trained models for cell types")
    result = evaluate_models_for_celltypes(adata, models, features)
    logging.info(f"Evaluated models for cell types")
    common_features = get_common_features(adata, models, features, importance_threshold=args.importance_threshold)
    logging.info(f"Selected common features")

    if args.results_csv is not None:
        logging.info(f"Saving results to {args.results_csv}")
        result.to_csv(args.results_csv, sep='\t')

    if args.common_features_csv is not None:
        logging.info(f"Saving common features to {args.common_features_csv}")
        common_features.to_csv(args.common_features_csv, sep='\t')

if __name__ == '__main__':
    main()
