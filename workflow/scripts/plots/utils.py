import logging
from matplotlib.colors import LinearSegmentedColormap

def set_logging(_log, log_file, log_level):
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

def get_marker_gene_values(adata, celltype, marker_genes):

    gene_markers = adata[adata.obs["celltype"] == celltype, marker_genes.index].to_df()
    # map id to names. using `marker_genes` because `gene_markers` are in a different order
    gene_markers.columns = gene_markers.columns.map(lambda x: marker_genes.loc[x, 'gene'])
    # sort columns by the order in `marker_genes`
    gene_markers = gene_markers.reindex(columns=marker_genes['gene'].values)

    return gene_markers

def get_cmap():
    return LinearSegmentedColormap.from_list('mycmap', ['#001AFF', 'white', '#FF0000'])