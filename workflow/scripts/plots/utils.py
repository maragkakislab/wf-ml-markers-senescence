import logging

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
    gene_markers.columns = marker_genes['gene']
    gene_markers.index = adata[adata.obs["celltype"] == celltype].obs['treatment']
    return gene_markers