rule anndata_from_excel:
    input:
        counts = lambda wilds: os.path.join(DATA_DIR, config["INPUT_COUNTS"][wilds.counts_type]["file"])
    output:
        anndata = os.path.join(ANALYSIS_DIR, "{counts_type}.h5ad"),
        # save as txt file for SenCID - it cannot load h5ad object created by newer versions of scanpy
        # index is Gene_name -> rows with NaNs are removed
        txt = os.path.join(ANALYSIS_DIR, "{counts_type}.for_SenCID.txt")
    params:
        var_columns = lambda wilds: config["INPUT_COUNTS"][wilds.counts_type]["var_columns"],
        index_column = lambda wilds: config["INPUT_COUNTS"][wilds.counts_type]["index_column"],
        log_level = LOG_LEVEL
    conda:
        "../envs/data.yaml"
    log:
        os.path.join(LOG_DIR, "{counts_type}.anndata_from_excel.log")
    shell:
        """
        python {workflow.basedir}/scripts/data/anndata_from_excel.py \
            --input-excel {input.counts:q} \
            --output-h5ad {output.anndata:q} \
            --output-txt {output.txt:q} \
            --var-columns {params.var_columns:q} \
            --index-column {params.index_column:q} \
            --log {log:q} \
            --log-level {LOG_LEVEL} \
        """

rule assemble_anndata:
    input:
        count_files = expand(os.path.join(DATA_DIR, "3_{cell_type}_all_results_annot.txt"), cell_type=CELL_TYPES)
    output:
        anndata = os.path.join(ANALYSIS_DIR, "SenCat.h5ad"),
        # save as txt file for SenCID - it cannot load h5ad object created by newer versions of scanpy
        # index is Gene_name -> rows with NaNs are removed
        txt = os.path.join(ANALYSIS_DIR, "SenCat.for_SenCID.txt")
    log:
        os.path.join(LOG_DIR, "assemble_anndata.log")
    conda:
        "../envs/data.yaml"
    shell:
        """
        python {workflow.basedir}/scripts/data/assemble_anndata.py \
            {input.count_files} \
            --output_h5ad {output.anndata} \
            --output_txt {output.txt} \
            --log {log} \
            --log-level {LOG_LEVEL} \
            2>&1 | tee {log}
        """

rule build_gene_id_table_from_ensembl:
    output:
        os.path.join(DATA_DIR, config["ASSEMBLY"], 'gene_ids_' + config["ARCHIVE_NAME"] + '.txt')
    params:
        link = 'http://' + config["ARCHIVE_NAME"] + '.ensembl.org/biomart/martservice?query=',
        xml = '<?xml version="1.0" encoding="UTF-8"?>',
        qopen = '<!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "1" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >',
        dopen = '<Dataset name = "' + config["ENSEMBL_ASSEMBLY_TO_SPECIES_NAME"][config["ASSEMBLY"]] + '_gene_ensembl" interface = "default" >',
        attr = "".join(['<Attribute name = "'+ a +'" />' for a in ['ensembl_gene_id', 'external_gene_name']]),
        dclose = '</Dataset>',
        qclose = '</Query>'
    shell:
        """
        wget -O {output} '{params.link}{params.xml}{params.qopen}{params.dopen}{params.attr}{params.dclose}{params.qclose}'
        """ 

rule anndata_from_txt:
    output:
        h5ad = os.path.join(ANALYSIS_DIR,  "Matts.h5ad"),
        # save as txt file for SenCID - it cannot load h5ad object created by newer versions of scanpy
        # index is Gene_name -> rows with NaNs are removed
        txt = os.path.join(ANALYSIS_DIR,  "Matts.for_SenCID.txt")
    input:
        counts = os.path.join(DATA_DIR, "23_all_sample_gene_counts.txt"),
        mapping = os.path.join(DATA_DIR, config["ASSEMBLY"], 'gene_ids_' + config["ARCHIVE_NAME"] + '.txt'),
        metadata = os.path.join(DATA_DIR, "23_coldata.xlsx")
    conda:
        "../envs/data.yaml"
    params:
        input_gene_col = "ENSG",
        mapping_gene_col = "Gene stable ID",
        mapping_name_col = "Gene name"
    log:
        os.path.join(LOG_DIR, "anndata_from_txt.log")
    shell:
        """
        python {workflow.basedir}/scripts/data/anndata_from_txt.py \
            --input-tsv {input.counts} \
            --mapping-tsv {input.mapping} \
            --metadata-excel {input.metadata} \
            --input-gene-col {params.input_gene_col} \
            --mapping-gene-col {params.mapping_gene_col:q} \
            --mapping-name-col {params.mapping_name_col:q} \
            --output-h5ad {output.h5ad} \
            --output-txt {output.txt} \
            --log {log} \
            --log-level {LOG_LEVEL} \
        """

rule normalize_counts:
    input:
        anndata = os.path.join(ANALYSIS_DIR, "{counts_file}.h5ad"),
    output:
        anndata = os.path.join(ANALYSIS_DIR, "{counts_file}.normalized.h5ad"),
    params:
        design = lambda wilds: config["NORMALIZATION_DESIGN"][wilds.counts_file],
    log:
        os.path.join(LOG_DIR, "{counts_file}.normalize_counts.log")
    conda:
        "../envs/pydeseq2.yaml"
    shell:
        """
        python {workflow.basedir}/scripts/data/normalize_counts.py \
            --input-h5ad {input.anndata} \
            --design {params.design:q} \
            --output-h5ad {output.anndata} \
            --log {log} \
            --log-level {LOG_LEVEL} \
            2>&1 | tee {log}
        """

rule prepare_common_sen_markers:
    output:
        os.path.join(ANALYSIS_DIR, "{sen_markers}.common_features.csv")
    params:
        markers = lambda wildcards: COMMON_SEN_MARKERS[wildcards.sen_markers]
    shell:
        """
        cp {params.markers} {output}   
        """