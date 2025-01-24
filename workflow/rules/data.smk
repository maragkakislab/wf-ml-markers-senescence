rule assemble_anndata:
    input:
        count_files = expand(os.path.join(DATA_DIR, "3_{cell_type}_all_results_annot.txt"), cell_type=CELL_TYPES)
    output:
        anndata = os.path.join(ANALYSIS_DIR, "counts.h5ad"),
        # save as txt file for SenCID - it cannot load h5ad object created by newer versions of scanpy
        # index is Gene_name -> rows with NaNs are removed
        txt = os.path.join(ANALYSIS_DIR, "counts.txt")
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