rule qc_plots:
    input:
        anndata = os.path.join(ANALYSIS_DIR, "{counts_file}.h5ad"),
    output:
        library_size_plot = os.path.join(PLOTS_DIR, "{counts_file}.library_size.pdf"),
        mito_plot = os.path.join(PLOTS_DIR, "{counts_file}.mito_plot.pdf"),
        ribo_plot = os.path.join(PLOTS_DIR, "{counts_file}.ribo_plot.pdf"),
        hb_plot = os.path.join(PLOTS_DIR, "{counts_file}.hb_plot.pdf"),
    log:
        os.path.join(LOG_DIR, "{counts_file}_qc_plots.log")
    conda:
        "../envs/data.yaml"
    shell:
        """
        python workflow/scripts/plots/qc.py \
            {input.anndata} \
            --library_size_path {output.library_size_plot} \
            --mito_path {output.mito_plot} \
            --ribo_path {output.ribo_plot} \
            --hb_path {output.hb_plot} \
            --log {log} \
            --log-level {LOG_LEVEL} \
            2>&1 | tee {log}
        """

rule preds_plots:
    input:
        rec = os.path.join(ANALYSIS_DIR, "{counts_type}.{predictio_method}_results.csv"),
        h5ad = os.path.join(ANALYSIS_DIR, "{counts_type}.h5ad"),
    output:
        pdf = os.path.join(PLOTS_DIR, "{counts_type}.{predictio_method}_results.pdf")
    log:
        os.path.join(LOG_DIR, "{counts_type}.{predictio_method}_preds_plots.log")
    conda:
        "../envs/data.yaml"
    shell:
        """
        python {workflow.basedir}/scripts/plots/predictions.py \
            --results-csv {input.rec} \
            --input-h5ad {input.h5ad} \
            --output-plot {output.pdf} \
            --log {log} \
            --log-level {LOG_LEVEL} \
        """

rule gene_markers:
    input:
        h5ad = os.path.join(ANALYSIS_DIR, "{counts_file_input}.h5ad"),
        csv = os.path.join(ANALYSIS_DIR, "{counts_file_markers}_common_features.csv"),
    output:
        plot = os.path.join(PLOTS_DIR, "{counts_file_input}_counts.{counts_file_markers}_gene_markers.pdf")
    conda:
        "../envs/data.yaml"
    log:
        os.path.join(LOG_DIR, "gene_markers.{counts_file_input}.{counts_file_markers}.log")
    shell:
        """
        python workflow/scripts/plots/gene_markers.py \
            --input-h5ad {input.h5ad} \
            --gene-markers-csv {input.csv} \
            --output-plot {output.plot} \
            --log {log} \
            --log-level {LOG_LEVEL} \
        """    

rule gene_markers_with_ml_marker:
    input:
        h5ad = os.path.join(ANALYSIS_DIR, "{counts_file_input}.h5ad"),
        csv = os.path.join(ANALYSIS_DIR, "{counts_file_markers}_common_features.csv"),
        ml_csv = os.path.join(ANALYSIS_DIR, "{counts_file_ml_markers}_tuned_common_features.csv"), 
    output:
        plot = os.path.join(PLOTS_DIR, "{counts_file_input}_counts.{counts_file_markers}_gene_markers.{counts_file_ml_markers}_ml_markers.pdf")
    conda:
        "../envs/data.yaml"
    log:
        os.path.join(LOG_DIR, "gene_markers.{counts_file_input}.{counts_file_markers}.{counts_file_ml_markers}.log")
    shell:
        """
        python workflow/scripts/plots/gene_markers_with_ml_marker.py \
            --input-h5ad {input.h5ad} \
            --gene-markers-csv {input.csv} \
            --ml-marker-csv {input.ml_csv} \
            --output-plot {output.plot} \
            --log {log} \
            --log-level {LOG_LEVEL} \
        """

rule ml_marker:
    input:
        h5ad = os.path.join(ANALYSIS_DIR, "{counts_file_input}.h5ad"),
        ml_csv = os.path.join(ANALYSIS_DIR, "{counts_file_ml_markers}_tuned_common_features.csv"), 
    output:
        plot = os.path.join(PLOTS_DIR, "{counts_file_input}_counts.{counts_file_ml_markers}_ml_markers.pdf")
    conda:
        "../envs/data.yaml"
    log:
        os.path.join(LOG_DIR, "ml_marker.{counts_file_input}.{counts_file_ml_markers}.log")
    shell:
        """
        python workflow/scripts/plots/ml_marker.py \
            --input-h5ad {input.h5ad} \
            --ml-marker-csv {input.ml_csv} \
            --output-plot {output.plot} \
            --log {log} \
            --log-level {LOG_LEVEL} \
        """