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
        rec = os.path.join(ANALYSIS_DIR, "{counts_file}.{predictio_method}_results.csv"),
        h5ad = os.path.join(ANALYSIS_DIR, "{counts_file}.h5ad"),
    output:
        pdf = os.path.join(PLOTS_DIR, "{counts_file}.{predictio_method}_results.pdf")
    conda:
        "../envs/data.yaml"
    shell:
        """
        python workflow/scripts/plots/predictions.py {input.rec} {output.pdf} {input.h5ad}
        """