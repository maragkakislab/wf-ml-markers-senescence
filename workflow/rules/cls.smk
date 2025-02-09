
rule ml_classifier:
    input:
        anndata = os.path.join(ANALYSIS_DIR, "SenCat.normalized.h5ad"),
    output:
        results_csv = os.path.join(ANALYSIS_DIR, "SenCat.classification_results.csv"),
        common_features_csv = os.path.join(ANALYSIS_DIR, "SenCat.common_features.csv"),
        tuned_common_features_csv = os.path.join(ANALYSIS_DIR, "SenCat.tuned_common_features.csv"),
        tuned_results_csv = os.path.join(ANALYSIS_DIR, "SenCat.tuned_classification_results.csv"),
        feat_select_plot = os.path.join(PLOTS_DIR, "SenCat.feat_select_plot.pdf"),
    params:
        importance_threshold = 0
    log:
        os.path.join(LOG_DIR, "ml_classifier.log")
    conda:
        "../envs/data.yaml"
    shell:
        """
        python {workflow.basedir}/scripts/cls/ml_classifier.py \
            {input.anndata} \
            --results_csv {output.results_csv} \
            --common_features_csv {output.common_features_csv} \
            --feat_select_plot {output.feat_select_plot} \
            --importance_threshold {params.importance_threshold} \
            --tuned_common_features_csv {output.tuned_common_features_csv} \
            --tuned_results_csv {output.tuned_results_csv} \
            --log {log} \
            --log-level {LOG_LEVEL} \
            2>&1 | tee {log}
        """

rule marker_classifier:
    input:
        h5ad = os.path.join(ANALYSIS_DIR, "{counts_file}.normalized.h5ad"),
        markers = os.path.join(ANALYSIS_DIR, "SenCat.{markers}.csv")
    output:
        results_csv = os.path.join(ANALYSIS_DIR, "{counts_file}.{markers}_results.csv")
    log:
        os.path.join(LOG_DIR, "{counts_file}.{markers}.marker_classifier.log")
    conda:
        "../envs/data.yaml"
    shell:
        """
        python {workflow.basedir}/scripts/cls/marker_classifier.py \
            --input-h5ad {input.h5ad} \
            --markers {input.markers} \
            --output-results-csv {output.results_csv} \
            --log {log} \
            --log-level {LOG_LEVEL} \
            2>&1 | tee {log}
        """

use rule get_SenCID_predictions from sencid_workflow as sencid_get_SenCID_predictions with:
    output:
        rec = os.path.join(ANALYSIS_DIR, "{counts_file}.SenCID_results.csv")
    input: 
        txt = os.path.join(ANALYSIS_DIR, "{counts_file}.for_SenCID.txt")
    params:
        denoising = 'f'