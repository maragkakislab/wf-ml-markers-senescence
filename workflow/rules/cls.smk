
rule ml_classifier:
    input:
        anndata = os.path.join(ANALYSIS_DIR, "{input_counts}.h5ad"),
    output:
        results_csv = os.path.join(ANALYSIS_DIR, "{input_counts}.{ml_classifier}.classification_results.csv"),
        common_features_csv = os.path.join(ANALYSIS_DIR, "{input_counts}_{ml_classifier}_common_features.csv"),
        tuned_common_features_csv = os.path.join(ANALYSIS_DIR, "{input_counts}_{ml_classifier}_tuned_common_features.csv"),
        tuned_results_csv = os.path.join(ANALYSIS_DIR, "{input_counts}.{ml_classifier}.tuned_classification_results.csv"),
    params:
        importance_threshold = 0,
        allowed_missing_samples = lambda wildcards: config["ML_CLASSIFIER"][wildcards.ml_classifier]["allowed_missing_samples"],
        num_features = lambda wildcards: config["ML_CLASSIFIER"][wildcards.ml_classifier]["num_features"],
        allowed_missing_models = lambda wildcards: config["ML_CLASSIFIER"][wildcards.ml_classifier]["allowed_missing_models"],
    log:
        os.path.join(LOG_DIR, "{input_counts}.{ml_classifier}.ml_classifier.log")
    conda:
        "../envs/data.yaml"
    shell:
        """
        python {workflow.basedir}/scripts/cls/ml_classifier.py \
            {input.anndata} \
            --results_csv {output.results_csv} \
            --common_features_csv {output.common_features_csv} \
            --importance_threshold {params.importance_threshold} \
            --tuned_common_features_csv {output.tuned_common_features_csv} \
            --tuned_results_csv {output.tuned_results_csv} \
            --allowed_missing_samples {params.allowed_missing_samples} \
            --num_features {params.num_features} \
            --allowed_missing_models {params.allowed_missing_models} \
            --log {log} \
            --log-level {LOG_LEVEL} \
            2>&1 | tee {log}
        """

rule marker_classifier:
    input:
        h5ad = os.path.join(ANALYSIS_DIR, "{input_counts}.h5ad"),
        markers = os.path.join(ANALYSIS_DIR, "{counts_type_markers}_{markers}.csv")
    output:
        results_csv = os.path.join(ANALYSIS_DIR, "{input_counts}.{counts_type_markers}_{markers}.results.csv")
    log:
        os.path.join(LOG_DIR, "{input_counts}.{counts_type_markers}_{markers}.marker_classifier.log")
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
        rec = os.path.join(ANALYSIS_DIR, "{counts_type}.SenCID_results.csv")
    input: 
        txt = os.path.join(ANALYSIS_DIR, "{counts_type}.for_SenCID.txt")
    params:
        denoising = 'f'