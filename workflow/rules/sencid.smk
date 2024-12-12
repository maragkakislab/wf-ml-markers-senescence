rule sencid_all:
    input:
        rec = os.path.join(ANALYSIS_DIR, "RecommendSID_Index_counts.txt")

rule get_SenCID_predictions:
    output:
        rec = os.path.join(ANALYSIS_DIR, "counts_SenCID_results.txt")
    input: 
        txt = os.path.join(ANALYSIS_DIR, "counts.txt")
    conda:
        "../envs/SenCID.yaml"
    params:
        output_dir = os.path.join(ANALYSIS_DIR),
        denoising = False
    shell:
        """
        SenCID --filepath {input.txt} --denoising {params.denoising} --fileclass txt --binarize t --output_dir {params.output_dir}
        """

rule plot_predictions:
    input:
        rec = os.path.join(ANALYSIS_DIR, "counts_SenCID_results.txt"),
    output:
        pdf = os.path.join(PLOTS_DIR, "SenCID_predictions.pdf")
    conda:
        "../envs/data.yaml"
    shell:
        """
        python workflow/scripts/plot_predictions.py {input.rec} {output.pdf}
        """