rule k_wise_intersection:
    input:
        lambda w: config["peak_sets"][w.name]
    output:
        "results/k_wise_intersection/{name}.{k}-wise-intersection.bed"
    script:
        "../scripts/k_wise_intersection.R"


rule bed_to_gtf:
    input:
        "{prefix}.bed"
    output:
        "{prefix}.gtf"
    script:
        "../scripts/bed_to_gtf.R"


rule mark_counts:
    input:
        gtf = lambda w: config["mark_count_jobs"][w.name]["gtf"],
        bams = lambda w: config["mark_count_jobs"][w.name]["bams"].values()
    output:
        counts_csv = "results/mark_counts/{name}.mark_counts.csv",
    params:
        colnames = lambda w: config["mark_count_jobs"][w.name]["bams"].keys()
    threads:
        config["featureCounts"]["threads"]
    script:
        "../scripts/mark_counts.R"


rule get_genebody:
    input:
        config["genome_gtf"]
    output:
        "results/features/genebody.bed"
    script:
        "../scripts/get_genebody.R"


rule assign_mark_to_feature:
    input:
        mark_counts = "results/mark_counts/{name}.mark_counts.csv",
        feature = "results/features/{feature}.bed",
    output:
        "results/mark_counts/{name}.mark_counts.assign-to-{feature}.csv",
    wildcard_constraints:
        feature = "[^-.]+"
    script:
        "../scripts/assign_mark_to_feature.R"


rule deseq2_contrast:
    input:
        counts = lambda w: config["de_analysis_jobs"]["contrast"][w.name]["counts"],
        coldata = config["coldata"],
        gene_info = config["gene_info"]
    output:
        "results/de_analysis/contrast/{name}.deseq2_contrast.csv"
    log:
        "results/de_analysis/contrast/{name}.deseq2_contrast.log"
    params:
        levels = lambda w: config["de_analysis_jobs"]["contrast"][w.name]["levels"]
    script:
        "../scripts/deseq2_contrast.R"


rule volcano_plot:
    input:
        "results/de_analysis/{type}/{name}.deseq2_{type}.csv"
    output:
        "results/de_analysis/{type}/{name}.deseq2_{type}.{p_column}-{p_th}.fc-{fc_th}.volcano.svg"
    log:
        "results/de_analysis/{type}/{name}.deseq2_{type}.{p_column}-{p_th}.fc-{fc_th}.volcano.log"
    script:
        "../scripts/volcano_plot.R"


rule de_res_enricher:
    input:
        de_res = "results/de_analysis/{type}/{name}.deseq2_{type}.csv",
        enrichment_data = [dataset.values() for dataset in config["enrichment_data"].values()]
    output:
        xlsx = "results/de_analysis/{type}/{name}.deseq2_{type}.{p_column}-{p_th}.fc-{fc_th}.enricher.xlsx",
        pdf = "results/de_analysis/{type}/{name}.deseq2_{type}.{p_column}-{p_th}.fc-{fc_th}.enricher.pdf"
    params:
        enrichment_data = config["enrichment_data"]
    log:
        "results/de_analysis/{type}/{name}.deseq2_{type}.{p_column}-{p_th}.fc-{fc_th}.enricher.log"
    script:
        "../scripts/de_res_enricher.R"

