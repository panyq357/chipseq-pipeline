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
        gtf = lambda w: config["counts"][w.name]["gtf"],
        bams = lambda w: config["counts"][w.name]["bams"].values()
    output:
        counts_csv = "results/counts/{name}.mark_counts.csv",
        fc_rds = "results/counts/{name}.mark_counts.fc.rds"
    params:
        colnames = lambda w: config["counts"][w.name]["bams"].keys(),
        isPairedEnd = config["featureCounts"]["isPairedEnd"]
    threads:
        20
    script:
        "../scripts/mark_counts.R"


rule get_promoter:
    input:
        config["gtf"]
    output:
        "results/features/promoter.bed"
    script:
        "../scripts/get_promoter.R"


rule assign_mark_to_feature:
    input:
        mark_counts = "results/counts/{name}.mark_counts.csv",
        feature = "results/features/{feature}.bed",
    output:
        mark_counts = "results/counts/{name}.mark_counts.to-{feature}.csv",
    wildcard_constraints:
        feature = "[^-.]+"
    script:
        "../scripts/assign_mark_to_feature.R"


rule deseq:
    input:
        counts = "results/counts/{name}.{infix}.csv"
    output:
        csv = "results/counts/{name}.{infix}.deseq.csv",
        rds = "results/counts/{name}.{infix}.deseq.rds"
    wildcard_constraints:
        name = "[^.]+"
    params:
        levels = lambda w: config["deseq"][w.name]["levels"],
        sample_groups = lambda w: config["deseq"][w.name]["sample_groups"]
    script:
        "../scripts/deseq.R"

