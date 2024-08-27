

rule peak_to_gene_macs:
    input:
        gtf = config["gtf"],
        peak = "results/macs_callpeak/{peak_type}/{peak}_peaks.xls"
    output:
        peak_to_gene = "results/peak_anno/{peak_type}/{peak}/{peak}.peak_to_gene_macs.csv",
        feature_bar = "results/peak_anno/{peak_type}/{peak}/{peak}.feature_bar.svg"
    params:
        anno = "gene_id"
    script:
        "../scripts/peak_to_gene.R"

