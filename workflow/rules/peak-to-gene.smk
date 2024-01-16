configfile: "config/peak-to-gene-config.yaml"


rule peak_to_gene_diffbind:
    input:
        gtf = config["peak-to-gene-gtf"],
        peak = "results/diffbind_analyze/{job}/{job}.report.rds"
    output:
        peak_to_gene = "results/peak_anno/{job}/{job}.peak_to_gene_diffbind.csv",
        feature_bar = "results/peak_anno/{job}/{job}.feature_bar.svg"
    params:
        anno = "gene_id"
    script:
        "../scripts/peak_to_gene.R"  # Can read .xls or .rds


rule peak_to_gene_macs:
    input:
        gtf = config["peak-to-gene-gtf"],
        peak = "results/macs_callpeak/{peak}/{peak}_peaks.xls"
    output:
        peak_to_gene = "results/peak_anno/{peak}/{peak}.peak_to_gene_macs.csv",
        feature_bar = "results/peak_anno/{peak}/{peak}.feature_bar.svg"
    params:
        anno = "gene_id"
    script:
        "../scripts/peak_to_gene.R"

