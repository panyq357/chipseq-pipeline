rule at_anno_macs:
    input:
        peak_to_gene = "results/peak_anno/{job}/{job}.peak_to_gene_macs.csv"
    output:
        xlsx = "results/peak_anno/{job}/{job}.at_anno_macs.xlsx",
        rds = "results/peak_anno/{job}/{job}.at_anno_macs.rds"
    script:
        "../scripts/at_anno_macs.R"


rule at_anno_diffbind:
    input:
        peak_to_gene = "results/peak_anno/{job}/{job}.peak_to_gene_diffbind.csv"
    output:
        xlsx = "results/peak_anno/{job}/{job}.at_anno_diffbind.xlsx",
        rds = "results/peak_anno/{job}/{job}.at_anno_diffbind.rds"
    params:
        fdr_th = 1,
        p_th = 0.05,
        fc_th = 1
    script:
        "../scripts/oryzabase_anno_diffbind.R"


