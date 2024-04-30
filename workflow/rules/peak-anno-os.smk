configfile: "config/peak-anno-os-config.yaml"


rule oryzabase_anno_macs:
    input:
        peak_to_gene = "results/peak_anno/{job}/{job}.peak_to_gene_macs.csv",
        oryzabase_xlsx = config["oryzabase-xlsx"],
        rap_anno = config["rap-anno"]
    output:
        enricher_xlsx = "results/peak_anno/{job}/{job}.oryzabase_anno_macs.xlsx",
        enricher_rds = "results/peak_anno/{job}/{job}.oryzabase_anno_macs.rds"
    params:
        oryzabase_sheet = config["oryzabase-sheet"],
        fold_enrichment_range = [5, 50]
    script:
        "../scripts/oryzabase_anno_macs.R"

