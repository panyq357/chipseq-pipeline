

rule oryzabase_anno_macs:
    input:
        peak_to_gene = "results/peak_anno/{peak_type}/{peak}/{peak}.peak_to_gene_macs.csv",
        oryzabase_xlsx = config["oryzabase"]["xlsx"],
        rap_anno = config["rap-anno"]
    output:
        enricher_xlsx = "results/peak_anno/{peak_type}/{peak}/{peak}.oryzabase_anno_macs.xlsx",
    params:
        oryzabase_sheet = config["oryzabase"]["sheet"],
    script:
        "../scripts/oryzabase_anno_macs.R"

