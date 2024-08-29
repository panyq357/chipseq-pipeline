
'''
Peaks overlap with what features?
'''
rule peak_to_gene:
    input:
        gtf = config["gtf"],
        peak = "results/macs_callpeak/{peak_type}/{peak}_peaks.xls"
    output:
        peak_to_gene = "results/peak_anno/{peak_type}/{peak}/{peak}.peak_to_gene.csv",
        feature_bar = "results/peak_anno/{peak_type}/{peak}/{peak}.feature_bar.svg"
    params:
        by = "gene_id",
        col_prefix = "Overlap."
    script:
        "../scripts/peak_to_gene.R"


'''
Highest peak in promoter/gene body?
'''
rule gene_to_peak:
    input:
        gtf = config["gtf"],
        peak = "results/macs_callpeak/{peak_type}/{peak}_peaks.{peak_type}Peak"
    output:
        gene_to_peak = "results/peak_anno/{peak_type}/{peak}/{peak}.gene_to_peak.csv"
    script:
        "../scripts/gene_to_peak.R"


rule oryzabase_anno:
    input:
        gene_to_peak = "results/peak_anno/{peak_type}/{peak}/{peak}.gene_to_peak.csv",
        oryzabase_xlsx = config["oryzabase"]["xlsx"],
        rap_anno = config["rap-anno"]
    output:
        anno_xlsx = "results/peak_anno/{peak_type}/{peak}/{peak}.min-signal-{signal}.oryzabase_anno.xlsx",
    params:
        oryzabase_sheet = config["oryzabase"]["sheet"],
        by = "promoter"
    script:
        "../scripts/oryzabase_anno_macs.R"

