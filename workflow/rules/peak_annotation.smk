rule ChIPseeker:
    input:
        gtf = config["gtf"],
        gene_info = config["gene_info"],
        peak = "results/peak_calling/{peak_type}/{peak_name}_peaks.{peak_type}Peak"
    output:
        csv = "results/peak_calling/{peak_type}/{peak_name}_peaks.ChIPseeker.csv",
        pdf = "results/peak_calling/{peak_type}/{peak_name}_peaks.ChIPseeker.pdf"
    log:
        "results/peak_calling/{peak_type}/{peak_name}_peaks.ChIPseeker.log"
    script:
        "../scripts/ChIPseeker.R"


rule enricher_promoter:
    input:
        csv = "results/peak_calling/{peak_type}/{peak_name}_peaks.ChIPseeker.csv",
        enrichment_data = [file for dataset in config["enrichment_data"].values() for file in dataset.values()]
    output:
        xlsx = "results/peak_calling/{peak_type}/{peak_name}_peaks.ChIPseeker.enricher-promoter.xlsx",
        pdf = "results/peak_calling/{peak_type}/{peak_name}_peaks.ChIPseeker.enricher-promoter.pdf"
    params:
        enrichment_data = config["enrichment_data"],
        annotation_grep = "Promoter"
    log:
        "results/peak_calling/{peak_type}/{peak_name}_peaks.ChIPseeker.enricher-promoter.log"
    script:
        "../scripts/peak_enricher.R"

