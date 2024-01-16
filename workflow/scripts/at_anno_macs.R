library(clusterProfiler)

config <- list(
    peak_to_gene = snakemake@input$peak_to_gene,
    oryzabase_xlsx = snakemake@input$oryzabase_xlsx,
    xlsx = snakemake@output$xlsx,
    rds = snakemake@output$rds
)

main <- function(){

    # All genes have annotation as universe.
    universe <- NULL

    peak <- read.csv(config$peak_to_gene, check.names=F)

    gene <- get_gene_id(peak)

    # Genes overlapped with diff peaks.
    enrich_go_res_list <- at_enrich_go(gene)
    enrich_kegg_res <- at_enrich_kegg(gene)

    names(enrich_go_res_list) <- sprintf("Enrich.GO.%s", names(enrich_go_res_list))

    all_res <- c(enrich_go_res_list, list(Enrich.KEGG=enrich_kegg_res))

    writexl::write_xlsx(lapply(all_res, as.data.frame), config$xlsx)
    saveRDS(all_res, config$rds)

}

get_gene_id <- function(peak, regex="AT[12345CM]G[0-9]{5}") {
    gene <- unlist(peak) |>
        stringr::str_extract_all(regex) |>
        unlist() |> unique() |> na.omit()
    return(gene)
}

at_enrich_go <- function(gene, universe=NULL) {
    ego <- list()
    for (ont in c("BP", "MF", "CC")) {
        ego[[ont]] <- clusterProfiler::enrichGO(
            gene = gene,
            universe = universe,
            OrgDb = org.At.tair.db::org.At.tair.db,
            keyType = 'TAIR',
            ont = ont,
            pvalueCutoff = 1,
            qvalueCutoff = 1
        )
    }
    return(ego)
}

at_enrich_kegg <- function(gene, universe=NULL) {
    kk <- clusterProfiler::enrichKEGG(
        gene = gene,
        universe = universe,
        organism = 'ath',
        pvalueCutoff = 1,
        qvalueCutoff = 1
    )
    return(kk)
}

at_gse_go <- function(l2fc) {
    gsego <- list()
    for (ont in c("BP", "MF", "CC")) {
        gsego[[ont]] <- clusterProfiler::gseGO(
            geneList = l2fc,
            OrgDb = org.At.tair.db::org.At.tair.db,
            keyType = 'TAIR',
            ont = ont,
            pvalueCutoff = 1
        )
    }
    return(gsego)
}

at_gse_kegg <- function(l2fc) {
    gsekk <- clusterProfiler::gseKEGG(
        geneList = l2fc,
        organism = 'ath',
        pvalueCutoff = 1
    )
    return(gsekk)
}

main()

