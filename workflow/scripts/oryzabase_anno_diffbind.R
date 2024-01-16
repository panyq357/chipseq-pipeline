library(clusterProfiler)


config <- list(
    peak_to_gene = snakemake@input$peak_to_gene,
    xlsx = snakemake@output$xlsx,
    rds = snakemake@output$rds,
    fdr_th = as.numeric(snakemake@params$fdr_th),
    p_th = as.numeric(snakemake@params$p_th),
    abs_l2fc_th = log2(as.numeric(snakemake@params$fc_th))
)


main <- function(){

    peak <- read.csv(config$peak_to_gene, check.names=F)

    pass_th <- with(peak, FDR < config$fdr_th & `p-value` < config$p_th)
    peak_list <- list(
        All  = subset(peak, pass_th & abs(Fold) >  config$abs_l2fc_th),
        Up   = subset(peak, pass_th & Fold      >  config$abs_l2fc_th),
        Down = subset(peak, pass_th & Fold      < -config$abs_l2fc_th)
    )

    # Genes overlapped with diff peaks.
    gene_list <- lapply(peak_list, get_gene_id)
    enrich_go_res_list <- lapply(gene_list, at_enrich_go)

    names(enrich_go_res_list$All) <- sprintf("Enrich.GO.%s", names(enrich_go_res_list$All))
    names(enrich_go_res_list$Up) <- sprintf("Enrich.GO.%s", names(enrich_go_res_list$Up))
    names(enrich_go_res_list$Down) <- sprintf("Enrich.GO.%s", names(enrich_go_res_list$Down))

    enrich_go_res_list <- unlist(enrich_go_res_list, use.names=TRUE)

    saveRDS(enrich_go_res_list, config$rds)

    writexl::write_xlsx(
        c(
            lapply(enrich_go_res_list, as.data.frame),
            setNames(peak_list, sprintf("Peak.%s", names(peak_list)))
        ),
        config$xlsx
    )

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

