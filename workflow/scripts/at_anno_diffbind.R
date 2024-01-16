library(clusterProfiler)


config <- list(
    peak_to_gene = snakemake@input$peak_to_gene,
    oryzabase_xlsx = snakemake@input$oryzabase_xlsx,
    enricher_xlsx = snakemake@output$enricher_xlsx,
    enricher_rds = snakemake@output$enricher_rds,
    fdr_th = as.numeric(snakemake@params$fdr_th),
    abs_l2fc_th = log2(as.numeric(snakemake@params$fc_th)),
    oryzabase_sheet = snakemake@params$oryzabase_sheet
)


main <- function(){

    # All genes have annotation as universe.
    universe <- NULL

    peak <- read.csv(config$peak_to_gene, check.names=F)

    peak_list <- list(
        all = subset(peak, FDR < config$fdr_th & abs(Fold) > config$abs_l2fc_th),
        up = subset(peak, FDR < config$fdr_th & Fold > config$abs_l2fc_th),
        down = subset(peak, FDR < config$fdr_th & Fold < -config$abs_l2fc_th)
    )

    # Genes overlapped with diff peaks.
    gene_list <- lapply(peak_list, peak_to_gene)
    enricher_res_list <- lapply(gene_list, oryzabase_enricher, universe, config$oryzabase_xlsx, config$oryzabase_sheet)

    names(enricher_res_list$all) <- sprintf("All_%s", names(enricher_res_list$all))
    names(enricher_res_list$up) <- sprintf("Up_%s", names(enricher_res_list$up))
    names(enricher_res_list$down) <- sprintf("Down_%s", names(enricher_res_list$down))

    saveRDS(enricher_res_list, config$enricher_rds)

    writexl::write_xlsx(
        c(
            lapply(enricher_res_list$all, as.data.frame), list(All_Peak_Used = peak_list$all),
            lapply(enricher_res_list$up, as.data.frame), list(Up_Peak_Used = peak_list$up),
            lapply(enricher_res_list$down, as.data.frame), list(Down_Peak_Used = peak_list$down)
        ),
        config$enricher_xlsx
    )

}


df_enricher <- function(gene, universe, onto_df, gene_col, term_col, name_col) {
    res <- clusterProfiler::enricher(
        gene = gene,
        universe = universe,
        TERM2GENE = onto_df[c(term_col, gene_col)],
        TERM2NAME = onto_df[c(term_col, name_col)],
        pvalueCutoff = 1,
        qvalueCutoff = 1
    )
    return(res)
}

df_GSEA <- function(ranked_gene_list, onto_df, gene_col, term_col, name_col) {
    res <- clusterProfiler::GSEA(
        geneList = ranked_gene_list,
        TERM2GENE = onto_df[c(term_col, gene_col)],
        TERM2NAME = onto_df[c(term_col, name_col)],
        pvalueCutoff = 1,
    )
    return(res)
}


oryzabase_enricher <- function(gene, universe, oryzabase_xlsx, sheet_list) {
    res_list <- list()
    for (sheet in sheet_list) {
        onto_df <- readxl::read_excel(oryzabase_xlsx, sheet=sheet)
        res <- df_enricher(gene, universe, onto_df, "GeneID", "OntoID", "Description")
        res_list[[sheet]] <- res
    }
    return(res_list)
}

oryzabase_GSEA <- function(ranked_gene_list, oryzabase_xlsx, sheet_list) {
    res_list <- list()
    for (sheet in sheet_list) {
        onto_df <- readxl::read_excel(oryzabase_xlsx, sheet=sheet)
        res <- df_GSEA(ranked_gene_list, onto_df, "GeneID", "OntoID", "Description")
        res_list[[sheet]] <- res
    }
    return(res_list)
}

peak_to_gene <- function(peak) {
    tx_matrix <- peak[12:length(peak)]
    gene <- unlist(tx_matrix) |>
        stringr::str_extract_all("Os[0-9]{2}t[0-9]{7}") |>
        unlist() |> unique() |> sub("t", "g", x=_)
    return(gene)
}

main()

