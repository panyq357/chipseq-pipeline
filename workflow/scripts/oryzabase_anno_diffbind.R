library(clusterProfiler)


config <- list(
    peak_to_gene = snakemake@input$peak_to_gene,
    oryzabase_xlsx = snakemake@input$oryzabase_xlsx,
    rap_anno = snakemake@input$rap_anno,

    enricher_xlsx = snakemake@output$enricher_xlsx,
    enricher_rds = snakemake@output$enricher_rds,

    oryzabase_sheet = snakemake@params$oryzabase_sheet,
    fdr_th = as.numeric(snakemake@params$fdr_th),
    p_th = as.numeric(snakemake@params$p_th),
    abs_l2fc_th = log2(as.numeric(snakemake@params$fc_th))
)


main <- function(){

    peak <- read.csv(config$peak_to_gene, check.names=F)

    pass_th <- with(peak, FDR < config$fdr_th & `p-value` < config$p_th)
    peak_list <- list(
        All  = subset(peak, pass_th & abs(Fold) >  config$abs_l2fc_th),
        Up   = subset(peak, pass_th &     Fold  >  config$abs_l2fc_th),
        Down = subset(peak, pass_th &     Fold  < -config$abs_l2fc_th)
    )

    # Genes overlapped with diff peaks.
    gene_list <- lapply(peak_list, peak_to_gene)
    enricher_res_list <- lapply(gene_list, oryzabase_enricher, universe=NULL, oryzabase_xlsx=config$oryzabase_xlsx, sheet_list=config$oryzabase_sheet)

    enricher_res_list <- unlist(enricher_res_list, use.names=TRUE)

    saveRDS(enricher_res_list, config$enricher_rds)

    writexl::write_xlsx(
        c(
            lapply(enricher_res_list, as.data.frame),
            setNames(peak_list, sprintf("Peak.%s", names(peak_list)))
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
    gene <- peak[12:length(peak)] |>
        unlist() |> strsplit(",") |> unlist() |> unique() |> na.omit()
    return(gene)
}


main()

