library(clusterProfiler)

config <- list(
    peak_to_gene = snakemake@input$peak_to_gene,
    oryzabase_xlsx = snakemake@input$oryzabase_xlsx,
    rap_anno = snakemake@input$rap_anno,

    enricher_xlsx = snakemake@output$enricher_xlsx,
    enricher_rds = snakemake@output$enricher_rds,

    oryzabase_sheet = snakemake@params$oryzabase_sheet,
    fold_enrichment_range = snakemake@params$fold_enrichment_range
)


main <- function(){

    # All genes have annotation as universe.
    universe <- NULL

    peak <- read.csv(config$peak_to_gene, check.names=F)
    peak <- subset(peak, fold_enrichment > config$fold_enrichment_range[[1]] & fold_enrichment < config$fold_enrichment_range[[2]])

    # Genes overlapped with diff peaks.
    gene <- peak_to_gene(peak)
    enricher_res_list <- oryzabase_enricher(gene, universe, config$oryzabase_xlsx, config$oryzabase_sheet)
    saveRDS(enricher_res_list, config$enricher_rds)

    rap_anno <- readr::read_tsv(config$rap_anno)
    gene_used <- subset(rap_anno, Locus_ID %in% gene, select = c("Locus_ID", "Oryzabase Gene Symbol Synonym(s)", "Oryzabase Gene Name Synonym(s)", "Description")) 
    gene_used <- gene_used[!duplicated(gene_used$Locus_ID), ]

    writexl::write_xlsx(c(lapply(enricher_res_list, as.data.frame), list(GeneUsed=gene_used)), config$enricher_xlsx)

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
    gene <- peak[13:length(peak)] |>
        unlist() |> strsplit(",") |> unlist() |> unique() |> na.omit()
    return(gene)
}

main()

