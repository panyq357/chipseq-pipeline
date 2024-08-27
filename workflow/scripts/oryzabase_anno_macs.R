library(clusterProfiler)

config <- list(
    gene_to_peak = snakemake@input$gene_to_peak,
    oryzabase_xlsx = snakemake@input$oryzabase_xlsx,
    rap_anno = snakemake@input$rap_anno,

    anno_xlsx = snakemake@output$anno_xlsx,

    oryzabase_sheet = snakemake@params$oryzabase_sheet,
    by = snakemake@params$by,

    min_signal = snakemake@wildcards$signal
)


main <- function(){

    # All genes have annotation as universe.
    universe <- NULL

    gene_to_peak <- read.csv(config$gene_to_peak, check.names=F)
    gene_to_peak <- gene_to_peak[c("gene_id", sprintf("%s.max_signal", config$by))]
    signal <- setNames(gene_to_peak[[2]], gene_to_peak[[1]]) |> sort(decreasing=T)

    gene <- names(signal[signal > config$min_signal])

    # Genes overlapped with diff peaks.
    enricher_res_list <- oryzabase_enricher(gene, universe, config$oryzabase_xlsx, config$oryzabase_sheet)

    rap_anno <- readr::read_tsv(config$rap_anno)

    gene_used <- data.frame(GeneID = gene, Signal = signal[gene])
    gene_used[c("Symbol", "Description")] <- rap_anno[match(gene_used$GeneID, rap_anno$Locus_ID), c("Oryzabase Gene Symbol Synonym(s)", "Description")] 

    writexl::write_xlsx(c(lapply(enricher_res_list, as.data.frame), list(GeneUsed=gene_used)), config$anno_xlsx)

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

main()

