log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

enrichment_data <- lapply(snakemake@params$enrichment_data, function(dataset) {
  lapply(dataset, readr::read_tsv, col_names=FALSE)
})

de_res <- readr::read_csv(snakemake@input$de_res)

# Use all gene detected in RNA-seq as background.
universe <- de_res$geneID

p_column <- snakemake@wildcards$p_column
p_th <- snakemake@wildcards$p_th |> as.numeric()
l2fc_th <- snakemake@wildcards$fc_th |> as.numeric() |> log2()

de_res$isDE <- "No"
de_res$isDE[de_res[[p_column]] < p_th & de_res$log2FoldChange >  l2fc_th] <- "Up"
de_res$isDE[de_res[[p_column]] < p_th & de_res$log2FoldChange < -l2fc_th] <- "Down"

gene_list <- list(
    All  = de_res$geneID[de_res$isDE != "No"],
    Up   = de_res$geneID[de_res$isDE == "Up"],
    Down = de_res$geneID[de_res$isDE == "Down"]
)


enrich_res <- lapply(gene_list, function(gene) {
  lapply(enrichment_data, function(dataset) {
    clusterProfiler::enricher(
      gene = gene,
      universe = universe,
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      TERM2GENE = dataset$term_to_gene,
      TERM2NAME = dataset$term_to_description
    )
  })
}) |> unlist()

c(lapply(enrich_res, as.data.frame), list(GeneUsed = de_res)) |>
  writexl::write_xlsx(snakemake@output$xlsx)

pdf(snakemake@output$pdf)
for (name in names(enrich_res)) {
  enrichplot::dotplot(enrich_res[[name]], title = name) |> print()
}
dev.off()

