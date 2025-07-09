log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

enrichment_data <- lapply(snakemake@params$enrichment_data, function(dataset) {
  lapply(dataset, readr::read_tsv, col_names=FALSE)
})

chipseeker_res <- readr::read_csv(snakemake@input$csv)

gene <- chipseeker_res[grep(snakemake@params$annotation_grep, chipseeker_res$annotation),]$geneId |>
  unique()

enrich_res <- lapply(enrichment_data, function(dataset) {
  clusterProfiler::enricher(
    gene = gene,
    universe = NULL,
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    TERM2GENE = dataset$term_to_gene,
    TERM2NAME = dataset$term_to_description
  )
})

c(lapply(enrich_res, as.data.frame), list(GeneUsed = chipseeker_res[match(gene, chipseeker_res$geneId), c("geneId", "symbol", "description")])) |>
  writexl::write_xlsx(snakemake@output$xlsx)

pdf(snakemake@output$pdf)
for (name in names(enrich_res)) {
  enrichplot::dotplot(enrich_res[[name]], title = name) |> print()
}
dev.off()
