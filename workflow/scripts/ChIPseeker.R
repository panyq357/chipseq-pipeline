log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")


library(ChIPseeker)
library(GenomicFeatures)


txdb <- txdbmaker::makeTxDbFromGFF(snakemake@input$gtf)

peak <- rtracklayer::import(snakemake@input$peak)

anno_res <- annotatePeak(peak = peak, TxDb = txdb)

gene_info <- readr::read_tsv(snakemake@input$gene_info) |> setNames(c("geneId", "symbol", "description"))

out_df <- as.data.frame(anno_res) |>
  dplyr::left_join(gene_info, by = "geneId")

readr::write_excel_csv(out_df, snakemake@output$csv)

pdf(snakemake@output$pdf)

  plotAvgProf2(
    peak = peak, TxDb = txdb,
    upstream = 1000, downstream = 1000,
    conf = 0.95
  ) |> print()

  plotPeakProf2(
    peak = peak, TxDb = txdb,
    upstream = rel(0.2), downstream = rel(0.2),
    conf = 0.95,
    by = "gene", type = "body",
    nbin = min(end(genes(txdb)) - start(genes(txdb))), ignore_strand = F
  ) |> print()

  plotAnnoPie(anno_res) |> print()

dev.off()
