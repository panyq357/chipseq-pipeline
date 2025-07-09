log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

txdb <- txdbmaker::makeTxDbFromGFF(snakemake@input$gtf)

peak <- rtracklayer::import(snakemake@input$peak)

anno_res <- ChIPseeker::annotatePeak(peak = peak, TxDb = txdb)

gene_info <- readr::read_tsv(snakemake@input$gene_info) |> setNames(c("geneId", "symbol", "description"))

out_df <- as.data.frame(anno_res) |>
  dplyr::left_join(gene_info, by = "geneId")

readr::write_excel_csv(out_df, snakemake@output$csv)

pdf(snakemake@output$pdf)
ChIPseeker::plotAvgProf2(peak = peak, TxDb = txdb, upstream = 1000, downstream = 1000, conf = 0.95) |> print()
ChIPseeker::plotAnnoPie(anno_res) |> print()
dev.off()
