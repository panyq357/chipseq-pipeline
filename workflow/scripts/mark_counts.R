fc <- Rsubread::featureCounts(
  files = snakemake@input$bams,
  annot.ext = snakemake@input$gtf,
  isGTFAnnotationFile = TRUE,
  GTF.featureType = "sequence_feature",
  GTF.attrType = "name",
  isPairedEnd = snakemake@params$isPairedEnd,
  nthreads = snakemake@threads
)

counts <- fc$counts
colnames(counts) <- snakemake@params$colnames

saveRDS(fc, snakemake@output$fc_rds)
write.csv(counts, snakemake@output$counts_csv)

