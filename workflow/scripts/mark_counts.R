# Use connection to output both STDOUT and STDERR to same file.
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")   # STDOUT
sink(log, type = "message")  # STDERR

fc <- Rsubread::featureCounts(
  files = snakemake@input$bams,
  annot.ext = snakemake@input$gtf,
  isGTFAnnotationFile = TRUE,
  GTF.featureType = "sequence_feature",
  GTF.attrType = "name",
  isPairedEnd = length(grep("paired", snakemake@input$bams)) == length(snakemake@input$bams),
  nthreads = snakemake@threads
)

counts <- fc$counts

colnames(counts) <- snakemake@params$colnames

write.csv(counts, snakemake@output$counts_csv)
