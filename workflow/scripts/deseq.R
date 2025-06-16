make_coldata <- function(sample_groups, levels) {
  group_vect <- c()
  for (group in names(sample_groups)) {
    group_vect <- c(group_vect, rep(group, length(sample_groups[[group]])))
  }
  coldata <- data.frame(Group=factor(group_vect, levels=levels))
  row.names(coldata) <- unlist(sample_groups)
  return(coldata)
}

coldata <- make_coldata(snakemake@params$sample_groups, snakemake@params$levels)

counts <- read.csv(snakemake@input$counts, row.names=1, check.names=F)

cts <- counts[,row.names(coldata)]

dds <- DESeq2::DESeqDataSetFromMatrix(cts, coldata, ~ Group)

dds <- DESeq2::DESeq(dds)

res_df <- as.data.frame(DESeq2::results(dds))

res_df <- cbind(res_df, DESeq2::counts(dds, normalized=T))

res_df <- res_df[order(res_df$padj),]

write.csv(res_df, snakemake@output$csv)

saveRDS(dds, snakemake@output$rds)
