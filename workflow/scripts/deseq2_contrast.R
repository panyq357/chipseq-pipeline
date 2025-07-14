log <- file(snakemake@log[[1]], open="wt")
sink(log)  # STDOUT
sink(log, type = "message")  # STDERR


library(DESeq2)


# Filter and order coldata.
coldata <- read.csv(snakemake@input$coldata, check.names = F)

levels <- snakemake@params$levels

for (col in names(levels)[length(levels):1]) {
  coldata <- subset(coldata, coldata[[col]] %in% snakemake@params$levels[[col]])
  coldata[[col]] <- factor(coldata[[col]], levels = snakemake@params$levels[[col]])
}

# Use first col to do analysis.
coldata <- coldata[order(coldata[[col]]),]


design <- as.formula(paste("~", col))


cts <- read.csv(snakemake@input$counts, check.names = F, row.names = 1)

# Save extra columns (e.g. OriginMark) for appending to final result.
extra_column <- cts[! names(cts) %in% coldata$sample]

cts <- cts[coldata$sample]

dds <- DESeqDataSetFromMatrix(cts, coldata, design)

# Pre-filtering.
dds <- dds[rowSums(counts(dds)) >= 10,]

# DESeq.
dds <- DESeq(dds)

# Use lfcShrink to get a better Log2FoldChange for downstream GSEA.
# colnames(coef(dds))[2] looks like: "Group_mutent_vs_wildtype".
res <- lfcShrink(dds, colnames(coef(dds))[2], type="apeglm")

res_df <- data.frame(res, check.names=F)

# Combine normalized counts.
nc <- counts(dds, normalized=TRUE)

# Add gene annotation
gene_info <- readr::read_tsv(snakemake@input$gene_info, col_names=c("geneID", "symbol", "description"))
gene_info <- gene_info[match(row.names(res_df), gene_info$geneID),]

res_df <- cbind(res_df, nc[match(row.names(res_df), row.names(nc)),], gene_info[-1], extra_column[row.names(res_df),, drop=FALSE])

# Sort by padj
res_df <- res_df[order(res_df$pvalue, decreasing=FALSE),]

## Output DESeq results.
tibble::rownames_to_column(res_df, var = "geneID") |>
  readr::write_csv(snakemake@output[[1]])
