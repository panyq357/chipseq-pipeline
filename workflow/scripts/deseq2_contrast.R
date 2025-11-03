# Use connection to output both STDOUT and STDERR to same file.
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")   # STDOUT
sink(log, type = "message")  # STDERR

library(DESeq2)


counts <- read.csv(snakemake@input$counts, check.names = FALSE, row.names = 1)

is_integer <- sapply(counts, class) == "integer"

count_columns <- counts[ is_integer]
extra_columns <- counts[!is_integer]

coldata <- read.csv(snakemake@input$coldata, check.names = FALSE)

levels <- snakemake@params$levels

# Remove samples not in levels
for (name in names(levels)) {
  coldata <- coldata[coldata[[name]] %in% levels[[name]], ]
}

# Change character to factor
for (col in names(coldata)) {
  if (class(coldata[[col]]) == "character") {
    coldata[[col]] <- factor(coldata[[col]], levels=unique(coldata[[col]]))
  }
}

dds <- DESeqDataSetFromMatrix(
  count_columns[coldata$sample],
  coldata,
  as.formula(paste("~", col))
)

# Pre-filtering.
dds <- dds[rowSums(counts(dds)) >= 10, ]

# DESeq.
dds <- DESeq(dds)

# Use lfcShrink to get a better Log2FoldChange for downstream GSEA.
# colnames(coef(dds))[2] looks like: "Group_mutent_vs_wildtype".
res <- lfcShrink(dds, colnames(coef(dds))[2], type="apeglm")

res_df <- data.frame(res, check.names=FALSE)

# Combine normalized counts.
nc <- counts(dds, normalized=TRUE)

# Add gene annotation
gene_info <- readr::read_tsv(snakemake@input$gene_info, col_names=c("geneID", "symbol", "description"))

res_df <- cbind(
  res_df,
  nc[row.names(res_df), ],
  gene_info[match(row.names(res_df), gene_info$geneID), -1],
  extra_columns[row.names(res_df), , drop=FALSE]
)

# Sort by padj
res_df <- res_df[order(res_df$pvalue, decreasing=FALSE), ]

## Output DESeq results.
tibble::rownames_to_column(res_df, var = "geneID") |> readr::write_csv(snakemake@output[[1]])
