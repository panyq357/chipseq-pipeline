counts <- read.csv(snakemake@input[[1]], row.names=1, check.names=FALSE)

is_integer <- sapply(counts, class) == "integer"

# Reuse counts, replace counts with cpm.
counts[is_integer] <- edgeR::cpm(counts[is_integer])

write.csv(counts, snakemake@output[[1]])
