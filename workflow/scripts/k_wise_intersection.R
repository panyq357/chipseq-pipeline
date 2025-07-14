peak_set <- lapply(snakemake@input, rtracklayer::import)

# Use combn to generate a list of all combination of k-wise intersection.
# Then apply intersect to all combinations.
overlap_list <- combn(peak_set, as.integer(snakemake@wildcards$k), simplify=F) |>
  lapply(function(x) Reduce(GenomicRanges::intersect, x))

# Take union of all combinations.
if (length(overlap_list) == 1) {
    overlap <- overlap_list[[1]]
} else {
    overlap <- Reduce(GenomicRanges::union, overlap_list)
}

GenomicRanges::mcols(overlap)$name <- paste(
  GenomicRanges::seqnames(overlap),
  GenomicRanges::start(overlap),
  GenomicRanges::end(overlap),
  sep = "-"
)

rtracklayer::export(overlap, snakemake@output[[1]])
