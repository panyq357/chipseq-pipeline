peaks <- lapply(snakemake@input, rtracklayer::import)

all_intersect <- function(gr_list) {
  if (length(gr_list) == 1) {
    return(gr_list[[1]])
  } else {
    return(Reduce(GenomicRanges::intersect, gr_list))
  }
}

all_union <- function(gr_list) {
  if (length(gr_list) == 1) {
    return(gr_list[[1]])
  } else {
    return(Reduce(GenomicRanges::union, gr_list))
  }
}

out <- combn(peaks, snakemake@params$k, all_intersect, simplify=FALSE) |> all_union()

GenomicRanges::mcols(out)$name <- paste(
  GenomicRanges::seqnames(out),
  GenomicRanges::start(out),
  GenomicRanges::end(out),
  sep = "-"
)

rtracklayer::export(out, snakemake@output[[1]])
