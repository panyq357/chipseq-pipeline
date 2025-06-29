counts <- read.csv(snakemake@input$mark_counts, check.names=F, row.names=1)

feature <- rtracklayer::import(snakemake@input$feature)

match_res <- stringr::str_match(row.names(counts), "([^-]+)-([^-]+)-([^-]+)")

gr <- GenomicRanges::makeGRangesFromDataFrame(
  data.frame(
    seqnames = match_res[,2],
    start = match_res[,3],
    end = match_res[,4],
    name = match_res[,1]
  ),
  keep.extra.columns = TRUE
)

hits <- GenomicRanges::findOverlaps(gr, feature)

counts_to_feature <- split(counts[hits@from,], feature$name[hits@to]) |>
  lapply(function(df) apply(df, 2, sum)) |>
  do.call(rbind, args=_) |>
  as.data.frame()

write.csv(counts_to_feature, snakemake@output[[1]])
