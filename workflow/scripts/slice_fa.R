library(GenomicRanges)

genome_fa <- Biostrings::readDNAStringSet(snakemake@input$genome)

gr <- read.table(snakemake@input$peak, header=T, check.names=F) |> makeGRangesFromDataFrame(keep.extra.columns=T)
gr <- gr[order(mcols(gr)[[snakemake@wildcards$by]], decreasing=T),]
gr <- gr[as.integer(snakemake@wildcards$from):as.integer(snakemake@wildcards$to),]

radius <- as.integer(snakemake@wildcards$radius)

if (snakemake@wildcards$mode == "summit") {
    summit <- mcols(gr)[["abs_summit"]]
    df <- data.frame(chr=seqnames(gr), start=summit-radius, end=summit+radius)
    fa_gr <- data.frame(chr=seqnames(gr), start=summit-radius, end=summit+radius) |>
        makeGRangesFromDataFrame()
}

out_fa <- genome_fa[fa_gr]
names(out_fa) <- with(df, sprintf("%s:%d-%d", chr, start, end))
Biostrings::writeXStringSet(out_fa, snakemake@output$fa)
