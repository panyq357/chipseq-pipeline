library(GenomicRanges)

genome_fa <- Biostrings::readDNAStringSet(snakemake@input$genome)

gr <- rtracklayer::import(snakemake@input$peak)
gr <- gr[order(mcols(gr)[[snakemake@wildcards$by]], decreasing=T),]
gr <- gr[as.integer(snakemake@wildcards$from):as.integer(snakemake@wildcards$to),]

radius <- as.integer(snakemake@wildcards$radius)

if (snakemake@wildcards$mode == "summit") {
    summit <- start(gr) + mcols(gr)[["peak"]]
    df <- data.frame(chr=seqnames(gr), start=summit-radius, end=summit+radius)
    fa_gr <- data.frame(chr=seqnames(gr), start=summit-radius, end=summit+radius) |>
        makeGRangesFromDataFrame()
}

out_fa <- genome_fa[fa_gr]
names(out_fa) <- with(df, sprintf("%s:%d-%d", chr, start, end))
Biostrings::writeXStringSet(out_fa, snakemake@output$fa)
