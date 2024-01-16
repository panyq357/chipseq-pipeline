library(Biostrings)
library(GenomicRanges)

config <- list(
    ref_genome = snakemake@input$ref_genome,
    peak_file = snakemake@input$peak_file,
    top = snakemake@params$top,
    radius = snakemake@params$radius,
    out_fasta = snakemake@output$out_fasta
)

main <- function() {

    genome <- readDNAStringSet(config$ref_genome)

    peak_gr <- read_macs_xls(config$peak_file)

    target_dna <- extract_fasta_from_summit(peak_gr, genome=genome, top_n=config$top, radius=config$radius)

    writeXStringSet(target_dna, config$out_fasta)
}

extract_fasta_from_summit <- function(peak_gr, genome, top_n=500, radius=50) {

    rank_by_pvalue <- order(mcols(peak_gr)[["-log10(pvalue)"]], decreasing=T)

    top_gr <- peak_gr[rank_by_pvalue,][1:top_n, ]

    top_summit_gr <- GRanges(
        seqnames=seqnames(top_gr),
        ranges=IRanges(start=top_gr$abs_summit, end=top_gr$abs_summit),
        name=top_gr$name
    )
    target_gr <- flank(top_summit_gr, radius, both=T)

    target_dna <- genome[target_gr]
    names(target_dna) <- sprintf("%s|%s:%s-%s", target_gr$name, seqnames(target_gr), start(target_gr), end(target_gr))

    return(target_dna)
}

read_macs_xls <- function (macs_xls_path) {
    df <- read.table(macs_xls_path, header=T, check.names=F)
    gr <- GenomicRanges::makeGRangesFromDataFrame(df, keep.extra.columns=TRUE)
    return(gr)
}

main()

