library(Biostrings)
library(GenomicRanges)


config <- list(
    ref_genome = snakemake@input$ref_genome,
    peak_xls = snakemake@input$peak_xls,

    ref_point = snakemake@params[[1]]$ref_point,
    rank_by = snakemake@params[[1]]$rank_by,
    top = snakemake@params[[1]]$top,
    radius = snakemake@params[[1]]$radius,
    p_th = snakemake@params[[1]]$p_th,
    q_th = snakemake@params[[1]]$q_th,
    fold_range = snakemake@params[[1]]$fold_range,

    out_fasta = snakemake@output$out_fasta
)


main <- function() {

    genome <- readDNAStringSet(config$ref_genome)

    peak_gr <- read_macs_xls(config$peak_xls)

    peak_gr <- subset(
        peak_gr,
        `-log10(pvalue)` > -log10(config$p_th) &
            `-log10(qvalue)` > -log10(config$q_th) &
            fold_enrichment > config$fold_range[[1]] &
            fold_enrichment < config$fold_range[[2]]
    )

    target_gr <- peak_gr[order(mcols(peak_gr)[[config$rank_by]], decreasing=T),]
    target_gr <- target_gr[1:ifelse(config$top > length(target_gr), length(target_gr), config$top),]

    if (config$ref_point == "center") {

        ranges(target_gr) <- IRanges(
            start = (end(target_gr) + start(target_gr)) / 2 - config$radius,
            end = (end(target_gr) + start(target_gr)) / 2 + config$radius
        )

    } else if (config$ref_point == "summit") {
        
        ranges(target_gr) <- IRanges(
            start = target_gr$abs_summit - config$radius,
            end = target_gr$abs_summit + config$radius
        )
    }

    target_dna <- genome[target_gr]
    names(target_dna) <- sprintf("%s|%s:%s-%s", target_gr$name, seqnames(target_gr), start(target_gr), end(target_gr))

    writeXStringSet(target_dna, config$out_fasta)
}


read_macs_xls <- function (macs_xls_path) {
    df <- read.table(macs_xls_path, header=T, check.names=F)
    gr <- GenomicRanges::makeGRangesFromDataFrame(df, keep.extra.columns=TRUE)
    return(gr)
}


main()

