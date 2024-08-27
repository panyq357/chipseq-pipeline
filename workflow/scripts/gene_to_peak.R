library(GenomicRanges)

config <- list(
    # input
    gtf = snakemake@input$gtf,
    peak = snakemake@input$peak,
    
    # output
    gene_to_peak = snakemake@output$gene_to_peak
)

main <- function() {
    
    gtf <- rtracklayer::import(config$gtf)

    gene <- subset(gtf, type == "gene")

    peak <- rtracklayer::import(config$peak)

    df <- data.frame(
        gene_id = gene$gene_id,
        promoter.max_signal = get_max_signal(promoters(gene), peak, mc.cores=config$threads),
        gene_body.max_signal = get_max_signal(gene, peak, mc.cores=config$threads)
    )

    write.csv(df, config$gene_to_peak, row.names=F)
}

get_max_signal <- function(gr, peak, ...) {

    hits <- findOverlaps(gr, peak)

    x <- split(peak$signalValue[hits@to], hits@from) |> lapply(max) |> unlist()

    values <- rep(0, length(gr))

    values[as.integer(names(x))] <- x

    return(values)
}

main()
