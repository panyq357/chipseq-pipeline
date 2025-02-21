library(GenomicRanges)

config <- list(
    # input
    gtf = snakemake@input$gtf,
    peak = snakemake@input$peak,
    
    # output
    gene_to_peak = snakemake@output$gene_to_peak
)

main <- function() {
    
    gtf <- rtracklayer::import(snakemake@input$gtf)
    peak <- rtracklayer::import(snakemake@input$peak)

    gene <- split(gtf, gtf$gene_id) |> range() |> unlist()

    df <- data.frame(
        gene_id = names(gene),
        promoter.max_signal = get_max_signal(promoters(gene), peak),
        gene_body.max_signal = get_max_signal(gene, peak)
    )

    write.csv(df, config$gene_to_peak, row.names=F)
}

get_max_signal <- function(gr, peak) {

    hits <- findOverlaps(gr, peak)

    x <- split(peak$signalValue[hits@to], hits@from) |> lapply(max) |> unlist()

    values <- rep(0, length(gr))

    values[as.integer(names(x))] <- x

    return(values)
}

main()
