# description: Use DiffBind to get consensus peakset and counts
# date: 2023-10-12
# author: panyq357

library(DiffBind)

config <- list(
    sample_sheet = snakemake@input$sample_sheet,
    counts_rds = snakemake@output$counts_rds,
    peakset_csv = snakemake@output$peakset_csv,
    threads = snakemake@threads
)

main <- function() {

    sample_sheet <- read.csv(config$sample_sheet)

    peaks <- dba(sampleSheet=sample_sheet)
        
    # Use less cores, and less memory
    peaks$config$cores <- config$threads

    # Time consumming, don't re-run.
    counts <- dba.count(peaks, summits=FALSE, bParallel=TRUE)
    counts <- dba.normalize(counts)

    peakset <- dba.peakset(counts, bRetrieve=TRUE)

    saveRDS(counts, config$counts_rds)
    write.csv(peakset, config$peakset_csv, row.names=F)
}

main()

