config <- list(
    peak = snakemake@input$peak,
    bam = snakemake@input$bam,
    fc = snakemake@output$fc,
    nthreads = snakemake@threads
)

main <- function() {

    peak <- rtracklayer::import(config$peak)

    temp_gtf <- tempfile()
    rtracklayer::export(peak, temp_gtf, format="gtf")

    fc <- Rsubread::featureCounts(
        files=config$bam,
        annot.ext=temp_gtf,
        isGTFAnnotationFile=TRUE,
        GTF.featureType="sequence_feature",
        GTF.attrType="name",
        isPairedEnd=TRUE,
        nthreads = config$nthreads
    )

    file.remove(temp_gtf)

    saveRDS(fc, config$fc)
}

main()

