config <- list(
    peak_info = data.frame(
        row.names=snakemake@params$peak_names,
        BED = snakemake@input$bed,
        BAM = snakemake@input$bam
    ),
    fc_stat = snakemake@output$fc_stat,
    plot = snakemake@output$plot
)

fc_stat <- apply(config$peak_info, 1, function(row) {

    peak <- rtracklayer::import(row["BED"])

    temp_gtf <- tempfile()
    rtracklayer::export(peak, temp_gtf, format="gtf")

    fc <- Rsubread::featureCounts(
        files=row["BAM"],
        annot.ext=temp_gtf,
        isGTFAnnotationFile=TRUE,
        GTF.featureType="sequence_feature",
        GTF.attrType="name",
        isPairedEnd=TRUE,
        nthreads = 10
    )

    file.remove(temp_gtf)

    return(setNames(fc$stat[[2]], fc$stat[[1]]))
})

write.csv(fc_stat, config$fc_stat)

png(config$plot)
apply(fc_stat, 2, function(col) {
    col["Assigned"] / sum(col)
}) |> barplot(ylim = c(0, 1))
dev.off()
