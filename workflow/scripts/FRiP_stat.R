peak <- rtracklayer::import(snakemake@input$peak)

temp_gtf <- tempfile()
rtracklayer::export(peak, temp_gtf, format="gtf")

fc <- Rsubread::featureCounts(
    files=snakemake@input$treatment,
    annot.ext=temp_gtf,
    isGTFAnnotationFile=TRUE,
    GTF.featureType="sequence_feature",
    GTF.attrType="name",
    isPairedEnd=TRUE,
    nthreads = snakemake@threads
)

file.remove(temp_gtf)

saveRDS(fc, snakemake@output$fc)

