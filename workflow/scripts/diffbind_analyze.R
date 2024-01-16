library(DiffBind)

config <- list(
    counts_rds = snakemake@input$counts_rds,
    fdr_th = snakemake@params$fdr_th,
    p_th = snakemake@params$p_th,
    ma_plot = snakemake@output$ma_plot,
    volcano = snakemake@output$volcano,
    report_rds = snakemake@output$report_rds,
    report_csv = snakemake@output$report_csv,
    fdr_peak = snakemake@output$fdr_peak,
    fdr_peak_up = snakemake@output$fdr_peak_up,
    fdr_peak_down = snakemake@output$fdr_peak_down
)

main <- function() {


    counts <- readRDS(config$counts_rds)

    model <- dba.contrast(counts)
    result <- dba.analyze(model)

    svg(config$ma_plot)
    dba.plotMA(result)
    dev.off()

    svg(config$volcano)
    dba.plotVolcano(result)
    dev.off()

    report <- dba.report(result, th=1)  # th: p-value threshold

    out_df <- as.data.frame(GenomicRanges::granges(report))
    out_df <- cbind(out_df, mcols(report))

    saveRDS(report, config$report_rds)
    write.csv(out_df, config$report_csv, row.names=F)

    # Prepare BED files for deeptools visualization
    bed_df <- granges_to_bed_df(report)

    is_fdr_and_p <- (report$FDR < config$fdr_th) & (report$`p-value` < config$p_th)
    is_up <- report$Fold > 0
    is_down <- report$Fold < 0

    df_to_bed(bed_df[is_fdr_and_p,], config$fdr_peak)
    df_to_bed(bed_df[is_fdr_and_p & is_up,], config$fdr_peak_up)
    df_to_bed(bed_df[is_fdr_and_p & is_down,], config$fdr_peak_down)
}

granges_to_bed_df <- function(report) {
    df <- data.frame(
        seqnames=seqnames(report),
        starts=start(report)-1,
        ends=end(report),
        names=rep(".", length(report)),
        scores=rep(".", length(report)),
        strands=strand(report)
    )
    return(df)
}

df_to_bed <- function(df, bed_path) {
    write.table(df, bed_path, quote=F, sep="\t", row.names=F, col.names=F)
}

main()

