library(ggplot2)

fc_stat <- lapply(snakemake@input$fc_list, readRDS) |>
    lapply(function(x) {
        df <- x$stat
        row.names(df) <- df[[1]]
        row <- df[2] |> t() |> as.data.frame()
        return(row)
    }) |>
    do.call(what=rbind, args=_)

fc_stat$FRiP <- fc_stat$Assigned / rowSums(fc_stat)
fc_stat$Peak <- sub("results/FRiP/(.+).fc.rds", "\\1", snakemake@input$fc_list)

pdf(snakemake@output$plot)
    print(
        ggplot(fc_stat) +
        geom_col(aes(y=Peak, x=FRiP)) +
        geom_text(aes(y=Peak, x=FRiP + 0.1, label=sprintf("%.2f%%", FRiP * 100))) +
        scale_x_continuous(limits=c(0, 1))
    )
dev.off()

