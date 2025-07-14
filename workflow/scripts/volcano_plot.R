log <- file(snakemake@log[[1]], open="wt")
sink(log)  # STDOUT
sink(log, type = "message")  # STDERR


library(ggplot2)


config <- list(
  de_res = snakemake@input[[1]],
  p_column = snakemake@wildcards$p_column,
  p_th = as.numeric(snakemake@wildcards$p_th),
  l2fc_th = log2(as.numeric(snakemake@wildcards$fc_th)),
  name = snakemake@wildcards$name,
  outfile = snakemake@output[[1]]
)


de_res <- read.csv(config$de_res)

de_res$color <- "none"

is_up <- de_res[[config$p_column]] < config$p_th & de_res$log2FoldChange > config$l2fc_th
is_down <- de_res[[config$p_column]] < config$p_th & de_res$log2FoldChange < -config$l2fc_th

de_res$color[is_up] <- "up"
de_res$color[is_down] <- "down"

svg(config$outfile, width = 4, height = 4)

print(
  ggplot(de_res) +
    geom_point(aes(x = log2FoldChange, y = -log10(.data[[config$p_column]]), color = color)) +
    theme_bw() +
    scale_color_manual(
      name = "",
      values = c("none" = "grey", "down" = "blue", "up" = "red"),
      labels = c(
        "none" = sprintf("none: %s", sum(de_res$color == "none")),
        "down" = sprintf("down: %s", sum(de_res$color == "down")),
        "up" = sprintf("up: %s", sum(de_res$color == "up"))
      )
    ) +
    geom_vline(xintercept = c(-config$l2fc_th, config$l2fc_th), linetype="dashed") +
    geom_hline(yintercept = -log10(config$p_th), linetype="dashed") +
    theme(legend.position = "top") +
    labs(title = config$name, y = sprintf("-log10(%s)", config$p_column))
)

dev.off()
