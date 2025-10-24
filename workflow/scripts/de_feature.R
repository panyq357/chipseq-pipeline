de_res <- readr::read_csv(snakemake@input$de_res)

p_column <- snakemake@wildcards$p_column
p_th <- as.numeric(snakemake@wildcards$p_th)
fc_th <- as.numeric(snakemake@wildcards$fc_th)

is_up   <- de_res[[p_column]] < p_th & de_res$log2FoldChange >  log2(fc_th)
is_down <- de_res[[p_column]] < p_th & de_res$log2FoldChange < -log2(fc_th)

feature <- rtracklayer::import(snakemake@input$feature)

feature_up <- feature[feature$name %in% de_res$geneID[is_up]]
feature_down <- feature[feature$name %in% de_res$geneID[is_down]]

rtracklayer::export(feature_up, snakemake@output$up)
rtracklayer::export(feature_down, snakemake@output$down)
