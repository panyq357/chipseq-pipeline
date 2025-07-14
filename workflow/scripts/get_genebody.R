library(GenomicRanges)

gtf <- rtracklayer::import(snakemake@input[[1]])

genebody <- subset(gtf, type == "gene")

mcols(genebody) <- mcols(genebody)["gene_id"] |> setNames("name")

rtracklayer::export(genebody, snakemake@output[[1]], format = "bed")
