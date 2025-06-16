gtf <- rtracklayer::import(snakemake@input[[1]])

promoter <- gtf |> subset(type == "gene") |> GenomicRanges::promoters()

GenomicRanges::mcols(promoter) <- data.frame(name = promoter$gene_id)

rtracklayer::export(promoter, snakemake@output[[1]])
