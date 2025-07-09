cbind(
  file = unlist(snakemake@input), 
  sapply(snakemake@input, function(file) {
    match_res <- readLines(file) |>
      stringr::str_match("(\\d+)\\s\\+\\s\\d+\\s([^(]+?)\\s*(\\(([^%]+%)?.*\\))?$")

    x <- setNames(match_res[,2], match_res[,3])
    names(x)[16] <- paste(names(x)[16], match_res[16, 4])
    for (i in c(14, 12, 8, 7)) {
      x <- append(x, setNames(match_res[i, 5], paste(names(x)[i], "percentage")), after = match(names(x)[i], names(x)))
    }
    return(x)
  }) |> t() |> as.data.frame()
) |> t() |>
  write.table(snakemake@output[[1]], sep="\t", col.names=FALSE)
