#!/usr/bin/env Rscript
# RMA normalize GSE22762 HGU133A (44) + HGU133B (44), extract CD200/CD70/TNFSF9/CD86.
suppressMessages({
  library(affy)
})

raw_dir <- "/home/zhh2/CLL-ScRNA-seqdatasets/extracted/GSE22762_RAW"
out <- "/home/zhh2/CLL-ScRNA-seqdatasets/NK_analysis/data"
GENES <- c("CD200","CD70","TNFSF9","CD86")

process_platform <- function(suffix, db_name, outfile) {
  suppressMessages(library(db_name, character.only=TRUE))
  cels <- list.files(raw_dir, pattern=paste0("_", suffix, "\\.CEL\\.gz$"), full.names=TRUE)
  cat("\n=== ", suffix, ": ", length(cels), " CEL files ===\n", sep="")
  eset <- rma(ReadAffy(filenames=cels))
  expr <- exprs(eset)
  colnames(expr) <- sub(paste0("_", suffix, "\\.CEL\\.gz$"), "", basename(cels))
  # probe -> symbol
  map <- select(get(db_name), keys=GENES, keytype="SYMBOL", columns="PROBEID")
  res <- list()
  for (g in GENES) {
    ps <- unique(map$PROBEID[map$SYMBOL==g]); ps <- ps[ps %in% rownames(expr)]
    if (length(ps) > 0) { res[[g]] <- colMeans(expr[ps,,drop=FALSE]); cat(g,": ",length(ps)," probes\n",sep="") }
    else cat(g,": no probes\n")
  }
  df <- as.data.frame(res); df$sample <- rownames(df)
  write.table(df, file.path(out, outfile), sep="\t", row.names=FALSE, quote=FALSE)
  cat("saved", outfile, "\n")
}

process_platform("A", "hgu133a.db", "GSE22762_expr_A.tsv")
process_platform("B", "hgu133b.db", "GSE22762_expr_B.tsv")
cat("DONE\n")
