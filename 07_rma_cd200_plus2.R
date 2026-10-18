#!/usr/bin/env Rscript
# RMA normalize GSE22762 CEL files, extract CD200 (and CD70/TNFSF9) expression.
# Uses HGU133plus2 (GPL570, 107 samples) — the majority platform with CD200.

suppressMessages({
  library(affy)
  library(hgu133plus2.db)
})

raw_dir <- "/home/zhh2/CLL-ScRNA-seqdatasets/extracted/GSE22762_RAW"
out <- "/home/zhh2/CLL-ScRNA-seqdatasets/NK_analysis/data"

# plus2 CEL files (no _A/_B suffix)
cels <- list.files(raw_dir, pattern="\\.CEL\\.gz$", full.names=TRUE)
cels <- cels[!grepl("_[AB]\\.CEL", cels)]
cat("processing", length(cels), "HGU133plus2 CEL files\n")

affy <- ReadAffy(filenames=cels)
eset <- rma(affy)
expr <- exprs(eset)   # probes x samples

# sample IDs from CEL filenames (GSMxxxx)
sids <- sub("\\.CEL\\.gz$", "", basename(cels))
colnames(expr) <- sids

# CD200 probes
probes <- select(hgu133plus2.db, keys=c("CD200","CD70","TNFSF9","CD86"),
                 keytype="SYMBOL", columns="PROBEID")
cat("probes found:\n"); print(probes)

res <- list()
for (g in c("CD200","CD70","TNFSF9","CD86")) {
  ps <- unique(probes$PROBEID[probes$SYMBOL==g])
  ps <- ps[ps %in% rownames(expr)]
  if (length(ps) > 0) {
    v <- colMeans(expr[ps, , drop=FALSE])   # average over probes
    res[[g]] <- v
    cat(sprintf("%s: %d probes, mean expr %.2f\n", g, length(ps), mean(v)))
  } else {
    cat(g, ": no probes\n")
  }
}
df <- as.data.frame(res)
df$sample <- rownames(df)
write.table(df, file.path(out, "GSE22762_expr_plus2.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
cat("saved GSE22762_expr_plus2.tsv\n")
