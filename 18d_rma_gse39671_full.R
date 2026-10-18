#!/usr/bin/env Rscript
# Step 18d — GSE39671 RMA, save FULL matrix + core genes + positive-control
# prognostic markers (LPL, CD38, ZAP70 with individual probes) for pipeline QC.
suppressMessages({
  library(affy)
  library(hgu133plus2.db)
  library(AnnotationDbi)
})

raw_dir <- "/home/zhh2/CLL-ScRNA-seqdatasets/GSE39671_RAW"
out <- "/home/zhh2/CLL-ScRNA-seqdatasets/NK_analysis/data"

cels <- list.files(raw_dir, pattern = "\\.CEL\\.gz$", full.names = TRUE)
affy <- ReadAffy(filenames = cels)
eset <- rma(affy)
expr <- exprs(eset)  # probes x samples
sids <- sub("\\.CEL\\.gz$", "", basename(cels))
sids <- sub("_.*$", "", sids)
colnames(expr) <- sids

# save full matrix (probes x samples)
write.table(expr, file.path(out, "GSE39671_expr_full_probes.tsv"),
            sep = "\t", quote = FALSE)
cat("saved full matrix:", nrow(expr), "probes x", ncol(expr), "samples\n")

# gene-averaged (core + controls)
GENES <- c("CD200", "CD70", "TNFSF9", "CD86", "ZAP70", "LPL", "CD38")
probes <- select(hgu133plus2.db, keys = GENES, keytype = "SYMBOL", columns = "PROBEID")
res <- list()
for (g in GENES) {
  ps <- unique(probes$PROBEID[probes$SYMBOL == g])
  ps <- ps[ps %in% rownames(expr)]
  res[[g]] <- colMeans(expr[ps, , drop = FALSE])
}
df <- as.data.frame(res); df$gsm <- rownames(df)
write.table(df[, c("gsm", GENES)], file.path(out, "GSE39671_expr.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("saved GSE39671_expr.tsv (gene-averaged)\n")

# individual ZAP70 probes for the bimodal-check
for (p in c("214032_at", "1555613_a_at")) {
  if (p %in% rownames(expr)) {
    v <- expr[p, ]
    cat(sprintf("ZAP70 %s: mean %.3f, sd %.3f, >20%%tile frac high\n", p, mean(v), sd(v)))
  }
}
# LPL probes
lplp <- unique(probes$PROBEID[probes$SYMBOL == "LPL"])
cat("LPL probes:", paste(lplp[lplp %in% rownames(expr)], collapse=", "), "\n")
