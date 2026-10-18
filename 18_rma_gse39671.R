#!/usr/bin/env Rscript
# Step 18 — GSE39671 (UCSD CLL, 130 samples, GPL570/HGU133plus2) RMA + gene extraction.
# External replication cohort #2 for the immune-evasion ligand survival finding.
# Endpoint: time-to-first-treatment (TTFT), NOT OS. (IGHV not in GEO annotation.)

suppressMessages({
  library(affy)
  library(hgu133plus2.db)
  library(AnnotationDbi)
})

raw_dir <- "/home/zhh2/CLL-ScRNA-seqdatasets/GSE39671_RAW"
out <- "/home/zhh2/CLL-ScRNA-seqdatasets/NK_analysis/data"

cels <- list.files(raw_dir, pattern = "\\.CEL\\.gz$", full.names = TRUE)
cat("processing", length(cels), "CEL files\n")

affy <- ReadAffy(filenames = cels)
eset <- rma(affy)
expr <- exprs(eset)  # probes x samples

# GSM accession from CEL filename (GSM977137_SDG_00021.CEL.gz -> GSM977137)
sids <- sub("\\.CEL\\.gz$", "", basename(cels))
sids <- sub("_.*$", "", sids)
colnames(expr) <- sids

GENES <- c("CD200", "CD70", "TNFSF9", "CD86", "ZAP70")
probes <- select(hgu133plus2.db, keys = GENES, keytype = "SYMBOL", columns = "PROBEID")
cat("probes found:\n"); print(probes)

res <- list()
for (g in GENES) {
  ps <- unique(probes$PROBEID[probes$SYMBOL == g])
  ps <- ps[ps %in% rownames(expr)]
  if (length(ps) > 0) {
    res[[g]] <- colMeans(expr[ps, , drop = FALSE])
    cat(sprintf("%s: %d probes, mean expr %.3f\n", g, length(ps), mean(res[[g]])))
  } else {
    cat(g, ": no probes\n")
  }
}
df <- as.data.frame(res)
df$gsm <- rownames(df)
df <- df[, c("gsm", GENES)]
write.table(df, file.path(out, "GSE39671_expr.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
cat("saved GSE39671_expr.tsv\n")
