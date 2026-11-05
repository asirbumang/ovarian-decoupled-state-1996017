
suppressMessages({library(Matrix); library(data.table)})
args <- commandArgs(trailingOnly = TRUE)
GSE <- args[1]
baseraw <- file.path("D:/子宫内膜异位症2/data_spatial", GSE, "raw")
outdir <- "D:/子宫内膜异位症2/analysis/L0_bulk/results/spatial/nmf_input"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
files <- list.files(baseraw, pattern = "\\.mtx\\.gz$")
for (fi in files) {
  feat_path <- file.path(baseraw, sub("\\.mtx\\.gz$", ".tsv.gz", sub("_matrix", "_features", fi)))
  if (!file.exists(feat_path)) next
  code <- sub("^GSM[0-9]+_", "", sub("_matrix\\.mtx\\.gz$", "", fi))
  mtx <- readMM(gzfile(file.path(baseraw, fi))); mtx <- as(mtx, "dgCMatrix")
  fea <- fread(feat_path, header = FALSE, data.table = FALSE)
  if (file.exists(file.path(baseraw, sub("\\.mtx\\.gz$", ".tsv.gz", sub("_matrix", "_barcodes", fi))))) {
    bar <- fread(file.path(baseraw, sub("\\.mtx\\.gz$", ".tsv.gz", sub("_matrix", "_barcodes", fi))), header = FALSE, data.table = FALSE)
    colnames(mtx) <- bar$V1
  }
  rownames(mtx) <- make.unique(fea$V2)
  mtx <- mtx[!duplicated(rownames(mtx)), ]
  cpm <- sweep(mtx, 2, Matrix::colSums(mtx) / 1e6, "/")
  logc <- log2(cpm + 1)
  v <- apply(logc, 1, var)
  top <- names(sort(v, decreasing = TRUE))[1:min(2000, length(v))]
  out <- as.matrix(logc[top, ])
  gz <- gzfile(file.path(outdir, paste0(GSE, "_", code, ".csv.gz")), "wt")
  write.csv(out, gz)
  close(gz)
  cat("  exported", code, nrow(out), "x", ncol(out), "\n")
  rm(mtx, cpm, logc, out); gc()
}
cat("DONE 28c", GSE, "\n")
