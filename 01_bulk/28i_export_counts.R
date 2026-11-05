
suppressMessages({library(Matrix); library(data.table)})
args <- commandArgs(trailingOnly = TRUE)
GSE <- args[1]
baseraw <- file.path("D:/子宫内膜异位症2/data_spatial", GSE, "raw")
outdir <- "D:/子宫内膜异位症2/analysis/L0_bulk/results/spatial/deconv_input"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
files <- list.files(baseraw, pattern = "\\.mtx\\.gz$")
for (fi in files) {
  feat_path <- file.path(baseraw, sub("\\.mtx\\.gz$", ".tsv.gz", sub("_matrix", "_features", fi)))
  if (!file.exists(feat_path)) next
  code <- sub("^GSM[0-9]+_", "", sub("_matrix\\.mtx\\.gz$", "", fi))
  mtx <- readMM(gzfile(file.path(baseraw, fi))); mtx <- as(mtx, "dgCMatrix")
  fea <- fread(feat_path, header = FALSE, data.table = FALSE)
  bar <- fread(file.path(baseraw, sub("\\.mtx\\.gz$", ".tsv.gz", sub("_matrix", "_barcodes", fi))), header = FALSE, data.table = FALSE)
  rownames(mtx) <- fea$V2; colnames(mtx) <- bar$V1
  s <- summary(mtx)
  dt <- data.table(g = rownames(mtx)[s$i], c = colnames(mtx)[s$j], v = s$x)
  fwrite(dt, file.path(outdir, paste0(GSE, "_", code, ".triplet.gz")), row.names = FALSE)
  fwrite(data.table(gene = rownames(mtx)), file.path(outdir, paste0(GSE, "_", code, ".genes.tsv")))
  fwrite(data.table(spot = colnames(mtx)), file.path(outdir, paste0(GSE, "_", code, ".spots.tsv")))
  cat("  exported", code, "nnz:", nrow(dt), "\n")
  rm(mtx, s, dt); gc()
}
cat("DONE 28i", GSE, "\n")
