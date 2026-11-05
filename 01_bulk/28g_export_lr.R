
suppressMessages({library(Matrix); library(data.table)})
args <- commandArgs(trailingOnly = TRUE)
GSE <- args[1]
baseraw <- file.path("D:/子宫内膜异位症2/data_spatial", GSE, "raw")
outdir <- "D:/子宫内膜异位症2/analysis/L0_bulk/results/spatial/lr_input"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
LR_GENES <- c("FASLG","FAS","TNF","TNFRSF1A","TNFRSF1B","TNFSF10","TNFRSF10B","IL1B","IL1R1","IL18","IL18R1",
  "IFNG","IFNGR1","IFNGR2","TGFB1","TGFBR1","TGFBR2","CSF1","CSF1R","SPP1","CD44","MDK","SDC2","SDC1",
  "CXCL12","CXCR4","CXCL16","CXCR6","CCL2","CCR2","CCL5","CCR5","CD274","PDCD1","NECTIN2","TIGIT",
  "HGF","MET","EGF","EGFR","VEGFA","KDR","IL6","IL6R","IL10","IL10RA","LGALS9","HAVCR2","IDO1","BSG")
files <- list.files(baseraw, pattern = "\\.mtx\\.gz$")
for (fi in files) {
  feat_path <- file.path(baseraw, sub("\\.mtx\\.gz$", ".tsv.gz", sub("_matrix", "_features", fi)))
  if (!file.exists(feat_path)) next
  code <- sub("^GSM[0-9]+_", "", sub("_matrix\\.mtx\\.gz$", "", fi))
  mtx <- readMM(gzfile(file.path(baseraw, fi))); mtx <- as(mtx, "dgCMatrix")
  fea <- fread(feat_path, header = FALSE, data.table = FALSE)
  rownames(mtx) <- make.unique(fea$V2)
  mtx <- mtx[!duplicated(rownames(mtx)), ]
  keep <- intersect(LR_GENES, rownames(mtx))
  cpm <- sweep(mtx[keep, ], 2, Matrix::colSums(mtx) / 1e6, "/")
  out <- as.matrix(log1p(cpm))
  bar <- fread(file.path(baseraw, sub("\\.mtx\\.gz$", ".tsv.gz", sub("_matrix", "_barcodes", fi))), header = FALSE, data.table = FALSE)
  colnames(out) <- bar$V1
  gz <- gzfile(file.path(outdir, paste0(GSE, "_", code, ".csv.gz")), "wt")
  write.csv(out, gz); close(gz)
  cat("  exported", code, length(keep), "genes", "\n")
}
cat("DONE 28g", GSE, "\n")
