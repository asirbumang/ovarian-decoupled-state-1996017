# ============================================================
# 28a: 空间集通用评分（26e 同口径；arg1=数据集名如 GSE274657）
# 每样本: counts → Seurat QC → NormalizeData → AddModuleScore(9) → decouple(数据集内 z)
# 输出: results/spatial/<GSE>_spot_scores.csv
# ============================================================
suppressMessages({library(Matrix); library(Seurat); library(data.table); library(hdf5r)})
set.seed(42)
GSE <- commandArgs(trailingOnly = TRUE)[1]
stopifnot(!is.na(GSE))
baseraw <- file.path("D:/子宫内膜异位症2/data_spatial", GSE, "raw")
outdir <- "D:/子宫内膜异位症2/analysis/L0_bulk/results/spatial"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

sets <- readRDS("D:/子宫内膜异位症2/pipeline/mitoxy_sets.RDS")
modules <- sets$modules; icc45 <- sets$icc45
mitoxy_genes <- unique(unlist(modules))
death_sets_clean <- list(
  pyroptosis  = setdiff(c("GSDMD","GSDMC","NLRP3","CASP1","CASP4","CASP5"), mitoxy_genes),
  necroptosis = setdiff(c("MLKL","RIPK3","RIPK1"), mitoxy_genes),
  apoptosis   = setdiff(c("CASP3","CASP7","CASP8","CASP9","BCL2","BAX"), mitoxy_genes),
  ferroptosis = setdiff(c("TFRC","NCOA4","FTH1","FTL","HMOX1","CHAC1","PTGS2","SLC3A2"), mitoxy_genes))
score_sets <- c(modules, list(ICC45 = icc45), list(pyroptosis = death_sets_clean$pyroptosis,
                 necroptosis = death_sets_clean$necroptosis,
                 apoptosis = death_sets_clean$apoptosis,
                 ferroptosis = death_sets_clean$ferroptosis))

# 样本文件模式：GSM###_<sample>_matrix.mtx.gz
files <- list.files(baseraw, pattern = "\\.mtx\\.gz$")
cat(GSE, "samples:", length(files), "\n")
res_list <- list()
for (fi in files) {
  samp_code <- sub("^GSM[0-9]+_", "", sub("_matrix\\.mtx\\.gz$", "", fi))
  cat("  ", samp_code, "\n", sep = "")
  mtx <- readMM(gzfile(file.path(baseraw, fi))); mtx <- as(mtx, "dgCMatrix")
  fea <- fread(file.path(baseraw, sub("\\.mtx\\.gz", ".tsv.gz", sub("_matrix", "_features", fi))), header = FALSE, data.table = FALSE)
  bar <- fread(file.path(baseraw, sub("\\.mtx\\.gz", ".tsv.gz", sub("_matrix", "_barcodes", fi))), header = FALSE, data.table = FALSE)
  rownames(mtx) <- make.unique(fea$V2); colnames(mtx) <- bar$V1
  if (anyDuplicated(rownames(mtx))) mtx <- mtx[!duplicated(rownames(mtx)), ]
  nf <- Matrix::colSums(mtx > 0); nc <- Matrix::colSums(mtx)
  mtx <- mtx[, nf >= 100 & nc >= 200]
  seu <- CreateSeuratObject(counts = mtx, min.cells = 0, min.features = 0)
  seu <- NormalizeData(seu, verbose = FALSE)
  si <- lapply(score_sets, function(g) intersect(g, rownames(seu)))
  seu <- AddModuleScore(seu, features = si, name = "sc", nbin = 24, ctrl = 100)
  cols <- grep("^sc[0-9]+$", colnames(seu@meta.data), value = TRUE)
  stopifnot(length(cols) == length(score_sets))
  df <- data.frame(sample = samp_code, spot = colnames(seu), stringsAsFactors = FALSE)
  for (k in seq_along(score_sets)) df[[names(score_sets)[k]]] <- seu@meta.data[[cols[k]]]
  ## 坐标（如存在）
  pos_file <- list.files(baseraw, pattern = paste0("^", sub("_matrix.*", "", fi), ".*tissue_positions.*\\.csv\\.gz$"), full.names = TRUE)
  if (length(pos_file)) {
    pos <- fread(pos_file[1], data.table = FALSE)
    lab <- c("barcode","in_tissue","array_row","array_col","pxl_row","pxl_col")
    idx <- match(colnames(pos)[1:6], lab)
    names(pos)[1:6] <- lab
    pos$barcode <- as.character(pos$barcode)
    map <- match(df$spot, pos$barcode)
    df$array_row <- pos$array_row[map]; df$array_col <- pos$array_col[map]
    df$pxl_row <- pos$pxl_row[map]; df$pxl_col <- pos$pxl_col[map]
  }
  res_list[[fi]] <- df
  rm(seu, mtx); gc()
}
all <- rbindlist(res_list, fill = TRUE)
all[, M1z := as.numeric(scale(M1_immune_activation))]
all[, M2z := as.numeric(scale(M2_metabolic_GSH))]
all[, M3z := as.numeric(scale(M3_mito_cytoskeleton))]
all[, decouple := M1z - (M2z + M3z) / 2]
all[, state_high := M1_immune_activation >= quantile(M1_immune_activation, .75) &
     M2_metabolic_GSH <= quantile(M2_metabolic_GSH, .25) &
     M3_mito_cytoskeleton <= quantile(M3_mito_cytoskeleton, .25)]
fwrite(all, file.path(outdir, paste0(GSE, "_spot_scores.csv")))
cat("saved:", file.path(outdir, paste0(GSE, "_spot_scores.csv")), " rows:", nrow(all), "\n")
cat("state_high:", sum(all$state_high), "\n")
cat("DONE 28a", GSE, "\n")
