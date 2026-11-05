# ============================================================
# 26e-2: GSE189843 Visium 恶性 spot 子集 → ER vs PR（样本级）
# 恶性定义：上皮标记(EPCAM/KRT8/KRT18/PAX8/WFDC2) 均分 > P50（源内）
# ============================================================
suppressMessages({library(Matrix); library(Seurat); library(data.table); library(ggplot2)
  library(ggpubr); library(patchwork); library(ggsci); library(ragg); library(hdf5r)})
set.seed(42)
sd <- "D:/子宫内膜异位症2/data_sc/GSE189843/raw"
outdir <- "D:/子宫内膜异位症2/analysis/L0_bulk/results/spatial"

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
EPI_MARK <- intersect(c("EPCAM","KRT8","KRT18","KRT19","PAX8","WFDC2","MUC1"), NULL)  # 占位不报错
EPI_MARK <- c("EPCAM","KRT8","KRT18","KRT19","PAX8","WFDC2","MUC1")

grp_map <- data.frame(gsm = sprintf("GSM57084%02d", 85:96),
                      sample = c(paste0("ER", 1:6), paste0("PR", 1:6)), stringsAsFactors = FALSE)
files <- list.files(sd, pattern = "matrix_.*\\.mtx\\.gz$")
res_list <- list()
for (fi in files) {
  gsm <- sub("_matrix_.*", "", fi)
  grp <- grp_map$sample[match(gsm, grp_map$gsm)]
  mtx <- readMM(gzfile(file.path(sd, fi))); mtx <- as(mtx, "dgCMatrix")
  fea <- fread(file.path(sd, sub("\\.mtx\\.gz", ".tsv.gz", sub("matrix_", "features_", fi))), header = FALSE, data.table = FALSE)
  bar <- fread(file.path(sd, sub("\\.mtx\\.gz", ".tsv.gz", sub("matrix_", "barcodes_", fi))), header = FALSE, data.table = FALSE)
  rownames(mtx) <- fea$V2; colnames(mtx) <- bar$V1
  if (anyDuplicated(rownames(mtx))) mtx <- mtx[!duplicated(rownames(mtx)), ]
  nf <- Matrix::colSums(mtx > 0); nc <- Matrix::colSums(mtx)
  mtx <- mtx[, nf >= 100 & nc >= 200]
  seu <- CreateSeuratObject(counts = mtx, min.cells = 0, min.features = 0)
  seu <- NormalizeData(seu, verbose = FALSE)
  ## 上皮标记均分（logNorm）
  epi_g <- intersect(EPI_MARK, rownames(seu))
  epi_score <- colMeans(GetAssayData(seu, layer = "data")[epi_g, , drop = FALSE])
  si <- lapply(score_sets, function(g) intersect(g, rownames(seu)))
  seu <- AddModuleScore(seu, features = si, name = "sc", nbin = 24, ctrl = 100)
  cols <- grep("^sc[0-9]+$", colnames(seu@meta.data), value = TRUE)
  df <- data.frame(gsm = gsm, sample = grp, spot = colnames(seu), epi = epi_score, stringsAsFactors = FALSE)
  for (k in seq_along(score_sets)) df[[names(score_sets)[k]]] <- seu@meta.data[[cols[k]]]
  res_list[[fi]] <- df
  rm(seu, mtx); gc()
}
all <- rbindlist(res_list, fill = TRUE)
cat("spots:", nrow(all), "\n")
## 全数据标准化 + 恶性 spot（上皮分数 P50 内源内→ 全样本统一 P50）
all[, M1z := as.numeric(scale(M1_immune_activation))]
all[, M2z := as.numeric(scale(M2_metabolic_GSH))]
all[, M3z := as.numeric(scale(M3_mito_cytoskeleton))]
all[, decouple := M1z - (M2z + M3z) / 2]
p50_epi <- quantile(all$epi, 0.5, na.rm = TRUE)
mal <- all[epi >= p50_epi, ]
cat("malignant spots:", nrow(mal), "\n")
samp <- mal[, .(n_spots = .N, med = median(decouple)), by = sample]
samp[, grp := ifelse(grepl("^ER", sample), "ER", "PR")]
cat("sample-level (malignant spots):\n"); print(as.data.frame(samp))
w <- tryCatch(wilcox.test(med ~ grp, data = as.data.frame(samp)), error = function(e) NULL)
cat("ER vs PR (malignant spots) p =", if (!is.null(w)) format(w$p.value, digits = 4) else NA, "\n")
cat("ER med:", round(median(samp$med[samp$grp == "ER"]), 3), " PR med:", round(median(samp$med[samp$grp == "PR"]), 3), "\n")
fwrite(as.data.frame(mal), file.path(outdir, "GSE189843_malignant_spot_scores.csv"))
fwrite(as.data.frame(samp), file.path(outdir, "GSE189843_malignant_sample_summary.csv"))
cat("DONE 26e-2\n")
