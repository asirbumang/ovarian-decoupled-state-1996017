# ============================================================
# 26i: GSE222556（Cell 2024 Luo 等；新辅助 PARPi vs 化疗，HRD HGSC）
# 67 个 10x 样本（201 文件）→ 逐样本 Seurat 评分 → decouple → 患者/样本汇总
# 注：样本码 P#_M#-# / P#_P#-# 的臂含义待论文表 S1 解码（已知 M=单药/P=铂类为候选）。
# ============================================================
suppressMessages({
  library(Matrix); library(Seurat); library(data.table); library(ggplot2)
  library(ggpubr); library(patchwork); library(ggsci); library(ragg)
})
set.seed(42)
rawdir <- "D:/子宫内膜异位症2/data_sc/GSE222556/raw"
outdir <- "D:/子宫内膜异位症2/analysis/L1_scRNA/results/module_scores"
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

## 细胞注释（官方；用于分层）
annot <- fread("D:/子宫内膜异位症2/data_sc/GSE222556/GSE222556_10x_cell_annotation.txt.gz", data.table = FALSE)
cat("annot cells:", nrow(annot), "\n")
## 样本 → GSM
files <- list.files(rawdir, pattern = "\\.mtx\\.gz$")
cat("samples:", length(files), "\n")

res_list <- list()
ensg_map <- fread("D:/子宫内膜异位症2/data_omics/ensg2symbol_gse222556.tsv", data.table = FALSE)
for (fi in files) {
  gsm <- sub("_matrix.*$", "", fi)
  sample_name <- sub("^GSM[0-9]+_", "", gsm)   # 如 P01_P1-1
  feat_f <- sub("_matrix\\.mtx\\.gz$", "_features.tsv.gz", fi)
  bar_f <- sub("_matrix\\.mtx\\.gz$", "_barcodes.tsv.gz", fi)
  mtx <- readMM(gzfile(file.path(rawdir, fi))); mtx <- as(mtx, "dgCMatrix")
  fea <- fread(file.path(rawdir, feat_f), header = FALSE, data.table = FALSE)
  bar <- fread(file.path(rawdir, bar_f), header = FALSE, data.table = FALSE)
  ## ENSG → symbol（mygene 映射）
  idx <- match(fea$V2, ensg_map$ENSG)
  syms <- ensg_map$Symbol[idx]
  rownames(mtx) <- ifelse(is.na(syms) | syms == "", fea$V2, syms)
  colnames(mtx) <- bar$V1
  if (anyDuplicated(rownames(mtx))) mtx <- mtx[!duplicated(rownames(mtx)), ]
  nf <- Matrix::colSums(mtx > 0); nc <- Matrix::colSums(mtx)
  keep <- nf >= 100 & nc >= 200
  mtx <- mtx[, keep]
  seu <- CreateSeuratObject(counts = mtx, min.cells = 0, min.features = 0)
  seu <- NormalizeData(seu, verbose = FALSE)
  si <- lapply(score_sets, function(g) intersect(g, rownames(seu)))
  seu <- AddModuleScore(seu, features = si, name = "sc", nbin = 24, ctrl = 100)
  cols <- grep("^sc[0-9]+$", colnames(seu@meta.data), value = TRUE)
  stopifnot(length(cols) == length(score_sets))
  df <- data.frame(sample = sample_name, cell = colnames(seu), stringsAsFactors = FALSE)
  for (k in seq_along(score_sets)) df[[names(score_sets)[k]]] <- seu@meta.data[[cols[k]]]
  res_list[[fi]] <- df
  cat("  ", sample_name, "cells:", ncol(seu), "\n", sep = "")
  rm(seu, mtx); gc()
}
all <- rbindlist(res_list, fill = TRUE)
cat("total cells:", nrow(all), "\n")

## 细胞注释挂接（barcode 前缀与 annot$Barcode 一致，如 ...-P20_M1-1）
annot$key <- annot$Barcode
all$key <- ifelse(grepl("-P", all$cell), sub("^(.*-P[0-9]+_.*)$", "\\1", all$cell), all$cell)
## 直接末段匹配：注释 Barcode 形如 "AAACCCAAGTTGCTGT-P20_M1-1"；细胞名同格式
matched <- annot$key %in% all$cell | all$cell %in% annot$key
cat("annot-sample overlap check:", sum(annot$key %in% all$cell), "/", nrow(annot), "\n")
all <- merge(all, annot[, c("key", "Major cell type", "Minor cell state")],
             by.x = "cell", by.y = "key", all.x = TRUE)
cat("matched annotations:", sum(!is.na(all$`Major cell type`)), "/", nrow(all), "\n")

## decouple + state（全部细胞合并标准化）
mkz <- function(x) as.numeric(scale(x))
all[, M1z := mkz(M1_immune_activation)]
all[, M2z := mkz(M2_metabolic_GSH)]
all[, M3z := mkz(M3_mito_cytoskeleton)]
all[, decouple := M1z - (M2z + M3z) / 2]
fwrite(all, file.path(outdir, "scores_GSE222556.csv"))

## 样本级汇总（患者/样本码）
samp <- all[, .(n_cells = .N, med_dec = median(decouple),
                pct_high = 100 * mean(decouple > quantile(decouple, 0.75))), by = sample]
samp[, patient := sub("_", "", sub("_.*", "", sample))]
samp[, code := sub("^[A-Z0-9]+_", "", sample)]
samp[, arm_tag := sub("^([A-Z]).*", "\\1", code)]
fwrite(samp, file.path(outdir, "GSE222556_sample_summary.csv"))
cat("\n== sample summary ==\n")
print(samp[, .(sample, n_cells, med_dec)], nrow = 30)

## 快速分组检验（arm_tag 为候选解码：P vs M；患者级比较需先解码——仅描述）
cat("\narm_tag medians: P =", round(median(samp$med_dec[samp$arm_tag == "P"]), 3),
    " M =", round(median(samp$med_dec[samp$arm_tag == "M"]), 3), "\n")
cat("cells P:", sum(samp$n_cells[samp$arm_tag == "P"]), " M:", sum(samp$n_cells[samp$arm_tag == "M"]), "\n")
cat("DONE 26i\n")
