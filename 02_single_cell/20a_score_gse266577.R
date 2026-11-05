# ============================================================
# scRNA P0 · GSE266577（Cancer Cell 2024；137,083 细胞；化疗前 vs IDS）
# 母版 27a 口径：counts → Seurat → NormalizeData → AddModuleScore(nbin=24, ctrl=100)
# 输出：analysis/L1_scRNA/results/module_scores/scores_GSE266577.csv
# ============================================================
suppressMessages({
  library(Matrix); library(Seurat); library(data.table)
})
set.seed(42)
outdir <- "D:/子宫内膜异位症2/analysis/L1_scRNA"
sdir <- file.path(outdir, "results", "module_scores")
dir.create(sdir, recursive = TRUE, showWarnings = FALSE)

sets <- readRDS("D:/子宫内膜异位症2/pipeline/mitoxy_sets.RDS")
modules <- sets$modules; icc45 <- sets$icc45
mitoxy_genes <- unique(unlist(modules))
death_sets_clean <- list(
  pyroptosis  = setdiff(c("GSDMD","GSDMC","NLRP3","CASP1","CASP4","CASP5"), mitoxy_genes),
  necroptosis = setdiff(c("MLKL","RIPK3","RIPK1"), mitoxy_genes),
  apoptosis   = setdiff(c("CASP3","CASP7","CASP8","CASP9","BCL2","BAX"), mitoxy_genes),
  ferroptosis = setdiff(c("TFRC","NCOA4","FTH1","FTL","HMOX1","CHAC1","PTGS2","SLC3A2"), mitoxy_genes))
score_sets <- c(modules, list(ICC45 = icc45),
                list(pyroptosis = death_sets_clean$pyroptosis,
                     necroptosis = death_sets_clean$necroptosis,
                     apoptosis = death_sets_clean$apoptosis,
                     ferroptosis = death_sets_clean$ferroptosis))

## 1. 读稀疏矩阵（58,037 × 137,083）
cat("reading matrix...\n")
t0 <- Sys.time()
counts <- readMM(gzfile("D:/子宫内膜异位症2/data_sc/GSE266577/GSE266577_counts_raw.mtx.gz"))
counts <- as(counts, "dgCMatrix")
feat <- fread("D:/子宫内膜异位症2/data_sc/GSE266577/GSE266577_seurat_features.txt.gz",
              header = FALSE, data.table = FALSE)
bar <- fread("D:/子宫内膜异位症2/data_sc/GSE266577/GSE266577_barcodes.txt.gz",
             header = FALSE, data.table = FALSE)
rownames(counts) <- feat$V1
colnames(counts) <- bar$V1
cat(sprintf("matrix %d x %d, %.1f min\n", nrow(counts), ncol(counts),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

## 2. 元数据 + 自定 QC（miqc.keep 全 FALSE，不能用）
meta <- fread("D:/子宫内膜异位症2/data_sc/GSE266577/GSE266577_metadata.txt.gz", data.table = FALSE)
rownames(meta) <- meta$cell_name
keep <- meta[colnames(counts), ]
keep$nFeature_RNA <- as.numeric(keep$nFeature_RNA)
keep$percent.mt <- as.numeric(keep$percent.mt)
keep$nCount_RNA <- as.numeric(keep$nCount_RNA)
qc <- keep$nFeature_RNA >= 200 & keep$percent.mt <= 20
cat("QC pass:", sum(qc), "/", nrow(keep), "\n")
counts <- counts[, qc]
meta2 <- keep[qc, ]

## 3. Seurat：Single 层评分对象（母版 27a 逐字）
seu <- CreateSeuratObject(counts = counts, min.cells = 0, min.features = 0)
seu <- NormalizeData(seu, verbose = FALSE)
gn <- rownames(seu)
score_sets_i <- lapply(score_sets, function(g) intersect(g, gn))
cat("gene coverage:\n"); print(sapply(score_sets_i, length))
miss <- unique(unlist(lapply(score_sets, function(g) setdiff(g, gn))))
if (length(miss)) cat("missing genes:", paste(miss, collapse = ","), "\n")

seu <- AddModuleScore(seu, features = score_sets_i, name = "sc", nbin = 24, ctrl = 100)
sc_cols <- grep("^sc[0-9]+$", colnames(seu@meta.data), value = TRUE)
stopifnot(length(sc_cols) == length(score_sets))

## 4. 落盘（cell, dataset, 9 分, 关键 meta）
out <- data.frame(cell = colnames(seu), dataset = "GSE266577", stringsAsFactors = FALSE)
out <- cbind(out, seu@meta.data[, sc_cols, drop = FALSE])
names(out) <- c("cell", "dataset", names(score_sets))
keep_cols <- intersect(c("cell_type", "treatment_stage", "publication_patient_code_final",
                         "publication_sample_code_final", "nCount_RNA", "nFeature_RNA", "percent.mt"),
                       colnames(meta2))
out <- cbind(out, meta2[colnames(seu), keep_cols, drop = FALSE])
fwrite(out, file.path(sdir, "scores_GSE266577.csv"))
cat("saved scores_GSE266577.csv:", nrow(out), "cells\n")

## 5. 摘要统计（先用官方注释做快检；注释审核按协议单列）
cat("\n== quick summary (official labels) ==\n")
## 逐细胞解耦分数（同 track0 定义：z(M1) - (z(M2)+z(M3))/2，跨细胞标准化）
out$M1z <- as.numeric(scale(out$M1_immune_activation))
out$M2z <- as.numeric(scale(out$M2_metabolic_GSH))
out$M3z <- as.numeric(scale(out$M3_mito_cytoskeleton))
out$decouple <- out$M1z - (out$M2z + out$M3z) / 2
out$state_high <- out$M1_immune_activation >= quantile(out$M1_immune_activation, .75) &
                  out$M2_metabolic_GSH <= quantile(out$M2_metabolic_GSH, .25) &
                  out$M3_mito_cytoskeleton <= quantile(out$M3_mito_cytoskeleton, .25)
fwrite(out, file.path(sdir, "scores_GSE266577.csv"))

## 恶性（上皮）细胞 × 治疗阶段
epi <- out[out$cell_type == "Epithelial cells", ]
cat("epithelial cells:", nrow(epi), " states:", sum(epi$state_high), "\n")
for (st in c("chemo-naive", "IDS")) {
  s <- epi$decouple[epi$treatment_stage == st]
  cat(sprintf("%-12s n=%5d  median decouple=%.3f  state_high=%.2f%%\n",
              st, length(s), median(s), 100 * mean(epi$state_high[epi$treatment_stage == st])))
}
w <- tryCatch(wilcox.test(decouple ~ treatment_stage, data = epi), error = function(e) NULL)
if (!is.null(w)) cat("Wilcoxon decouple chemo-naive vs IDS (epithelial) p =", format(w$p.value, digits = 4), "\n")

## 各细胞类型 decouple 中位（化疗前 vs IDS）
cat("\ncell-type x stage median decouple:\n")
summ <- aggregate(decouple ~ cell_type + treatment_stage, data = out, FUN = median)
print(as.data.frame(summ)[order(summ$cell_type), ])
cat("\nDONE GSE266577 scoring\n")
