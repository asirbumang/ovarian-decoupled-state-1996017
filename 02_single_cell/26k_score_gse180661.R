# ============================================================
# 26k: GSE180661（Nature 2022, Vázquez-García；治疗前 160 位点/42 例）
# 子集：92,000 细胞 × 46 位点（每 site ≤2000 分层抽样，文档化）
# 逐细胞 AddModuleScore → decouple → 细胞类型 × 位点（卵巢 vs 网膜等）
# ============================================================
suppressMessages({library(Matrix); library(Seurat); library(data.table); library(ggplot2)
  library(ggpubr); library(patchwork); library(ggsci); library(ragg)})
set.seed(42)
subdir <- "D:/子宫内膜异位症2/data_sc/GSE180661/subset"
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

cat("loading subset...\n")
mtx <- readMM(gzfile(file.path(subdir, "matrix.mtx.gz"))); mtx <- as(mtx, "dgCMatrix")
fea <- fread(file.path(subdir, "features.tsv.gz"), header = FALSE, data.table = FALSE)
bar <- fread(file.path(subdir, "barcodes.tsv.gz"), header = FALSE, data.table = FALSE)
rownames(mtx) <- fea$V1; colnames(mtx) <- bar$V1
cat("matrix:", dim(mtx), "\n")

seu <- CreateSeuratObject(counts = mtx, min.cells = 0, min.features = 0)
seu <- NormalizeData(seu, verbose = FALSE)
si <- lapply(score_sets, function(g) intersect(g, rownames(seu)))
cat("coverage:", paste(sprintf("%s=%d", names(si), sapply(si, length)), collapse = " | "), "\n")
seu <- AddModuleScore(seu, features = si, name = "sc", nbin = 24, ctrl = 100)
cols <- grep("^sc[0-9]+$", colnames(seu@meta.data), value = TRUE)
df <- data.frame(cell = colnames(seu), stringsAsFactors = FALSE)
for (k in seq_along(score_sets)) df[[names(score_sets)[k]]] <- seu@meta.data[[cols[k]]]

meta <- fread(file.path(subdir, "cells_subset.csv"), data.table = FALSE)
meta$cell <- meta$cell_id
df <- merge(df, meta, by = "cell", all.x = TRUE)
df <- as.data.table(df)
df[, M1z := as.numeric(scale(M1_immune_activation))]
df[, M2z := as.numeric(scale(M2_metabolic_GSH))]
df[, M3z := as.numeric(scale(M3_mito_cytoskeleton))]
df[, decouple := M1z - (M2z + M3z) / 2]
fwrite(df, file.path(outdir, "scores_GSE180661_subset.csv"))

cat("\n== cell types ==\n")
print(table(df$cell_type, useNA = "ifany"))
cat("\n== tumor_supersite ==\n")
print(table(df$tumor_supersite, useNA = "ifany"))
cat("therapy:", paste(names(table(df$therapy)), collapse = " | "), "\n")

## 恶性细胞 × 位点（Cancer 细胞类型 或 cluster_label 含 Cancer）
cancer_ct <- unique(c(df$cell_type[grepl("Cancer|Epithelial|Malignant", df$cell_type)]))
cat("\ncancer cell types:", paste(cancer_ct, collapse = ", "), "\n")
canc <- df[df$cell_type %in% cancer_ct, ]
cat("cancer cells:", nrow(canc), "\n")
if (nrow(canc) > 100) {
  cat("cancer decouple by supersite:\n")
  s <- canc[, .(n = .N, med = median(decouple)), by = tumor_supersite]
  print(as.data.frame(s))
  ## 卵巢 vs 网膜
  sub2 <- canc[tumor_supersite %in% c("Ovary", "Omentum"), ]
  # 实际层级：tumor_site / tumor_subsite 更细
  w <- tryCatch(wilcox.test(decouple ~ tumor_site, data = as.data.frame(sub2)), error = function(e) NULL)
  if (!is.null(w)) cat("cancer decouple by exact tumor_site p =", format(w$p.value, digits = 4), "\n")
}
## 全细胞类型 × 解耦（描述）
cat("\ncell_type x decouple (median):\n")
s2 <- df[, .(n = .N, med = median(decouple)), by = cell_type]
print(as.data.frame(s2)[order(-s2$n), ], nrow = 20)
cat("DONE 26k\n")
