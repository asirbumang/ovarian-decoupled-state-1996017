# ============================================================
# 26a: 卵巢 scRNA（GSE266577）→ 6 谱系伪bulk 参考（母版 28b 参考同构）
# 谱系映射（官方 22 类 → 6）：Epithelial / Stromal(Fib,Meso) / Endothelial /
#   Myeloid(Mφ,Mono,DC,Mast) / T_NK(T,NK,ILC,CD16±) / B_Plasma(B,Plasma)
# ============================================================
suppressMessages({
  library(Matrix); library(Seurat); library(data.table)
})
set.seed(20260831)
oc <- "D:/子宫内膜异位症2/analysis/L1_scRNA"
dir.create(file.path(oc, "results", "track2"), recursive = TRUE, showWarnings = FALSE)

## 1. 载入计数 + 元数据（与 20a 完全一致的 QC）
counts <- readMM(gzfile("D:/子宫内膜异位症2/data_sc/GSE266577/GSE266577_counts_raw.mtx.gz"))
counts <- as(counts, "dgCMatrix")
feat <- fread("D:/子宫内膜异位症2/data_sc/GSE266577/GSE266577_seurat_features.txt.gz", header = FALSE, data.table = FALSE)
bar <- fread("D:/子宫内膜异位症2/data_sc/GSE266577/GSE266577_barcodes.txt.gz", header = FALSE, data.table = FALSE)
rownames(counts) <- feat$V1; colnames(counts) <- bar$V1
meta <- fread("D:/子宫内膜异位症2/data_sc/GSE266577/GSE266577_metadata.txt.gz", data.table = FALSE)
rownames(meta) <- meta$cell_name
keep <- meta[colnames(counts), ]
qc <- as.numeric(keep$nFeature_RNA) >= 200 & as.numeric(keep$percent.mt) <= 20
counts <- counts[, qc]; meta2 <- keep[qc, ]
cat("cells:", ncol(counts), "\n")

## 2. 谱系映射
map_lin <- c(
  "Epithelial cells" = "Epithelial",
  "Fibroblasts" = "Stromal", "Mesothelial cells" = "Stromal",
  "Endothelial cells" = "Endothelial",
  "Macrophages" = "Myeloid", "Classical monocytes" = "Myeloid",
  "Migratory DCs" = "Myeloid", "pDC" = "Myeloid", "DC1" = "Myeloid", "DC2" = "Myeloid",
  "Mast cells" = "Myeloid",
  "Tem/Trm cytotoxic T cells" = "T_NK", "Tcm/Naive helper T cells" = "T_NK",
  "Regulatory T cells" = "T_NK", "Type 17 helper T cells" = "T_NK",
  "NK cells" = "T_NK", "CD16- NK cells" = "T_NK", "CD16+ NK cells" = "T_NK", "ILC" = "T_NK",
  "Memory B cells" = "B_Plasma", "Naive B cells" = "B_Plasma", "Plasma cells" = "B_Plasma")
meta2$lineage <- map_lin[meta2$cell_type]
cat("lineage counts:\n"); print(table(meta2$lineage))
stopifnot(!any(is.na(meta2$lineage)))

## 3. 伪bulk：每谱系 colSums（基因 × 谱系），CPM
groups <- factor(meta2$lineage, levels = c("Epithelial", "Stromal", "Endothelial", "Myeloid", "T_NK", "B_Plasma"))
pb <- sapply(levels(groups), function(ct) {
  Matrix::rowSums(counts[, meta2$lineage == ct, drop = FALSE])
})
rownames(pb) <- rownames(counts)
pb <- pb[, levels(groups)]
ref_cpm <- sweep(pb, 2, as.numeric(Matrix::colSums(pb)), "/") * 1e6
saveRDS(list(counts = pb, cpm = ref_cpm, lineage = levels(groups), meta = meta2),
        file.path(oc, "results", "track2", "pb_ovarian_celltype_counts.RDS"))
cat("pseudobulk saved:", dim(pb), "\n")
print(round(head(pb[, 1:3], 5)))
cat("DONE 26a reference\n")
