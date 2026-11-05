# ============================================================
# 26f: GSE266577 逐类型 marker（mean logCPM + pct_in）+ 规范标记面板导出
# 供：DeepSeek 注释三审（第一遍）+ 后续报告
# ============================================================
suppressMessages({library(Matrix); library(data.table)})
oc <- "D:/子宫内膜异位症2/analysis/L1_scRNA"
outdir <- file.path(oc, "results", "annotation_review")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

## 1. 载入
counts <- readMM(gzfile("D:/子宫内膜异位症2/data_sc/GSE266577/GSE266577_counts_raw.mtx.gz"))
counts <- as(counts, "dgCMatrix")
feat <- fread("D:/子宫内膜异位症2/data_sc/GSE266577/GSE266577_seurat_features.txt.gz", header = FALSE, data.table = FALSE)
bar <- fread("D:/子宫内膜异位症2/data_sc/GSE266577/GSE266577_barcodes.txt.gz", header = FALSE, data.table = FALSE)
rownames(counts) <- feat$V1; colnames(counts) <- bar$V1
meta <- fread("D:/子宫内膜异位症2/data_sc/GSE266577/GSE266577_metadata.txt.gz", data.table = FALSE)
rownames(meta) <- meta$cell_name
keep <- meta[colnames(counts), ]
qc <- as.numeric(keep$nFeature_RNA) >= 200 & as.numeric(keep$percent.mt) <= 20
counts <- counts[, qc]
ct <- keep$cell_type[qc]
cat("cells:", ncol(counts), "\n")

## 2. 全基因 × 类型 均值（logCPM）与 pct_in
cts <- sort(unique(ct))
mean_by <- matrix(NA, nrow(counts), length(cts), dimnames = list(rownames(counts), cts))
pct_by <- mean_by
for (k in seq_along(cts)) {
  sub <- counts[, ct == cts[k], drop = FALSE]
  lib <- Matrix::colSums(sub)
  cpm <- sweep(sub, 2, lib / 1e6, "/")
  mean_by[, k] <- Matrix::rowMeans(cpm)
  pct_by[, k] <- Matrix::rowMeans(sub > 0)
  cat("  ", cts[k], "cells:", ncol(sub), "\n", sep = "")
}
saveRDS(list(mean_cpm = mean_by, pct = pct_by, ct = cts), file.path(outdir, "GSE266577_type_means.RDS"))
cat("saved means\n")

## 3. top 差异 marker（mean 排名 vs 其他类型最大均值）
spec <- log(mean_by + 1) - apply(log(mean_by + 1), 1, function(x) sort(x, decreasing = TRUE)[2])
top_markers <- lapply(colnames(spec), function(ct) {
  s <- spec[, ct]
  g <- names(head(sort(s, decreasing = TRUE), 15))
  rbindlist(lapply(g, function(gene) data.table(ct = ct, gene = gene,
    mean_cpm = mean_by[gene, ct], pct_in = pct_by[gene, ct],
    sdiff = s[gene])))
})
top_markers <- rbindlist(top_markers)
fwrite(top_markers, file.path(outdir, "GSE266577_top_markers_per_type.csv"))

## 4. 规范标记面板表（审核用）
canon <- list(
  "Epithelial cells" = c("EPCAM","KRT8","KRT18","KRT19","PAX8","MUC1","WFDC2","MUC16","KRT7","MSLN"),
  "Fibroblasts" = c("COL1A1","COL1A2","DCN","PDGFRA","LUM","GREM1","FAP"),
  "Mesothelial cells" = c("MSLN","CALB2","WT1","UPK3B","KRT8"),
  "Endothelial cells" = c("PECAM1","VWF","CDH5","CLDN5","ENG","PLVAP"),
  "Macrophages" = c("CD68","CD163","CSF1R","ITGAM","C1QA","C1QB","MRC1","TYROBP","FCGR3A"),
  "Classical monocytes" = c("CD14","FCN1","S100A8","S100A9","LYZ","VCAN"),
  "Migratory DCs" = c("CCR7", "LAMP3", "CCL19", "CCL21"),
  "DC1" = c("XCR1","CLEC9A","C1QC","BATF3"),
  "DC2" = c("CD1C","FCER1A","CLEC10A"),
  "pDC" = c("LILRA4","GZMB","CLEC4C","TCF4"),
  "Mast cells" = c("TPSAB1","TPSB2","CPA3","KIT","MS4A2"),
  "Tem/Trm cytotoxic T cells" = c("CD8A","CD8B","GZMA","GZMB","GZMK","PRF1","NKG7","IFNG","CD3D"),
  "Tcm/Naive helper T cells" = c("CD4","IL7R","SELL","CCR7","LEF1","CD3E"),
  "Regulatory T cells" = c("FOXP3","CTLA4","IL2RA","IKZF2","TNFRSF18"),
  "Type 17 helper T cells" = c("RORC","IL17A","IL23R"),
  "NK cells" = c("NKG7","KLRD1","GNLY","NCAM1","KLRF1","GZMB","PRF1"),
  "CD16- NK cells" = c("SELL","TGFB1","CCL4","KLRD1"),
  "CD16+ NK cells" = c("FCGR3A","PRF1","GZMB","CX3CR1"),
  "ILC" = c("IL7R","KIT","GATA3","PTGDR2","TRDC"),
  "Memory B cells" = c("CD19","MS4A1","CD79A","CD79B","CD27","CD38"),
  "Naive B cells" = c("CD19","MS4A1","CD79A","CD79B","IGHM"),
  "Plasma cells" = c("MZB1","SDC1","CD38","XBP1","IGHG1","JCHAIN"))
rows <- list()
for (ct in names(canon)) {
  for (g in intersect(canon[[ct]], rownames(counts))) {
    rows[[length(rows) + 1]] <- data.frame(cell_type = ct, gene = g,
      mean_cpm = mean_by[g, ct], pct_in = pct_by[g, ct])
  }
}
panel <- rbindlist(rows)
fwrite(panel, file.path(outdir, "GSE266577_canonical_panel.csv"))
cat("canonical panel:", nrow(panel), "entries\n")
cat("DONE 26f\n")
