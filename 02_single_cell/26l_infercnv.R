# ============================================================
# 26l: GSE266577 恶性肿瘤细胞确认（inferCNV）
# 输入：稀疏计数（137k）；参考 = 免疫/基质（按比例抽样 20k）
# 输出：inferCNV 运行目录 + 上皮细胞 per-cell CNV 分数 → malignant 标记
# ============================================================
suppressMessages({
  library(Matrix); library(data.table); library(infercnv)
})
set.seed(2026)
oc <- "D:/子宫内膜异位症2/analysis/L1_scRNA"
cndir <- file.path(oc, "results", "cnv")
dir.create(cndir, recursive = TRUE, showWarnings = FALSE)

## 1. 计数 + meta
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

## 2. 选择：上皮全部 + 参考（免疫/基质）20k
epi <- which(meta2$cell_type == "Epithelial cells")
ref_types <- c("Tem/Trm cytotoxic T cells", "Tcm/Naive helper T cells", "Regulatory T cells",
               "Macrophages", "Fibroblasts", "B cells")
ref_cells <- which(meta2$cell_type %in% c("Tem/Trm cytotoxic T cells", "Tcm/Naive helper T cells",
                                          "Regulatory T cells", "Macrophages", "Fibroblasts",
                                          "Memory B cells", "Naive B cells"))
if (length(ref_cells) > 20000) ref_cells <- sample(ref_cells, 20000)
sel <- c(epi, ref_cells)
cat("epithelial:", length(epi), " ref:", length(ref_cells), "\n")
sub <- counts[, sel]
mt <- meta2[sel, ]
ann <- data.frame(cell = colnames(sub), group = ifelse(mt$cell_type == "Epithelial cells", "Epithelial", "Ref"))
write.table(ann, file.path(cndir, "anns.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

## 3. gene order（if 文件缺失则报错）
go <- fread("D:/子宫内膜异位症2/data_omics/gene_order_hg38.tsv", data.table = FALSE)
genes_in <- intersect(rownames(sub), go$gene)
sub <- sub[genes_in, ]
go <- go[match(genes_in, go$gene), ]
order_df <- data.frame(gene = genes_in, chr = go$chr, start = go$start, end = go$end)
write.table(order_df, file.path(cndir, "genes.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
cat("genes for infercnv:", nrow(order_df), "\n")
saveRDS(sub, file.path(cndir, "sub_counts.RDS"))

## 4. 运行（无 HMM 快版本；denoise + 标准参数）
obj <- CreateInfercnvObject(raw_counts_matrix = sub,
                            annotations_file = file.path(cndir, "anns.txt"),
                            gene_order_file = file.path(cndir, "genes.txt"),
                            ref_group_names = c("Ref"))
obj <- infercnv::run(obj, cutoff = 0.1, out_dir = file.path(cndir, "run"),
                     cluster_by_groups = TRUE, denoise = TRUE, scale_data = TRUE,
                     analysis_mode = "subclusters", HMM = FALSE, num_threads = 4)
cat("DONE 26l inferCNV run\n")
