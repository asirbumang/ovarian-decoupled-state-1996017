
suppressMessages({library(Matrix); library(Seurat); library(data.table); library(ggplot2); library(ggsci); library(ragg)})
set.seed(42)
sdir <- "D:/子宫内膜异位症2/analysis/L1_scRNA/results/module_scores"
d <- "D:/子宫内膜异位症2/data_sc/GSE293459/merged10x"
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

cat("loading counts...\n")
mtx <- readMM(file.path(d, "matrix.mtx")); mtx <- as(mtx, "dgCMatrix")
fea <- fread(file.path(d, "features.tsv"), header = FALSE, data.table = FALSE)
bar <- fread(file.path(d, "barcodes.tsv"), header = FALSE, data.table = FALSE)
rownames(mtx) <- fea$V1; colnames(mtx) <- bar$V1
cat("matrix:", dim(mtx), "\n")
## QC（snRNA：>200 基因 & <5% mt）
if (anyDuplicated(rownames(mtx))) mtx <- mtx[!duplicated(rownames(mtx)), ]
nf <- Matrix::colSums(mtx > 0); nc <- Matrix::colSums(mtx)
mt_genes <- grep("^MT-", rownames(mtx), value = TRUE); 
mt_frac <- Matrix::colSums(mtx[mt_genes, , drop = FALSE]) / (nc + 1e-9)
keep <- nf >= 200 & mt_frac < 0.05
mtx <- mtx[, keep]
cat("after QC:", ncol(mtx), "\n")
seu <- CreateSeuratObject(counts = mtx, min.cells = 0, min.features = 0)
seu <- NormalizeData(seu, verbose = FALSE)
si <- lapply(score_sets, function(g) intersect(g, rownames(seu)))
cat("coverage:", paste(sprintf("%s=%d", names(si), sapply(si, length)), collapse = " | "), "\n")
seu <- AddModuleScore(seu, features = si, name = "sc", nbin = 24, ctrl = 100)
cols <- grep("^sc[0-9]+$", colnames(seu@meta.data), value = TRUE)
df <- data.frame(cell = colnames(seu), stringsAsFactors = FALSE)
for (k in seq_along(score_sets)) df[[names(score_sets)[k]]] <- seu@meta.data[[cols[k]]]
meta <- fread("D:/子宫内膜异位症2/data_sc/GSE293459/meta_cells.csv", data.table = FALSE)
df <- merge(df, meta[, c("cell", "sample", "sample_type")], by = "cell", all.x = TRUE)
df <- as.data.table(df)
df[, M1z := as.numeric(scale(M1_immune_activation))]
df[, M2z := as.numeric(scale(M2_metabolic_GSH))]
df[, M3z := as.numeric(scale(M3_mito_cytoskeleton))]
df[, decouple := M1z - (M2z + M3z) / 2]
fwrite(df, file.path(sdir, "scores_GSE293459.csv"))
cat("cells:", nrow(df), "\n")
## C vs T（样本级 + 细胞级；R=1303_R 单列）
samp <- df[, .(n_cells = .N, med = median(decouple)), by = .(sample, sample_type)]
print(as.data.frame(samp))
w <- wilcox.test(med ~ sample_type, data = as.data.frame(samp[sample_type %in% c("C", "T")]))
cat("SAMPLE-level C vs T: p=", format(w$p.value, digits = 4), " med:", median(samp$med[samp$sample_type == "C"]),
    " vs ", median(samp$med[samp$sample_type == "T"]), "\n", sep = "")
wc <- wilcox.test(decouple ~ sample_type, data = as.data.frame(df[sample_type %in% c("C", "T")]))
cat("CELL-level C vs T: p=", format(wc$p.value, digits = 4), "\n")
cat("DONE 27b\n")
