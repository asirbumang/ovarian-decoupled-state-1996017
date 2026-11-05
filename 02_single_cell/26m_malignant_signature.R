
suppressMessages({library(Matrix); library(data.table)})
set.seed(2026)
oc <- "D:/子宫内膜异位症2/analysis/L1_scRNA"
sdir <- file.path(oc, "results", "module_scores")
cndir <- file.path(oc, "results", "cnv")
epi <- fread(file.path(sdir, "scores_GSE266577_epithelial.csv"), data.table = FALSE)
cnv <- fread(file.path(cndir, "percell_cnv_score_v2.csv"), data.table = FALSE)
cnv <- cnv[cnv$malignant_like, ]
mal <- merge(epi, cnv[, c("cell", "cnv_top")], by = "cell")
cat("malignant cells:", nrow(mal), "\n")
## 计数子集
counts <- readMM(gzfile("D:/子宫内膜异位症2/data_sc/GSE266577/GSE266577_counts_raw.mtx.gz"))
counts <- as(counts, "dgCMatrix")
feat <- fread("D:/子宫内膜异位症2/data_sc/GSE266577/GSE266577_seurat_features.txt.gz", header = FALSE, data.table = FALSE)
bar <- fread("D:/子宫内膜异位症2/data_sc/GSE266577/GSE266577_barcodes.txt.gz", header = FALSE, data.table = FALSE)
rownames(counts) <- feat$V1; colnames(counts) <- bar$V1
sub <- counts[, mal$cell]
rm(counts); gc()
q25 <- quantile(mal$decouple, .25); q75 <- quantile(mal$decouple, .75)
hi <- mal$cell[mal$decouple >= q75]; lo <- mal$cell[mal$decouple <= q25]
cat("hi:", length(hi), " lo:", length(lo), "\n")
cpm_hi <- Matrix::rowMeans(sweep(sub[, hi], 2, Matrix::colSums(sub[, hi]) / 1e6, "/"))
cpm_lo <- Matrix::rowMeans(sweep(sub[, lo], 2, Matrix::colSums(sub[, lo]) / 1e6, "/"))
lhi <- log2(cpm_hi + 1); llo <- log2(cpm_lo + 1)
d <- data.frame(gene = rownames(sub), mean_hi = lhi, mean_lo = llo, stringsAsFactors = FALSE)
d$diff <- d$mean_hi - d$mean_lo
d <- d[d$mean_hi > 0.3 | d$mean_lo > 0.3, ]
d <- d[order(-abs(d$diff)), ]
sig_up <- d$gene[d$diff > 0][1:200]
sig_dn <- d$gene[d$diff < 0][1:200]
writeLines(sig_up, file.path(cndir, "state_sig_up_mal.txt"))
writeLines(sig_dn, file.path(cndir, "state_sig_dn_mal.txt"))
cat("mal sig up:", length(sig_up), " dn:", length(sig_dn), "\n")
cat("overlap with original up:", length(intersect(sig_up, readLines(file.path(cndir, "state_sig_up.txt")))), "\n")
cat("DONE 26m-mal\n")
