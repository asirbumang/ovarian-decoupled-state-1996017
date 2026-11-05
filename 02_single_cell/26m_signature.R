
suppressMessages({library(Matrix); library(data.table)})
set.seed(2026)
oc <- "D:/子宫内膜异位症2/analysis/L1_scRNA"
sdir <- file.path(oc, "results", "module_scores")

## 1. 上皮层内分数（含 decouple, 26b 产出）
epi <- fread(file.path(sdir, "scores_GSE266577_epithelial.csv"), data.table = FALSE)
cat("epi cells:", nrow(epi), "\n")

## 2. 计数（只需上皮细胞）
counts <- readMM(gzfile("D:/子宫内膜异位症2/data_sc/GSE266577/GSE266577_counts_raw.mtx.gz"))
counts <- as(counts, "dgCMatrix")
feat <- fread("D:/子宫内膜异位症2/data_sc/GSE266577/GSE266577_seurat_features.txt.gz", header = FALSE, data.table = FALSE)
bar <- fread("D:/子宫内膜异位症2/data_sc/GSE266577/GSE266577_barcodes.txt.gz", header = FALSE, data.table = FALSE)
rownames(counts) <- feat$V1; colnames(counts) <- bar$V1
sub <- counts[, epi$cell]
cat("subset:", dim(sub), "\n")
rm(counts); gc()

## 3. 组：top/bottom quartile（decouple，层内已有）
q25 <- quantile(epi$decouple, .25); q75 <- quantile(epi$decouple, .75)
hi <- epi$cell[epi$decouple >= q75]; lo <- epi$cell[epi$decouple <= q25]
cat("hi:", length(hi), " lo:", length(lo), "\n")

cpm_hi <- Matrix::rowMeans(sweep(sub[, hi], 2, Matrix::colSums(sub[, hi]) / 1e6, "/"))
cpm_lo <- Matrix::rowMeans(sweep(sub[, lo], 2, Matrix::colSums(sub[, lo]) / 1e6, "/"))
lhi <- log2(cpm_hi + 1); llo <- log2(cpm_lo + 1)
d <- data.frame(gene = rownames(sub), mean_hi = lhi, mean_lo = llo, stringsAsFactors = FALSE)
d$diff <- d$mean_hi - d$mean_lo
## 过滤低表达（两组都 < 0.5 dim）
d <- d[d$mean_hi > 0.3 | d$mean_lo > 0.3, ]
d <- d[order(-abs(d$diff)), ]
sig_up <- d$gene[d$diff > 0][1:200]
sig_dn <- d$gene[d$diff < 0][1:200]
cat("sig up:", length(sig_up), " dn:", length(sig_dn), "\n")
fwrite(d, file.path(oc, "results", "cnv", "state_signature_full.csv"))
writeLines(sig_up, file.path(oc, "results", "cnv", "state_sig_up.txt"))
writeLines(sig_dn, file.path(oc, "results", "cnv", "state_sig_dn.txt"))
cat("DONE 26m\n")
