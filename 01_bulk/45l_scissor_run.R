
suppressMessages({library(Scissor); library(data.table)})
basedir <- "D:/子宫内膜异位症2/analysis/L0_bulk/results/track2/mime"
sc <- readRDS(file.path(basedir, "scissor_scRNA.rds"))
e <- as.matrix(sc$expr); rownames(e) <- sc$genes; colnames(e) <- sc$cells
## Scissor: 二者均需 基因×样本（基因为行）；抽样 60k 控内存（保 state 比例：用 state 标签缺失时随机）
set.seed(42)
sel_cells <- sample(colnames(e), 60000)
e <- e[, sel_cells]
cat("scRNA (genes x cells):", dim(e), "\n")
bt <- fread(file.path(basedir, "tcga_exp.csv.gz"), data.table = FALSE)
rownames(bt) <- bt[[1]]; bt <- bt[, -1]
cl <- fread(file.path(basedir, "tcga_clin.csv"))
samples <- intersect(cl$sample, colnames(bt))
bt <- as.matrix(bt[, samples])
cl <- cl[cl$sample %in% samples, ]
surv <- data.frame(OS.time = cl$futime, OS = cl$fustat)
rownames(surv) <- cl$sample
cat("bulk:", dim(bt), " surv n:", nrow(surv), " events:", sum(surv$OS), "\n")
common <- intersect(rownames(e), rownames(bt))
cat("common genes:", length(common), "\n")
e <- e[common, ]; bt <- bt[common, ]
e <- matrix(as.numeric(e), nrow = nrow(e), dimnames = dimnames(e))
bt <- matrix(as.numeric(bt), nrow = nrow(bt), dimnames = dimnames(bt))
## 最小兼容补丁：V5（class 双值 if + data layer）
txt <- paste(deparse(Scissor), collapse = "\n")
txt <- sub('Seurat_tmp <- ScaleData\\(Seurat_tmp, verbose = F\\)',
           'Seurat_tmp <- Seurat::NormalizeData(Seurat_tmp, verbose = F); Seurat_tmp <- ScaleData(Seurat_tmp, verbose = F)', txt)
txt <- sub('class\\(sc_dataset\\) == "Seurat"', 'inherits(sc_dataset, "Seurat")', txt)
txt <- sub('class\\(bulk_dataset\\) == "Seurat"', 'inherits(bulk_dataset, "Seurat")', txt)
eval(parse(text = paste0("patched <- ", txt)))
cat("running Scissor (alpha=0.05, cox; V5-compat)...\n")
res <- patched(bulk_dataset = bt, sc_dataset = e, phenotype = surv,
               family = "cox", alpha = 0.05, Save_file = file.path(basedir, "scissor_inputs.RData"))
cat("pos:", length(res$Scissor_pos), " neg:", length(res$Scissor_neg), "\n")
saveRDS(res, file.path(basedir, "scissor_res_alpha05.rds"))
if (length(res$Scissor_pos) > 5) {
  res2 <- patched(bulk_dataset = bt, sc_dataset = e, phenotype = surv,
                  family = "cox", alpha = 0.01, Save_file = file.path(basedir, "scissor_inputs01.RData"))
  cat("alpha01 pos:", length(res2$Scissor_pos), " neg:", length(res2$Scissor_neg), "\n")
  saveRDS(res2, file.path(basedir, "scissor_res_alpha01.rds"))
}
cat("DONE 45l\n")
