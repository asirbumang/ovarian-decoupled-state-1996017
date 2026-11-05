
suppressMessages({library(infercnv); library(data.table)})
rd <- "D:/子宫内膜异位症2/analysis/L1_scRNA/results/cnv/run"
out <- "D:/子宫内膜异位症2/analysis/L1_scRNA/results/cnv"
obj <- readRDS(file.path(rd, "15_tumor_subclusters.leiden.infercnv_obj"))
cells <- colnames(obj@expr.data)
cat("expr:", dim(obj@expr.data), "\n")
vals <- as.matrix(obj@expr.data)
rm(obj); gc()
absm <- abs(vals)
cnv_all <- colMeans(absm)
v <- apply(vals, 1, var)
top_genes <- names(sort(v, decreasing = TRUE))[1:3000]
cnv_top <- colMeans(absm[top_genes, ])
obj2 <- readRDS(file.path(rd, "run.final.infercnv_obj"))
grp <- rep("Ref", length(cells)); names(grp) <- cells
grp[cells[unlist(obj2@observation_grouped_cell_indices)]] <- "Epithelial"
rm(obj2); gc()
res <- data.frame(cell = cells, cnv_all = cnv_all, cnv_top = cnv_top, group = grp)
q95 <- quantile(res$cnv_top[res$group == "Ref"], 0.95)
q99 <- quantile(res$cnv_top[res$group == "Ref"], 0.99)
cat("ref q95/q99 (cnv_top):", round(q95, 3), round(q99, 3), "\n")
res$malignant_like <- res$cnv_top > q99
epi <- res[res$group == "Epithelial", ]
cat("Epithelial:", nrow(epi), " malignant_like:", sum(epi$malignant_like),
    paste0("(", round(100*mean(epi$malignant_like), 1), "%)\n"))
fwrite(res, file.path(out, "percell_cnv_score_final.csv"))
cat("\ncnv_top by group:\n")
print(aggregate(cnv_top ~ group, data = res, FUN = function(x) round(quantile(x, c(.5, .9, .99)), 3)))
cat("DONE 26n final\n")
