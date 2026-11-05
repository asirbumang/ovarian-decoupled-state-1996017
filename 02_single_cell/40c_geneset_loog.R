
suppressMessages({library(Matrix); library(data.table); library(Seurat); library(openxlsx)})
set.seed(2026)
mods <- list()
for (line in readLines("D:/子宫内膜异位症2/pipeline/mitoxy_four_modules.gmt")) {
  p <- strsplit(line, "\t")[[1]]
  mods[[p[1]]] <- setdiff(p[3:length(p)], "")
}
mods <- mods[grep("MITOXY_", names(mods))]
cat("modules:", sapply(mods, length), "\n")
## ICC45（JTM 2026 的 45 基因——从 sets rds 拿）
sets <- readRDS("D:/子宫内膜异位症2/pipeline/mitoxy_sets.RDS")
icc <- sets$icc45
cat("ICC45 n:", length(icc), "\n")

## 载入 GSE266577 上皮 20,000 采样（b9_input triplet）
d <- "D:/子宫内膜异位症2/analysis/L1_scRNA/results/cnv/b9_input"
genes <- readLines(file.path(d, "GSE266577_epi.genes.tsv"))
dt <- fread(file.path(d, "GSE266577_epi.triplet.gz"))
gidx <- data.table::chmatch  # noop
setkeyv(as.data.table(data.frame(g = genes, gi = seq_along(genes))), "g")
gene_key <- data.table(g = genes, gi = seq_along(genes)); setkey(gene_key, "g")
dt[, gi := match(g, genes)]
dt[, ci := match(c, unique(c))]
set.seed(2026)
cells <- dt[, unique(ci)]; sel <- sample(cells, min(20000, length(cells)))
sub <- dt[ci %in% sel]
sub[, ci := match(ci, unique(ci))]
C <- Matrix::sparseMatrix(i = sub$gi, j = sub$ci, x = sub$v,
                          dims = c(length(genes), length(unique(sub$ci))))
rownames(C) <- genes
cat("matrix:", dim(C), "\n")
seu <- CreateSeuratObject(counts = C, min.cells = 0, min.features = 0)
seu <- NormalizeData(seu, verbose = FALSE)

score_one <- function(gset) {
  g <- intersect(gset, rownames(seu))
  if (length(g) < 3) return(rep(NA, ncol(seu)))
  seu <<- AddModuleScore(seu, features = list(g), name = "s", nbin = 24, ctrl = 100)
  seu@meta.data$s1
}
full <- lapply(mods, score_one)
names(full) <- names(mods)
fM1 <- full[["MITOXY_M1_immune_activation"]]; fM2 <- full[["MITOXY_M2_metabolic_GSH"]]; fM3 <- full[["MITOXY_M3_mito_cytoskeleton"]]
z <- function(x) as.numeric(scale(x))
full_dec <- z(fM1) - (z(fM2) + z(fM3)) / 2
icc_s <- score_one(icc)
icc_dec <- z(icc_s) - (z(fM2) + z(fM3)) / 2  # ICC45 只有免疫层→ 用同 M2/M3
cat("dual-set: cor(自建M1, ICC45) =", round(cor(fM1, icc_s), 3), "\n")
cat("dual-set: cor(自建decouple, ICC45-decouple) =", round(cor(full_dec, icc_dec), 3), "\n")
cat("Venn: 自建四模块 ∩ ICC45 =", length(intersect(unique(unlist(mods)), icc)), "\n")

## LOOG：逐基因剔除（每模块），报告 decouple 与原版 Spearman
l1 <- unlist(mods[["MITOXY_M1_immune_activation"]]); l2 <- unlist(mods[["MITOXY_M2_metabolic_GSH"]]); l3 <- unlist(mods[["MITOXY_M3_mito_cytoskeleton"]])
res <- data.frame()
for (gene in unique(c(l1, l2, l3))) {
  m1 <- full[[1]]; m2 <- full[[2]]; m3 <- full[[3]]
  if (gene %in% l1) { m1 <- score_one(setdiff(l1, gene)) } else if (gene %in% l2) { m2 <- score_one(setdiff(l2, gene)) } else { m3 <- score_one(setdiff(l3, gene)) }
  d <- z(m1) - (z(m2) + z(m3)) / 2
  r <- cor(d, full_dec, method = "spearman")
  res <- rbind(res, data.frame(gene = gene, module = ifelse(gene %in% l1, "M1", ifelse(gene %in% l2, "M2", "M3")), sp = round(r, 4)))
}
cat("\nLOOG min/median/max sp:", round(min(res$sp), 4), round(median(res$sp), 4), round(max(res$sp), 4), "\n")
print(res[order(res$sp), ][1:8, ])
fwrite(res, "D:/子宫内膜异位症2/analysis/L1_scRNA/results/cnv/geneset_loog.csv")
cat("DONE 40c\n")
