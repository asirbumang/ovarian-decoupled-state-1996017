suppressMessages({library(data.table)})
d <- fread("D:/子宫内膜异位症2/analysis/L0_bulk/results/survival_meta/scores_GSE63885.csv", data.table = FALSE)
cat("n:", nrow(d), "\n")
for (v in c("brca1", "chemo")) {
  d2 <- d[!is.na(d[[v]]) & d[[v]] != "", ]
  if (nrow(d2) >= 10) {
    d2$grp <- ifelse(grepl("mutation", d2[[v]]), "mut", "no")
    w <- tryCatch(wilcox.test(decouple ~ grp, data = d2), error = function(e) NULL)
    cat(sprintf("%s: %s: n=%d p=%.3f (mut med %.3f vs no med %.3f)\n", v, v, nrow(d2),
                if (!is.null(w)) w$p.value else NA,
                median(d2$decouple[d2$grp == "mut"]), median(d2$decouple[d2$grp == "no"])))
  }
}
## 残灶 R0 vs R1/R2
d$rt <- d$`residual tumor size (r0 < 1 cm;  1 < r1 < 5 cm; r2 > 5cm)`
if (is.null(d$rt)) d$rt <- NA
cat("residual col missing, skip\n")
