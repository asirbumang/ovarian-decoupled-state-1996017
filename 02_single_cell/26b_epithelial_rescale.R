# ============================================================
# 26b: GSE266577 恶性（上皮）细胞专属再标定 —— 状态/解耦在"恶性层内"重新定义
# 问题：全细胞 z 标准化以免疫细胞为主，上皮细胞 M1 整体偏低 → 层内定义更准确
# ============================================================
suppressMessages({library(data.table)})
sdir <- "D:/子宫内膜异位症2/analysis/L1_scRNA/results/module_scores"
sc <- fread(file.path(sdir, "scores_GSE266577.csv"), data.table = FALSE)
epi <- sc[sc$cell_type == "Epithelial cells", ]
cat("epithelial cells:", nrow(epi), "\n")

## 层内 z 标准化（M1/M2/M3 在上皮细胞间标准化）
epi$M1z <- as.numeric(scale(epi$M1_immune_activation))
epi$M2z <- as.numeric(scale(epi$M2_metabolic_GSH))
epi$M3z <- as.numeric(scale(epi$M3_mito_cytoskeleton))
epi$decouple <- epi$M1z - (epi$M2z + epi$M3z) / 2
## 层内状态（P75/P25/P25 与 P60/P40/P40）
epi$state_strict <- epi$M1_immune_activation >= quantile(epi$M1_immune_activation, .75) &
                    epi$M2_metabolic_GSH <= quantile(epi$M2_metabolic_GSH, .25) &
                    epi$M3_mito_cytoskeleton <= quantile(epi$M3_mito_cytoskeleton, .25)
epi$state_relax <- epi$M1_immune_activation >= quantile(epi$M1_immune_activation, .60) &
                   epi$M2_metabolic_GSH <= quantile(epi$M2_metabolic_GSH, .40) &
                   epi$M3_mito_cytoskeleton <= quantile(epi$M3_mito_cytoskeleton, .40)
fwrite(epi, file.path(sdir, "scores_GSE266577_epithelial.csv"))

cat("\n== 层内状态频率 ==\n")
cat("strict:", sum(epi$state_strict), "=", round(100 * mean(epi$state_strict), 2), "%\n")
cat("relax :", sum(epi$state_relax), "=", round(100 * mean(epi$state_relax), 2), "%\n")

cat("\n== 层内解耦 × 化疗阶段 ==\n")
for (st in c("chemo-naive", "IDS")) {
  s <- epi$decouple[epi$treatment_stage == st]
  cat(sprintf("%-12s n=%5d median=%.3f  IQR=(%.3f, %.3f)  strict=%.2f%%\n",
              st, length(s), median(s), quantile(s, .25), quantile(s, .75),
              100 * mean(epi$state_strict[epi$treatment_stage == st])))
}
w <- wilcox.test(decouple ~ treatment_stage, data = epi)
cat("Wilcoxon (层内 z) p =", format(w$p.value, digits = 4), "\n")
## 患者级（同一患者化疗前/后均有上皮细胞者）
epi <- as.data.table(epi)
pairs <- epi[, .(med = median(decouple)), by = .(patient = publication_patient_code_final, treatment_stage)]
w2 <- wilcox.test(med ~ treatment_stage, data = as.data.frame(pairs))
cat("\n== 患者级（中位解耦，患者×阶段）==\n")
for (st in c("chemo-naive", "IDS")) {
  cat(sprintf("%-12s 患者数=%d 中位=%.3f\n", st, sum(pairs$treatment_stage == st), median(pairs$med[pairs$treatment_stage == st])))
}
cat("患者级 Wilcoxon p =", format(w2$p.value, digits = 4), "\n")
## 配对（同患者前后均有的子集；患者编码 Sxxx）
pw <- dcast(pairs, patient ~ treatment_stage, value.var = "med")
pw <- pw[!is.na(`chemo-naive`) & !is.na(IDS), ]
cat("paired patients:", nrow(pw), "\n")
if (nrow(pw) >= 3) {
  wp <- wilcox.test(pw$IDS, pw$`chemo-naive`, paired = TRUE)
  cat("paired Wilcoxon 化疗前 vs IDS p =", format(wp$p.value, digits = 4), "\n")
}
fwrite(as.data.frame(epi), file.path(sdir, "scores_GSE266577_epithelial.csv"))
cat("DONE 26b\n")
