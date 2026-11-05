# ============================================================
# 轨道0 补充统计（06 核验文档引用）：分层中位OS/反应/残留病灶/年龄/ICC45 双轨
# 输入：results/track0_tcga_scores_states.csv
# ============================================================
suppressMessages({library(data.table); library(survival); library(ggplot2)})
d <- fread("D:/子宫内膜异位症2/analysis/L0_bulk/results/track0_tcga_scores_states.csv", data.table = FALSE)
cat("n =", nrow(d), "\n")

## 1) 三分位分层：n / 事件 / 中位 OS / 95%CI
d$grp <- factor(d$decouple_tier, levels = c("Low", "Mid", "High"))
sd <- d[!is.na(d$os_days) & d$os_days > 0, ]
cat("\n== tier (OS cohort n=", nrow(sd), ") ==\n")
for (g in c("Low", "Mid", "High")) {
  sub <- sd[sd$grp == g, ]
  q <- survfit(Surv(os_days, event) ~ 1, data = sub)
  m <- summary(q)$table[c("median", "0.95LCL", "0.95UCL")]
  cat(sprintf("%-5s n=%3d events=%3d  median OS=%5.0f [%5.0f-%5.0f] days\n",
              g, nrow(sub), sum(sub$event), m[1], m[2], m[3]))
}
lr <- survdiff(Surv(os_days, event) ~ grp, data = sd)
cat("logrank p =", format(pchisq(lr$chisq, length(lr$n) - 1, lower.tail = FALSE), digits = 3), "\n")
## 趋势检验：把 tier 编码为 1/2/3 用 score 检验（Cox 单变量）
sd$tier_num <- as.numeric(sd$grp)
ct <- coxph(Surv(os_days, event) ~ tier_num, data = sd)
cat("trend cox HR(per tier) =", round(exp(coef(ct)), 3), " p =", round(summary(ct)$coef[5], 4), "\n")

## 2) 化疗反应 4 组：n / decouple 中位 / Kruskal p
cat("\n== response groups ==\n")
print(table(d$outcome_bin, useNA = "ifany"))
krus <- kruskal.test(decouple ~ outcome_bin, data = d)
cat("Kruskal p =", format(krus$p.value, digits = 3), "\n")
for (g in c("CR", "PR", "SD", "PD")) {
  s <- d$decouple[d$outcome_bin == g]
  cat(sprintf("%-3s n=%3d  median=%.3f  IQR=(%.3f,%.3f)\n", g, length(s), median(s), quantile(s, .25), quantile(s, .75)))
}

## 3) 残留病灶（ov:tumor_residual_disease 分级：无宏观 → 1-10 → 11-20 → >20 mm）
cat("\n== residual tumor ==\n")
sub <- d[!is.na(d$residual_tumor) & d$residual_tumor != "", ]
ord <- c("No Macroscopic disease", "1-10 mm", "11-20 mm", ">20 mm")
sub$rt <- factor(sub$residual_tumor, levels = ord)
sub <- sub[!is.na(sub$rt), ]  # 剔除历史遗留 R0/RX 等非规范值
cat("counts:\n"); print(table(sub$rt))
cat("Kruskal p =", format(kruskal.test(decouple ~ rt, data = sub)$p.value, digits = 3), "\n")
for (g in ord) {
  s <- sub$decouple[sub$rt == g]
  cat(sprintf("%-24s n=%3d median=%.3f\n", g, length(s), median(s)))
}
sub$rt_ord <- as.numeric(sub$rt)
cat("Spearman r(ordinal residual, decouple) =", round(suppressWarnings(cor(sub$decouple, sub$rt_ord, method = "spearman", use = "complete.obs")), 3),
    " p =", format(suppressWarnings(cor.test(sub$decouple, sub$rt_ord, method = "spearman")$p.value), digits = 3), "\n")
sub$rt2 <- ifelse(sub$rt_ord == 1, "No_Macro", "Residual")
w2 <- wilcox.test(decouple ~ rt2, data = sub)
cat("Wilcoxon 无宏观 vs 有残留 p =", format(w2$p.value, digits = 3), "\n")

## 4) 年龄
cat("\n== age ==\n")
cat("Spearman r(decouple, age) =", round(suppressWarnings(cor(d$decouple, d$age_at_dx, method = "spearman", use = "complete.obs")), 3),
    " p =", format(suppressWarnings(cor.test(d$decouple, d$age_at_dx, method = "spearman")$p.value), digits = 3), "\n")

## 5) 组织学
cat("\n== histology ==\n")
print(table(d$histological_type, useNA = "ifany"))

## 6) ICC45 双轨对照：解耦分数 vs ICC45（免疫检查点/细胞毒）—— 不是"免疫热即可"论证
cat("\n== ICC45 vs modules ==\n")
ic <- readRDS("D:/子宫内膜异位症2/analysis/L0_bulk/results/track0_tcga_stats.RDS")
for (m in c("M1_immune_activation", "M2_metabolic_GSH", "M3_mito_cytoskeleton", "M4_lipid_perox_rupture")) {
  ct2 <- suppressWarnings(cor.test(d$ICC45, d[[m]], method = "spearman"))
  cat(sprintf("%-24s r=%.3f p=%.2e\n", m, ct2$estimate, ct2$p.value))
}

## 7) 多因素 Cox：decouple + 年龄 + 残留病灶（R0 vs R1+）
cat("\n== multivariable cox ==\n")
sd$age <- sd$age_at_dx
sd$rt2 <- ifelse(sd$residual_tumor == "No Macroscopic disease", 0,
                 ifelse(sd$residual_tumor %in% c("1-10 mm", "11-20 mm", ">20 mm"), 1, NA))
mx <- tryCatch(coxph(Surv(os_days, event) ~ decouple + age + rt2, data = sd), error = function(e) NULL)
if (!is.null(mx)) {
  sm <- summary(mx)
  print(round(sm$coefficients[, c("coef", "se(coef)", "Pr(>|z|)")], 4))
  cat("n used =", mx$n, "\n")
} else cat("cox failed\n")
cat("\nDONE supp stats\n")
