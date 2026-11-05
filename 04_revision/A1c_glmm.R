# A1c: patient random-intercept logistic models (lme4) for carrier enrichment and ICC
suppressMessages({library(lme4); library(data.table)})
d <- fread("D:/子宫内膜异位症2/analysis/L1_scRNA/results/module_scores/scores_GSE266577.csv",
           data.table = FALSE)
d$patient <- as.character(d$publication_patient_code_final)
d$state <- as.integer(d$state_high)
d$plasma <- as.integer(d$cell_type == "Plasma cells")
d$fibro <- as.integer(d$cell_type == "Fibroblasts")

out <- list()
for (v in c("plasma", "fibro")) {
  f <- as.formula(paste0("state ~ ", v, " + (1|patient)"))
  m <- glmer(f, data = d, family = binomial, nAGQ = 0,
             control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
  b <- fixef(m)[[v]]; se <- sqrt(vcov(m)[v, v])
  vp <- as.numeric(VarCorr(m)$patient[1, 1])
  icc <- vp / (vp + pi^2 / 3)
  out[[v]] <- list(OR = exp(b), CI_low = exp(b - 1.96 * se), CI_high = exp(b + 1.96 * se),
                   p = 2 * pnorm(-abs(b / se)), var_patient = vp, ICC = icc,
                   n = nrow(d), n_patient = length(unique(d$patient)))
}
cat(jsonlite::toJSON(out, auto_unbox = TRUE, pretty = TRUE))
writeLines(jsonlite::toJSON(out, auto_unbox = TRUE, pretty = TRUE),
           "D:/子宫内膜异位症2/R2_revision/04_analysis/results/A1c_glmm.json")
