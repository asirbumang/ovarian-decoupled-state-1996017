
suppressMessages({library(data.table); library(AnnotationDbi); library(hgu133plus2.db); library(openxlsx)})
lines <- readLines(gzfile("D:/子宫内膜异位症2/data_bulk/GSE26193_series_matrix.txt.gz", "rt"), warn = FALSE)
beg <- which(startsWith(lines, "!series_matrix_table_begin"))[1]
end <- which(startsWith(lines, "!series_matrix_table_end"))[1]
hdr <- strsplit(sub('"', "", lines[beg + 1]), "\t")[[1]]
dat <- read.table(text = lines[(beg + 2):(end - 1)], sep = "\t", quote = "\"", header = FALSE,
                  row.names = 1, check.names = FALSE)
colnames(dat) <- gsub('"', "", hdr[-1])
char <- list()
for (i in seq_len(beg)) {
  if (startsWith(lines[i], "!Sample_characteristics_ch1")) {
    parts <- gsub('"', "", strsplit(lines[i], "\t")[[1]][-1])
    if (length(parts) != ncol(dat)) next
    k <- sub(":.*$", "", parts[1]); v <- sub("^[^:]+: ?", "", parts)
    char[[k]] <- v
  }
}
print(names(char))
syms <- mapIds(hgu133plus2.db, keys = rownames(dat), keytype = "PROBEID", column = "SYMBOL", multiVals = "first")
ok <- !is.na(syms); e <- as.matrix(dat[ok, ]); rownames(e) <- syms[ok]
e <- e[!duplicated(rownames(e)), ]
## 临床：OS 时间/状态（GSE26193: "overall survival in months"? / "vital status"）
keys_ch <- list("time" = c("overall survival", "survival", "time"), "status" = c("vital status", "status", "event"))
grp <- rep(NA, ncol(e)); tm <- rep(NA, ncol(e))
for (k in names(char)) {
  kl <- tolower(k)
  if (grepl("os event|event", kl)) grp <- as.numeric(char[[k]])
  if (grepl("os time|time", kl)) tm <- as.numeric(char[[k]]) * 365.25
}
cat("e:", dim(e), " status:", table(grp, useNA = "ifany"), " time nonNA:", sum(!is.na(tm)), "\n")
clin <- data.frame(sample = colnames(e), futime = tm, fustat = grp)
fwrite(data.table(gene = rownames(e), e), "D:/子宫内膜异位症2/analysis/L0_bulk/results/track2/mime/GSE26193_exp.csv.gz")
fwrite(clin, "D:/子宫内膜异位症2/analysis/L0_bulk/results/track2/mime/GSE26193_clin.csv")
cat("DONE 45e\n")
