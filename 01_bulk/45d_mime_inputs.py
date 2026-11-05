# -*- coding: utf-8 -*-
"""45d_mime_inputs.py: Mime 输入准备（TCGA-OV 训练 + GSE26193 测试；基因×样本 + futime/fustat）。"""
import numpy as np, pandas as pd, sys, os, gzip
sys.stdout.reconfigure(encoding='utf-8')
OUT = r"D:\子宫内膜异位症2\analysis\L0_bulk\results\track2\mime"
os.makedirs(OUT, exist_ok=True)

## ---- TCGA-OV 训练 ----
mx = pd.read_csv(r"D:\子宫内膜异位症2\data_omics\tcga_ov_counts_matrix.csv.gz", index_col=0)
clin = pd.read_csv(r"D:\子宫内膜异位症2\data_omics\tcga_ov_clinical.csv")
clin["sample"] = clin["barcode"].str.replace(r"-\d{2}[A-Z]$", "", regex=True)
clin["futime"] = np.where(clin["vital_status"] == "Dead", clin["days_to_death"], clin["days_to_last_followup"])
clin["fustat"] = (clin["vital_status"] == "Dead").astype(int)
clin = clin[clin["futime"].notna() & (clin["futime"] > 0)]
samples = [c for c in clin["sample"] if c in mx.columns]
E = mx[samples]
E = E.apply(lambda r: np.log2(r / r.sum() * 1e6 + 1), axis=1)
C = clin.set_index("sample").loc[samples, ["futime", "fustat"]]
E.to_csv(os.path.join(OUT, "tcga_exp.csv.gz"))
C.to_csv(os.path.join(OUT, "tcga_clin.csv"))
print("TCGA train:", E.shape, "events:", int(C["fustat"].sum()), flush=True)

## ---- GSE26193 测试（GPL570 log2；03b 同源加载）----
def parse_sm(gz):
    with gzip.open(gz, "rt", encoding="utf-8", errors="replace") as fh:
        lines = fh.readlines()
    beg = [i for i, l in enumerate(lines) if l.startswith("!series_matrix_table_begin")][0]
    dat = pd.read_csv(gz, sep="\t", skiprows=beg + 1, header=None, index_col=0)
    dat = dat.drop([])
    return dat
## 用 03b 已保存的表达/临床（若存在）
p1 = r"D:\子宫内膜异位症2\analysis\L0_bulk\results\survival_meta\GSE26193_expr.csv.gz"
p2 = r"D:\子宫内膜异位症2\analysis\L0_bulk\results\survival_meta\GSE26193_clin.csv"
if os.path.exists(p1):
    E2 = pd.read_csv(p1, index_col=0)
    C2 = pd.read_csv(p2, index_col=0)
    print("GSE26193 prep from saved:", E2.shape, list(C2.columns), flush=True)
else:
    print("saved GSE26193 not found; preparing from series matrix...", flush=True)
    import gzip as gz2
    import io
    with gz2.open(r"D:\子宫内膜异位症2\data_bulk\GSE26193_series_matrix.txt.gz", "rt", encoding="utf-8", errors="replace") as fh:
        txt = fh.read().splitlines()
    beg = [i for i, l in enumerate(txt) if l.startswith("!series_matrix_table_begin")][0]
    end = [i for i, l in enumerate(txt) if l.startswith("!series_matrix_table_end")][0]
    hdr = txt[beg + 1].replace('"', "").split("\t")
    df = pd.read_csv(io.StringIO("\n".join(txt[beg + 2:end - 1])), sep="\t", index_col=0, header=None, quotechar='"')
    df.columns = hdr[1:]
    ## 探针→symbol（hgu133plus2）
    import warnings
    warnings.filterwarnings("ignore")
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    pandas2ri.activate()
    ro.r('suppressMessages({library(AnnotationDbi); library(hgu133plus2.db)})')
    syms = ro.r('mapIds(hgu133plus2.db, keys=rownames(exp_df), keytype="PROBEID", column="SYMBOL", multiVals="first")')
    print("mapper ok")
print("DONE 45d")
