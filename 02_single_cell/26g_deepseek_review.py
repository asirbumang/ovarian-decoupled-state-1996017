# -*- coding: utf-8 -*-
"""
GSE266577（Cancer Cell 2024 官方标注）注释三审 第一轮：DeepSeek ×2
Call A: 标签-标记一致性（canonical panel）
Call B: 粒度/合并拆分建议（top markers 特异性）
"""
import csv, json, sys, time
import requests

sys.stdout.reconfigure(encoding="utf-8")
CFG = json.load(open(r"D:\子宫内膜异位症2\scripts\deepseek_config.json", encoding="utf-8"))
URL = CFG["base_url"].rstrip("/") + "/chat/completions"
HEADERS = {"Authorization": f"Bearer {CFG['api_key']}", "Content-Type": "application/json"}
OUT = r"D:\子宫内膜异位症2\analysis\L1_scRNA\results\annotation_review"

def call_api(messages, max_tokens=4000, temperature=0.1, retries=5):
    body = {"model": CFG["model"], "messages": messages,
            "temperature": temperature, "max_tokens": max_tokens}
    wait = 2
    for i in range(retries):
        try:
            r = requests.post(URL, headers=HEADERS, json=body, timeout=600)
        except Exception as e:
            if i == retries - 1:
                raise
            print("  retry", i, repr(e)[:80]); time.sleep(wait); wait *= 2; continue
        if r.status_code == 200:
            return r.json()["choices"][0]["message"]["content"]
        print("  api err", r.status_code, r.text[:150]); time.sleep(wait); wait *= 2
    raise RuntimeError("API failed")

# ---- 载入表 ----
panel = []   # ct, gene, mean, pct
with open(rf"{OUT}\GSE266577_canonical_panel.csv", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        panel.append(f"{row['cell_type']}|{row['gene']}|{float(row['mean_cpm']):.1f}|{float(row['pct_in']):.2f}")
top = []
with open(rf"{OUT}\GSE266577_top_markers_per_type.csv", encoding="utf-8") as f:
    rows = list(csv.DictReader(f))
# 每类型 top8
seen = {}
for row in rows:
    seen.setdefault(row["ct"], []).append(f"{row['gene']}({float(row['sdiff']):.2f})")
top8 = "\n".join(f"{ct}: " + ", ".join(v[:8]) for ct, v in seen.items())

CT_LIST = ", ".join(sorted(seen))

prompt_base = (
    "你是卵巢癌(高级别浆液性/HGSC) 单细胞转录组注释专家。以下来自 10x scRNA-seq 数据集 "
    "GSE266577（Cancer Cell 2024, Launonen et al.；137,060 细胞；29 例 HGSC，"
    "化疗前 vs 手术后(IDS)配对）的官方 22 类注释。请严格依标记证据评估，不要臆造。"

    f"\n\n细胞类型清单: {CT_LIST}"
)

prompt_A = prompt_base + (
    "\n\n【Call A: 一致性核查】下面是'规范标记基因 × 细胞类型'表达表 (mean_cpm=该类型平均 CPM 近似, pct_in=该类型内表达细胞比例):\n"
    + "\n".join(panel[:]) +
    "\n\n请逐类型给出: ①该标签与标记是否一致（一致/部分一致/不一致）; ②若不一致, 更合理的类型是什么; "
    "③每个类型给出置信度(高/中/低)。最后列出: 需要人工复核的类型+原因。"
    "格式: 类型 | 判定 | 合理标签 | 置信度。中文回复，不超过 1200 字。"
)

prompt_B = prompt_base + (
    "\n\n【Call B: 粒度核查】下面是每类型 top8 特异性 marker (括号内为与次高类型的 logCPM 差):\n"
    + top8 +
    "\n\n请评估: ①是否有多类型实为同一种群(该合并)? ②是否有单类型混有多个种群(该拆分)? "
    "③NK 亚群(CD16±)、DC 亚群(DC1/DC2/pDC/迁移DC)、T 亚群(Tem/Trm vs Tcm/Naive vs Treg vs Th17)、"
    "B 亚群(记忆/幼稚/浆细胞) 的划分是否可信? "
    "④若我们要做'肿瘤微环境组成分析'，建议的最终谱系层标签集(≤12 类)是什么? "
    "中文回复，不超过 1000 字。"
)

results = {}
for tag, prompt in (("A_consistency", prompt_A), ("B_granularity", prompt_B)):
    print(f"== DeepSeek call {tag} ==", flush=True)
    out = call_api([
        {"role": "system", "content": "你是严谨的 scRNA-seq 注释审稿人，只依据给出的证据，明确指出可疑点。"},
        {"role": "user", "content": prompt}])
    results[tag] = out
    with open(rf"{OUT}\deepseek_review_{tag}.md", "w", encoding="utf-8") as f:
        f.write(out)
    print(out[:400], "\n", flush=True)
    time.sleep(2)

print("DONE deepseek review x2")
