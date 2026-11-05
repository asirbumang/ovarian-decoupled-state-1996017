# -*- coding: utf-8 -*-
import pandas as pd

d = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L0_bulk\results\track0_tcga_scores_states.csv")
p = pd.read_csv(r"D:\子宫内膜异位症2\analysis\L0_bulk\results\track0_tcga_discriminator_partial.csv")
print("n samples:", len(d))
print("strict:", int(d.state_high.sum()), " relax:", int(d.state_relax.sum()),
      " pct:", round(100 * d.state_high.mean(), 3), "%")
print("tier counts:", dict(d.decouple_tier.value_counts()))
print("response counts:", dict(d.outcome_bin.value_counts(dropna=False)))
print("decouple median CR/PD/PR/SD:")
for g in ["CR", "PR", "SD", "PD"]:
    print("  ", g, round(d.loc[d.outcome_bin == g, "decouple"].median(), 3))
print("partial table:")
print(p.to_string(index=False))
# 与图对照：figure says p=0.019 / 0.65 / 0.0025 / 0.00065 与 partial_p 一致
print("partial_p rounded:", [f"{x:.4g}" for x in p.partial_p])
print("partial_r rounded:", [f"{x:.3f}" for x in p.partial_r])
# OS 队列 tier 计数
od = d[(d.os_days.notna()) & (d.os_days > 0)]
print("OS cohort n:", len(od), " tier:", dict(od.decouple_tier.value_counts()))
