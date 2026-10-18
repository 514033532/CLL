#!/usr/bin/env python3
"""
Step 18c — GSE39671 (130 CLL) time-to-first-treatment survival replication.
Genes: TNFSF9, CD86, CD200, CD70 (+ ZAP70 as positive-control pipeline check).
Endpoint: TTFT (years), event=1 treated.
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test

OUT = "/home/zhh2/CLL-ScRNA-seqdatasets/NK_analysis"

expr = pd.read_csv(f"{OUT}/data/GSE39671_expr.tsv", sep="\t")
clin = pd.read_csv(f"{OUT}/data/GSE39671_clinical.tsv", sep="\t")

m = clin.merge(expr, on="gsm", how="inner")
print(f"merged n={len(m)}, events={m['event'].sum()}")

GENES = ["TNFSF9", "CD86", "CD200", "CD70", "ZAP70"]

# z-score each gene (within this cohort) for comparable Cox HRs
rows = []
for g in GENES:
    m[g + "_z"] = (m[g] - m[g].mean()) / m[g].std(ddof=0)

for g in GENES:
    z = m[g + "_z"].astype(float)
    hi = z >= z.median()
    t = m["ttft_years"].astype(float)
    e = m["event"].astype(int)
    lr = logrank_test(t[hi], t[~hi], e[hi], e[~hi])
    cph = CoxPHFitter()
    d = m[["ttft_years", "event", g + "_z"]].rename(columns={g + "_z": "x"})
    cph.fit(d, "ttft_years", "event")
    hr = float(cph.hazard_ratios_["x"])
    lo, hi2 = float(cph.confidence_intervals_.loc["x"].iloc[0]), float(cph.confidence_intervals_.loc["x"].iloc[1])
    p = float(cph.summary.loc["x", "p"])
    # direction: median TTFT in high vs low group
    rows.append(dict(gene=g, logrank_p=lr.p_value, cox_HR=hr, cox_CI_low=lo,
                     cox_CI_high=hi2, cox_p=p,
                     med_ttft_high=t[hi].median(), med_ttft_low=t[~hi].median()))

res = pd.DataFrame(rows)
res.to_csv(f"{OUT}/results/GSE39671_survival_ttft.tsv", sep="\t", index=False)
print(res.round(4).to_string(index=False))

# merged table for reproducibility
m.to_csv(f"{OUT}/results/GSE39671_merged.tsv", sep="\t", index=False)

# ---- figure: KM for the 4 core genes + ZAP70 ----
fig, axes = plt.subplots(1, 5, figsize=(18, 3.6))
for i, g in enumerate(GENES):
    ax = axes[i]
    z = m[g + "_z"].astype(float)
    hi = z >= z.median()
    t = m["ttft_years"].astype(float); e = m["event"].astype(int)
    lr = logrank_test(t[hi], t[~hi], e[hi], e[~hi])
    kmf = KaplanMeierFitter()
    for grp, lbl, c in [(True, "high", "#c0392b"), (False, "low", "#2e7d32")]:
        kmf.fit(t[hi == grp], e[hi == grp], label=f"{lbl} (n={(hi == grp).sum()})")
        kmf.plot_survival_function(ax=ax, color=c)
    ax.set_xlabel("TTFT (years)"); ax.set_ylabel("treatment-free prob.")
    ax.set_ylim(0, 1.05)
    ax.set_title(f"{g}\nlog-rank p={lr.p_value:.3f}")
    ax.legend(frameon=False, fontsize=7)
    ax.text(-0.16, 1.12, "ABCDE"[i], transform=ax.transAxes, fontsize=14, fontweight="bold")
plt.tight_layout()
fig.text(0.99, 0.01, "GSE39671 (UCSD CLL, n=130, TTFT)", ha="right", va="bottom", fontsize=7, color="gray")
plt.savefig(f"{OUT}/figures/GSE39671_survival.png", bbox_inches="tight", dpi=200)
plt.savefig(f"{OUT}/figures/GSE39671_survival.pdf", bbox_inches="tight")
plt.close()
print("saved GSE39671_survival.png/.pdf + results/GSE39671_survival_ttft.tsv")
