#!/usr/bin/env python3
"""Final GSE39671 survival figure: 4 core genes + LPL positive control."""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

OUT = "/home/zhh2/CLL-ScRNA-seqdatasets/NK_analysis"
expr = pd.read_csv(f"{OUT}/data/GSE39671_expr.tsv", sep="\t")
clin = pd.read_csv(f"{OUT}/data/GSE39671_clinical.tsv", sep="\t")
m = clin.merge(expr, on="gsm", how="inner")
t = m["ttft_years"].astype(float); e = m["event"].astype(int)

GENES = ["TNFSF9", "CD86", "CD200", "CD70", "LPL"]
TITLES = {"TNFSF9": "TNFSF9 (4-1BBL)", "CD86": "CD86", "CD200": "CD200",
          "CD70": "CD70", "LPL": "LPL (positive control)"}

fig, axes = plt.subplots(1, 5, figsize=(19, 3.6))
for i, g in enumerate(GENES):
    ax = axes[i]
    z = (m[g] - m[g].mean()) / m[g].std(ddof=0)
    hi = z >= z.median()
    lr = logrank_test(t[hi], t[~hi], e[hi], e[~hi])
    kmf = KaplanMeierFitter()
    for grp, lbl, c in [(True, "high", "#c0392b"), (False, "low", "#2e7d32")]:
        kmf.fit(t[hi == grp], e[hi == grp], label=f"{lbl} (n={(hi == grp).sum()})")
        kmf.plot_survival_function(ax=ax, color=c)
    ax.set_xlabel("TTFT (years)"); ax.set_ylabel("treatment-free prob.")
    ax.set_ylim(0, 1.05)
    ax.set_title(f"{TITLES[g]}\nlog-rank p={lr.p_value:.3f}")
    ax.legend(frameon=False, fontsize=7)
    ax.text(-0.16, 1.12, "ABCDE"[i], transform=ax.transAxes, fontsize=14, fontweight="bold")
plt.tight_layout()
fig.text(0.99, 0.01, "GSE39671 (UCSD CLL, n=130, TTFT endpoint)", ha="right", va="bottom", fontsize=7, color="gray")
plt.savefig(f"{OUT}/figures/GSE39671_survival.png", bbox_inches="tight", dpi=200)
plt.savefig(f"{OUT}/figures/GSE39671_survival.pdf", bbox_inches="tight")
plt.close()
print("saved final GSE39671_survival figure")

# combined results table (core + controls)
rows = []
for g in ["TNFSF9", "CD86", "CD200", "CD70", "ZAP70", "LPL", "CD38"]:
    z = (m[g] - m[g].mean()) / m[g].std(ddof=0)
    hi = z >= z.median()
    lr = logrank_test(t[hi], t[~hi], e[hi], e[~hi])
    from lifelines import CoxPHFitter
    cph = CoxPHFitter(); d = m[["ttft_years", "event"]].copy(); d["x"] = z.values
    cph.fit(d, "ttft_years", "event")
    rows.append(dict(gene=g, role="control" if g in ("ZAP70", "LPL", "CD38") else "core",
                     logrank_p=round(lr.p_value, 4),
                     cox_HR=round(float(cph.hazard_ratios_["x"]), 3),
                     cox_p=round(float(cph.summary.loc["x", "p"]), 4),
                     med_ttft_high=round(t[hi].median(), 2), med_ttft_low=round(t[~hi].median(), 2)))
pd.DataFrame(rows).to_csv(f"{OUT}/results/GSE39671_survival_ttft.tsv", sep="\t", index=False)
print("saved results/GSE39671_survival_ttft.tsv")
