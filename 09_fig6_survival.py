#!/usr/bin/env python3
"""
Step 17 — Fig 6: bulk survival validation (GSE22762, 151 samples).
A: TNFSF9 KM, B: CD86 KM, C: CD200 KM, D: Cox HR forest plot.
Dataset source (GSE22762) labeled on the figure.
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines import CoxPHFitter

OUT = "/home/zhh2/CLL-ScRNA-seqdatasets/NK_analysis"
plt.rcParams.update({"font.size": 9, "axes.labelsize": 10, "axes.labelweight": "bold",
                     "axes.titleweight": "bold", "xtick.labelsize": 8, "ytick.labelsize": 8,
                     "legend.fontsize": 8, "figure.dpi": 200, "pdf.fonttype": 42})

m = pd.read_csv(f"{OUT}/results/GSE22762_merged_151.tsv", sep="\t")
print(f"n={len(m)}, events={m['event'].sum()}")

def add_panel(ax, L): ax.text(-0.16, 1.12, L, transform=ax.transAxes, fontsize=16, fontweight="bold", va="top", ha="left")

fig, axes = plt.subplots(1, 4, figsize=(15, 3.8))

# A/B/C: KM curves for TNFSF9, CD86, CD200
km_info = [("TNFSF9", "TNFSF9 (4-1BBL)", "#c0392b"),
           ("CD86", "CD86", "#1f6fb2"),
           ("CD200", "CD200", "#7b2e3e")]
for i, (gene, title, col) in enumerate(km_info):
    ax = axes[i]
    z = m[gene + "_z"].astype(float)
    hi = z >= z.median()
    lr = logrank_test(m["survival_time"][hi], m["survival_time"][~hi],
                      m["event"][hi], m["event"][~hi])
    kmf = KaplanMeierFitter()
    for grp, lbl, c in [(True, "high", "#c0392b"), (False, "low", "#2e7d32")]:
        kmf.fit(m["survival_time"][hi == grp], m["event"][hi == grp],
                label=f"{lbl} (n={(hi==grp).sum()})")
        kmf.plot_survival_function(ax=ax, color=c)
    ax.set_xlabel("OS (days)"); ax.set_ylabel("survival probability")
    ax.set_ylim(0, 1.05)
    ax.set_title(f"{title}\nlog-rank p={lr.p_value:.3f}")
    ax.legend(frameon=False)
    add_panel(ax, "ABCD"[i])

# D: Cox HR forest plot
ax = axes[3]
genes = ["TNFSF9", "CD86", "CD70", "CD200"]
hrs, lows, highs, ps = [], [], [], []
for g in genes:
    cph = CoxPHFitter()
    d = m[["survival_time", "event", g + "_z"]].rename(columns={g + "_z": "x"})
    cph.fit(d, "survival_time", "event")
    hr = float(cph.hazard_ratios_["x"])
    ci = cph.confidence_intervals_.loc["x"]
    lo = float(ci.iloc[0]); hi = float(ci.iloc[1])
    hrs.append(hr); lows.append(lo); highs.append(hi); ps.append(float(cph.summary.loc["x","p"]))
y = np.arange(len(genes))[::-1]
for i, (hr, lo, hi, p) in enumerate(zip(hrs, lows, highs, ps)):
    col = "#c0392b" if p < 0.05 else "#888"
    ax.plot([lo, hi], [y[i], y[i]], color=col, lw=2)
    ax.plot(hr, y[i], "o", color=col, ms=6)
    sig = "***" if p < 0.001 else ("**" if p < 0.01 else ("*" if p < 0.05 else "n.s."))
    ax.text(2.55, y[i], f"HR={hr:.2f} {sig}", va="center", fontsize=8, color=col, ha="left")
ax.axvline(1.0, color="k", ls="--", lw=0.8)
ax.set_yticks(y); ax.set_yticklabels(genes, fontweight="bold")
ax.set_xlabel("Cox hazard ratio (per SD)")
ax.set_ylabel("Gene")
ax.set_title("Cox univariate (OS)")
ax.set_xlim(0.3, 3.8)
add_panel(ax, "D")

plt.tight_layout()
fig.text(0.99, 0.01, "GSE22762", ha="right", va="bottom", fontsize=7, color="gray")
plt.savefig(f"{OUT}/figures/F6_survival.png", bbox_inches="tight")
plt.savefig(f"{OUT}/figures/F6_survival.pdf", bbox_inches="tight")
plt.close()
print("saved F6_survival.png/.pdf")
