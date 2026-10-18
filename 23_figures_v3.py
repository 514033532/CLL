#!/usr/bin/env python3
"""v3 figures with GEO source labels: Fig2 (CD70-CD27 permutation), Fig3 (CD200R1), Fig5 (survival replication)."""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test

OUT = "/home/zhh2/CLL-ScRNA-seqdatasets/NK_analysis"
plt.rcParams.update({"font.size": 9, "axes.labelsize": 10, "axes.labelweight": "bold",
                     "axes.titleweight": "bold", "xtick.labelsize": 8, "ytick.labelsize": 8,
                     "legend.fontsize": 8, "figure.dpi": 200, "pdf.fonttype": 42})

def panel(ax, L): ax.text(-0.18, 1.12, L, transform=ax.transAxes, fontsize=15, fontweight="bold", va="top")
def style(ax, x=None, y=None, t=None):
    if x: ax.set_xlabel(x)
    if y: ax.set_ylabel(y)
    if t: ax.set_title(t, fontsize=10)
    for s in ["top", "right"]: ax.spines[s].set_visible(False)

# ================= Fig 2: CD70-CD27 permutation test =================
fig, ax = plt.subplots(figsize=(6, 3.8))
interactions = ["B→B", "B→CD4 T", "B→CD8 T", "B→NK", "NK→CD4 T"]
lr_means = [0.534, 0.541, 0.458, 0.175, 0.600]
pvals = [0.000, 0.177, 0.887, 1.000, 0.101]
colors = ["#c0392b" if p < 0.05 else "#b0b0b0" for p in pvals]
x = np.arange(len(interactions))
bars = ax.bar(x, lr_means, color=colors, edgecolor="k", lw=0.6, width=0.62)
for i, (v, p) in enumerate(zip(lr_means, pvals)):
    sig = "p=%.0e" % p if p < 0.001 else f"p={p:.3f}"
    ax.text(i, v + 0.02, sig, ha="center", fontsize=7.5, color="#333")
ax.set_xticks(x); ax.set_xticklabels(interactions, fontsize=8)
ax.set_ylabel("CD70→CD27 interaction score (lr_means)")
ax.set_title("CD70–CD27 in CLL LN: only B→B is significant")
ax.axhline(0, color="k", lw=0.8)
# legend for significance
import matplotlib.patches as mpatches
ax.legend(handles=[mpatches.Patch(color="#c0392b", label="significant (p<0.05)"),
                   mpatches.Patch(color="#b0b0b0", label="not significant")], frameon=False, loc="upper right", fontsize=7.5)
fig.text(0.99, 0.01, "GSE161610 | CellPhoneDB, 1,000 permutations (OmniPath consensus)", ha="right", va="bottom", fontsize=6.5, color="gray")
plt.tight_layout()
plt.savefig(f"{OUT}/figures/Fig2_CD70CD27_permutation.png", bbox_inches="tight")
plt.savefig(f"{OUT}/figures/Fig2_CD70CD27_permutation.pdf", bbox_inches="tight")
plt.close()

# ================= Fig 3: CD200R1 minimal expression =================
fig, ax = plt.subplots(figsize=(6, 3.8))
ct = ["B", "CD4 T", "CD8 T", "NK", "Myeloid"]
vals = [0.002, 0.028, 0.027, 0.000, 0.038]
ax.bar(ct, vals, color="#1f6fb2", edgecolor="k", lw=0.6, width=0.6)
for i, v in enumerate(vals):
    ax.text(i, v + 0.0015, f"{v:.3f}", ha="center", fontsize=8)
ax.set_ylabel("CD200R1 mean expression (log1p)")
ax.set_title("CD200R1 is essentially absent on CLL LN immune cells")
ax.set_ylim(0, 0.05)
fig.text(0.99, 0.01, "GSE161610 (lymph node)", ha="right", va="bottom", fontsize=6.5, color="gray")
plt.tight_layout()
plt.savefig(f"{OUT}/figures/Fig3_CD200R1_expression.png", bbox_inches="tight")
plt.savefig(f"{OUT}/figures/Fig3_CD200R1_expression.pdf", bbox_inches="tight")
plt.close()

# ================= Fig 5: survival replication =================
m22762 = pd.read_csv(f"{OUT}/results/GSE22762_merged_151.tsv", sep="\t")
m39671 = pd.read_csv(f"{OUT}/results/GSE39671_merged.tsv", sep="\t")
# add LPL_z from the expression file (has LPL, not in merged)
_expr = pd.read_csv(f"{OUT}/data/GSE39671_expr.tsv", sep="\t")[["gsm", "LPL"]]
m39671 = m39671.merge(_expr, on="gsm", how="left")
m39671["LPL_z"] = (m39671["LPL"] - m39671["LPL"].mean()) / m39671["LPL"].std(ddof=0)

fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.8))

# A: GSE22762 TNFSF9 KM (OS)
ax = axes[0]
z = m22762["TNFSF9_z"].astype(float); hi = z >= z.median()
t = m22762["survival_time"].astype(float); e = m22762["event"].astype(int)
lr = logrank_test(t[hi], t[~hi], e[hi], e[~hi])
kmf = KaplanMeierFitter()
for grp, lbl, c in [(True, "high", "#c0392b"), (False, "low", "#2e7d32")]:
    kmf.fit(t[hi == grp], e[hi == grp], label=f"{lbl} (n={(hi == grp).sum()})")
    kmf.plot_survival_function(ax=ax, color=c)
ax.set_xlabel("OS (days)"); ax.set_ylabel("survival prob."); ax.set_ylim(0, 1.05)
ax.set_title(f"GSE22762 TNFSF9\nlog-rank p={lr.p_value:.3f}")
ax.legend(frameon=False, fontsize=7); panel(ax, "A")

# B: GSE39671 TNFSF9 KM (TTFT)
ax = axes[1]
z = m39671["TNFSF9_z"].astype(float); hi = z >= z.median()
t = m39671["ttft_years"].astype(float); e = m39671["event"].astype(int)
lr = logrank_test(t[hi], t[~hi], e[hi], e[~hi])
kmf = KaplanMeierFitter()
for grp, lbl, c in [(True, "high", "#c0392b"), (False, "low", "#2e7d32")]:
    kmf.fit(t[hi == grp], e[hi == grp], label=f"{lbl} (n={(hi == grp).sum()})")
    kmf.plot_survival_function(ax=ax, color=c)
ax.set_xlabel("TTFT (years)"); ax.set_ylabel("treatment-free prob."); ax.set_ylim(0, 1.05)
ax.set_title(f"GSE39671 TNFSF9\nlog-rank p={lr.p_value:.3f}")
ax.legend(frameon=False, fontsize=7); panel(ax, "B")

# C: GSE39671 LPL (positive control)
ax = axes[2]
z = m39671["LPL_z"].astype(float); hi = z >= z.median()
t = m39671["ttft_years"].astype(float); e = m39671["event"].astype(int)
lr = logrank_test(t[hi], t[~hi], e[hi], e[~hi])
kmf = KaplanMeierFitter()
for grp, lbl, c in [(True, "high", "#c0392b"), (False, "low", "#2e7d32")]:
    kmf.fit(t[hi == grp], e[hi == grp], label=f"{lbl} (n={(hi == grp).sum()})")
    kmf.plot_survival_function(ax=ax, color=c)
ax.set_xlabel("TTFT (years)"); ax.set_ylabel("treatment-free prob."); ax.set_ylim(0, 1.05)
ax.set_title(f"GSE39671 LPL (control)\nlog-rank p={lr.p_value:.3f}")
ax.legend(frameon=False, fontsize=7); panel(ax, "C")

fig.text(0.99, 0.01, "A: GSE22762 (OS, n=151) | B,C: GSE39671 (TTFT, n=130)", ha="right", va="bottom", fontsize=6.5, color="gray")
plt.tight_layout()
plt.savefig(f"{OUT}/figures/Fig5_survival_replication.png", bbox_inches="tight")
plt.savefig(f"{OUT}/figures/Fig5_survival_replication.pdf", bbox_inches="tight")
plt.close()
print("saved Fig2/Fig3/Fig5")
