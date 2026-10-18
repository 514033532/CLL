#!/usr/bin/env python3
"""
Reframed figures for TNFSF9-focused manuscript:
Fig1 (prognosis, core): TNFSF9/CD86/CD200 KM + Cox forest (GSE22762, 151)
Fig2 (spatial + program): TNFSF9/CD200/CD70 PB-LN + program score (GSE161610)
Fig3 (contrast + cross-cancer): CD200 progression + CLL vs CML
"""
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test
import matplotlib.patches as mpatches

OUT = "/home/zhh2/CLL-ScRNA-seqdatasets/NK_analysis"
plt.rcParams.update({"font.size": 9, "axes.labelsize": 10, "axes.labelweight": "bold",
                     "axes.titleweight": "bold", "xtick.labelsize": 8, "ytick.labelsize": 8,
                     "legend.fontsize": 8, "figure.dpi": 200, "pdf.fonttype": 42})
TISSUE_COLORS = {"PB": "#c9d4e3", "LN": "#1f3a5f"}

def add_panel(ax, L): ax.text(-0.16, 1.12, L, transform=ax.transAxes, fontsize=16, fontweight="bold", va="top", ha="left")
def style(ax, x=None, y=None, t=None):
    if x: ax.set_xlabel(x)
    if y: ax.set_ylabel(y)
    if t: ax.set_title(t)
    for s in ["top", "right"]: ax.spines[s].set_visible(False)

m = pd.read_csv(f"{OUT}/results/GSE22762_merged_151.tsv", sep="\t")
g16 = sc.read_h5ad(f"{OUT}/data/gse161610_processed.h5ad")
def col(g): return np.asarray(g16.layers["log1p"][:, list(g16.var_names).index(g)].toarray()).ravel()
b16 = (g16.obs["celltype"] == "B").values

# ================= Fig 1: prognosis (core) =================================
fig, axes = plt.subplots(1, 4, figsize=(16, 4))
for i, (gene, title) in enumerate([("TNFSF9", "TNFSF9 (4-1BBL)"),
                                   ("CD86", "CD86"),
                                   ("CD200", "CD200")]):
    ax = axes[i]
    z = m[gene + "_z"].astype(float)
    hi = z >= z.median()
    lr = logrank_test(m["survival_time"][hi], m["survival_time"][~hi], m["event"][hi], m["event"][~hi])
    kmf = KaplanMeierFitter()
    for grp, lbl, c in [(True, "high", "#c0392b"), (False, "low", "#2e7d32")]:
        kmf.fit(m["survival_time"][hi == grp], m["event"][hi == grp], label=f"{lbl} (n={(hi==grp).sum()})")
        kmf.plot_survival_function(ax=ax, color=c)
    ax.set_xlabel("OS (days)"); ax.set_ylabel("survival"); ax.set_ylim(0, 1.05)
    ax.set_title(f"{title}\nlog-rank p={lr.p_value:.4f}")
    ax.legend(frameon=False)
    add_panel(ax, "ABC"[i])
# D: forest plot
ax = axes[3]
genes = ["TNFSF9", "CD86", "CD70", "CD200"]
hrs, lows, highs, ps = [], [], [], []
for g in genes:
    cph = CoxPHFitter()
    d = m[["survival_time", "event", g + "_z"]].rename(columns={g + "_z": "x"})
    cph.fit(d, "survival_time", "event")
    hr = float(cph.hazard_ratios_["x"]); ci = cph.confidence_intervals_.loc["x"]
    hrs.append(hr); lows.append(float(ci.iloc[0])); highs.append(float(ci.iloc[1])); ps.append(float(cph.summary.loc["x","p"]))
y = np.arange(len(genes))[::-1]
for i, (hr, lo, hi, p) in enumerate(zip(hrs, lows, highs, ps)):
    col_ = "#c0392b" if p < 0.05 else "#888"
    ax.plot([lo, hi], [y[i], y[i]], color=col_, lw=2)
    ax.plot(hr, y[i], "o", color=col_, ms=6)
    sig = "***" if p < 0.001 else ("**" if p < 0.01 else ("*" if p < 0.05 else "n.s."))
    ax.text(2.5, y[i], f"HR={hr:.2f} {sig}", va="center", fontsize=8, color=col_, ha="left")
ax.axvline(1.0, color="k", ls="--", lw=0.8)
ax.set_yticks(y); ax.set_yticklabels(genes, fontweight="bold")
ax.set_xlabel("Cox HR (per SD)"); ax.set_ylabel("Gene"); ax.set_title("Cox univariate (OS)")
ax.set_xlim(0.3, 3.6)
add_panel(ax, "D")
fig.text(0.99, 0.01, "GSE22762 (n=151)", ha="right", va="bottom", fontsize=7, color="gray")
plt.tight_layout(); plt.savefig(f"{OUT}/figures/TF1_prognosis.png", bbox_inches="tight")
plt.savefig(f"{OUT}/figures/TF1_prognosis.pdf", bbox_inches="tight"); plt.close()

# ================= Fig 2: spatial + program ================================
fig, axes = plt.subplots(1, 4, figsize=(14, 4))
for i, g in enumerate(["TNFSF9", "CD200", "CD70"]):
    ax = axes[i]
    data = [col(g)[b16 & (g16.obs["tissue"] == "PB").values],
            col(g)[b16 & (g16.obs["tissue"] == "LN").values]]
    vp = ax.violinplot(data, positions=[0, 1], showmeans=True, showextrema=False, widths=0.7)
    for pc, c in zip(vp["bodies"], [TISSUE_COLORS["PB"], TISSUE_COLORS["LN"]]):
        pc.set_facecolor(c); pc.set_alpha(0.75); pc.set_edgecolor("k")
    ax.set_xticks([0, 1]); ax.set_xticklabels(["PB", "LN"], fontweight="bold")
    style(ax, "Tissue", f"{g} (log1p)", g)
    if i == 0:
        ax.legend([mpatches.Patch(color=TISSUE_COLORS["PB"]), mpatches.Patch(color=TISSUE_COLORS["LN"])],
                  ["PB", "LN"], frameon=False, loc="upper left")
    add_panel(ax, "ABC"[i])
# D: program score PB vs LN
ax = axes[3]
prog_idx = [list(g16.var_names).index(g) for g in ["CD200", "CD70", "TNFSF9"]]
prog = np.asarray(g16.layers["log1p"][:, prog_idx].mean(axis=1)).ravel()
data = [prog[b16 & (g16.obs["tissue"] == "PB").values], prog[b16 & (g16.obs["tissue"] == "LN").values]]
vp = ax.violinplot(data, positions=[0, 1], showmeans=True, showextrema=False, widths=0.7)
for pc, c in zip(vp["bodies"], [TISSUE_COLORS["PB"], TISSUE_COLORS["LN"]]):
    pc.set_facecolor(c); pc.set_alpha(0.75); pc.set_edgecolor("k")
ax.set_xticks([0, 1]); ax.set_xticklabels(["PB", "LN"], fontweight="bold")
style(ax, "Tissue", "program score", "CD200+CD70+TNFSF9")
add_panel(ax, "D")
fig.text(0.99, 0.01, "GSE161610 (5 patients, paired PB/LN)", ha="right", va="bottom", fontsize=7, color="gray")
plt.tight_layout(); plt.savefig(f"{OUT}/figures/TF2_spatial.png", bbox_inches="tight")
plt.savefig(f"{OUT}/figures/TF2_spatial.pdf", bbox_inches="tight"); plt.close()

# ================= Fig 3: contrast + cross-cancer ==========================
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
# A: CD200 progression (diagnostic, not prognostic)
ax = axes[0]
stages = ["Reactive\nLN", "FL\nLN", "CLL\nPB", "CLL\nLN", "Richter"]
vals = [0.065, 0.027, 0.132, 0.341, 0.367]
cols = ["#9aa0a6", "#17a2b8", "#6cb2eb", "#d84b6a", "#7b2e3e"]
ax.bar(range(5), vals, color=cols, alpha=0.85, edgecolor="k", lw=0.5)
for i, v in enumerate(vals): ax.text(i, v + 0.012, f"{v:.2f}", ha="center", fontsize=8)
ax.set_xticks(range(5)); ax.set_xticklabels(stages, fontsize=8)
style(ax, "Disease stage", "CD200 (mean log1p)", "CD200: diagnostic, rises with progression")
add_panel(ax, "A")
# B: cross-cancer CLL vs CML
ax = axes[1]
ligands = ["TNFSF9", "CD200", "CD70", "PVR"]
cll = [0.396, 0.341, 0.140, 0.000]
cml = [0.074, 0.000, 0.009, 0.031]
x = np.arange(len(ligands)); w = 0.38
ax.bar(x - w/2, cll, w, label="CLL B cells", color="#d84b6a", edgecolor="k", lw=0.5)
ax.bar(x + w/2, cml, w, label="CML tumor", color="#1f6fb2", edgecolor="k", lw=0.5)
ax.set_xticks(x); ax.set_xticklabels(ligands, fontweight="bold")
ax.set_xlabel("Ligand"); ax.set_ylabel("mean expression")
ax.set_title("Immune-evasion ligands: CLL vs CML")
ax.legend(frameon=False, loc="upper right")
add_panel(ax, "B")
fig.text(0.99, 0.01, "A: multi-cohort | B: GSE161610 vs PMID 37919606", ha="right", va="bottom", fontsize=7, color="gray")
plt.tight_layout(); plt.savefig(f"{OUT}/figures/TF3_contrast.png", bbox_inches="tight")
plt.savefig(f"{OUT}/figures/TF3_contrast.pdf", bbox_inches="tight"); plt.close()

print("saved TF1/TF2/TF3")
