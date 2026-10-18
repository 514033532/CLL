#!/usr/bin/env python3
"""Rebuild 6 large multi-panel figures (大图套小图) with recovered data + GEO labels."""
import warnings; warnings.filterwarnings("ignore")
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

OUT = "/home/zhh2/CLL-ScRNA-seqdatasets/NK_analysis"
plt.rcParams.update({"font.size": 8.5, "axes.labelsize": 9, "axes.labelweight": "bold",
                     "axes.titleweight": "bold", "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
                     "legend.fontsize": 7, "figure.dpi": 300, "pdf.fonttype": 42, "ps.fonttype": 42,
                     "font.family": "sans-serif",
                     "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
                     "axes.linewidth": 1.0, "xtick.major.width": 1.0, "ytick.major.width": 1.0,
                     "lines.linewidth": 1.0, "patch.linewidth": 1.0})
CT_COLORS = {"B": "#3b7dd8", "CD4_T": "#57a05a", "CD8_T": "#e8912a", "NK": "#9b5de5",
             "Myeloid": "#d84b6a", "Other": "#9aa0a6"}
TISSUE_COLORS = {"PB": "#c9d4e3", "LN": "#1f3a5f"}

import os
JOUT = f"{OUT}/figures_journal"
os.makedirs(JOUT, exist_ok=True)
def save_all(fig, name):
    for ext in ["png", "pdf", "eps", "tif"]:
        kw = dict(dpi=300, pil_kwargs={"compression": "tiff_lzw"}) if ext == "tif" else {}
        fig.savefig(f"{JOUT}/{name}.{ext}", bbox_inches="tight", **kw)

def panel(ax, L): ax.text(-0.06, 1.16, L, transform=ax.transAxes, fontsize=12, fontweight="bold", va="top", ha="right")
def style(ax, x=None, y=None, t=None):
    if x: ax.set_xlabel(x)
    if y: ax.set_ylabel(y)
    if t: ax.set_title(t, fontsize=9)
    for s in ["top", "right"]: ax.spines[s].set_visible(False)

g16 = sc.read_h5ad(f"{OUT}/data/gse161610_processed.h5ad")
def col(d, g): return np.asarray(d.layers["log1p"][:, list(d.var_names).index(g)].toarray()).ravel()

# ============ Fig 1: Atlas + immune composition (4 panels) ============
fig, axes = plt.subplots(2, 2, figsize=(7.09, 6.30))
ax = axes[0,0]
for c in sorted(g16.obs["celltype"].unique()):
    m = g16.obs["celltype"]==c
    ax.scatter(g16.obsm["X_umap"][m,0], g16.obsm["X_umap"][m,1], s=1.5, alpha=0.5, color=CT_COLORS.get(c,"#9aa0a6"), rasterized=True, label=c)
handles = [plt.Line2D([0],[0], marker="o", color="w", markerfacecolor=CT_COLORS.get(c,"#9aa0a6"), markersize=5, label=c) for c in sorted(g16.obs["celltype"].unique())]
ax.legend(handles=handles, frameon=False, fontsize=6, loc="best", ncol=2, handletextpad=0.2, columnspacing=0.8)
style(ax, "UMAP1", "UMAP2", "Cell types"); panel(ax, "a")
ax = axes[0,1]
for t in ["PB","LN"]:
    m = g16.obs["tissue"]==t
    ax.scatter(g16.obsm["X_umap"][m,0], g16.obsm["X_umap"][m,1], s=1.5, alpha=0.5, color=TISSUE_COLORS[t], rasterized=True)
style(ax, "UMAP1", "UMAP2", "Tissue"); panel(ax, "b")
ax = axes[1,0]
prop = pd.read_csv(f"{OUT}/results/gse161610_proportions.csv", index_col=0)
prop = prop[["B","CD4_T","CD8_T","Myeloid","NK"]]
prop.plot(kind="bar", stacked=True, ax=ax, color=[CT_COLORS[c] for c in prop.columns], edgecolor="k", lw=1.0)
ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False, fontsize=6.5, title="cell type")
style(ax, "Tissue", "fraction", "Composition (GSE161610)"); panel(ax, "c")
ax = axes[1,1]
comp = pd.read_csv(f"{OUT}/results/immune_composition.csv")
sub = comp[comp.subset.isin(["NK_CD56dim","NK_CD56bright","NK_adaptive","Mono_classical","Mono_nonclassical","Mono_intermediate","cDC1","cDC2","pDC"])]
# healthy vs CLL PB vs CLL LN vs Richter
groups = [("GSE213516","healthy","Healthy"),("GSE161610","PB","CLL PB"),("GSE161610","LN","CLL LN"),("GSE205704","Richter","Richter")]
vals=[]
for ds,g,lab in groups:
    s = sub[(sub.dataset==ds)&(sub.group==g)]
    tot = s["total"].iloc[0] if len(s) else 1
    vals.append(s["n"].sum()/tot*100)
ax.bar(range(4), vals, color="#1f3a5f", edgecolor="k", lw=1.0)
for i,v in enumerate(vals): ax.text(i, v+0.2, f"{v:.1f}", ha="center", fontsize=7)
ax.set_xticks(range(4)); ax.set_xticklabels([g[2] for g in groups], fontsize=7)
style(ax, None, "% innate immune", "Innate compartment"); panel(ax, "d")
fig.text(0.99, 0.01, "GSE161610, GSE213516, GSE205704", ha="right", va="bottom", fontsize=6.5, color="gray")
plt.tight_layout()
save_all(fig, "Fig1_atlas_composition"); plt.close()

# ============ Fig 2: Tissue ligand landscape (4 panels) ============
b16 = (g16.obs["celltype"]=="B").values
fig, axes = plt.subplots(2, 2, figsize=(7.09, 5.51))
for i,(g,ttl) in enumerate([("CD200","CD200"),("CD70","CD70"),("TNFSF9","TNFSF9 (4-1BBL)"),("CD86","CD86")]):
    ax = axes[i//2, i%2]
    data = [col(g16,g)[b16 & (g16.obs["tissue"]=="PB").values], col(g16,g)[b16 & (g16.obs["tissue"]=="LN").values]]
    vp = ax.violinplot(data, positions=[0,1], showmeans=True, showextrema=False, widths=0.7)
    for pc,c in zip(vp["bodies"], [TISSUE_COLORS["PB"], TISSUE_COLORS["LN"]]): pc.set_facecolor(c); pc.set_alpha(0.75); pc.set_edgecolor("k")
    ax.set_xticks([0,1]); ax.set_xticklabels(["PB","LN"]); ax.set_ylabel(f"{g} (log1p)"); ax.set_title(ttl)
    for j,t in enumerate(["PB","LN"]):
        m = col(g16,g)[b16 & (g16.obs["tissue"]==t).values].mean()
        ax.text(j, m+0.05, f"{m:.2f}", ha="center", fontsize=7.5)
    panel(ax, "abcd"[i])
fig.text(0.99, 0.01, "GSE161610 (5 patients, paired PB/LN)", ha="right", va="bottom", fontsize=6.5, color="gray")
plt.tight_layout()
save_all(fig, "Fig2_ligand_landscape"); plt.close()

# ============ Fig 3: Receptor wiring (4 panels) ============
fig, axes = plt.subplots(2, 2, figsize=(7.09, 5.51))
# a: CD70-CD27 permutation
ax = axes[0,0]
interactions = ["B→B", "B→CD4 T", "B→CD8 T", "B→NK"]
lr = [0.534, 0.541, 0.458, 0.175]; pv = [0.000, 0.177, 0.887, 1.000]
colors = ["#c0392b" if p<0.05 else "#b0b0b0" for p in pv]
ax.bar(range(4), lr, color=colors, edgecolor="k", lw=1.0)
for i,(v,p) in enumerate(zip(lr,pv)):
    ax.text(i, v+0.02, f"p={p:.2f}", ha="center", fontsize=7)
ax.set_xticks(range(4)); ax.set_xticklabels(interactions, fontsize=7.5)
style(ax, None, "CD70→CD27 score", "CD70–CD27 (only B→B significant)"); panel(ax, "a")
# b: CD27 receptor expression
ax = axes[0,1]
ct_keys = ["B","CD4_T","CD8_T","NK","Myeloid"]; ct_lbl = ["B","CD4 T","CD8 T","NK","Myeloid"]
cd27_data = [col(g16,"CD27")[g16.obs["celltype"]==c] for c in ct_keys]
vp = ax.violinplot(cd27_data, showmeans=True, showextrema=False)
for pc in vp["bodies"]: pc.set_facecolor("#1f6fb2"); pc.set_alpha(0.55); pc.set_edgecolor("k"); pc.set_linewidth(1.0)
ax.set_xticks(range(1, len(ct_lbl)+1)); ax.set_xticklabels(ct_lbl, fontsize=7.5)
for i,d in enumerate(cd27_data): ax.text(i+1, d.mean()+0.03, f"{d.mean():.2f}", ha="center", fontsize=7)
style(ax, None, "CD27 (log1p)", "CD27 receptor expression (per cell)"); panel(ax, "b")
# c: CD200R1 expression (near zero)
ax = axes[1,0]
cd200r1 = [0.002,0.028,0.027,0.000,0.038]
ax.bar(ct_lbl, cd200r1, color="#d84b6a", edgecolor="k", lw=1.0)
for i,v in enumerate(cd200r1): ax.text(i, v+0.002, f"{v:.3f}", ha="center", fontsize=7)
ax.set_ylim(0, 0.05)
style(ax, None, "CD200R1 (log1p)", "CD200R1 essentially absent"); panel(ax, "c")
# d: CD200 ligand (high) vs CD200R1 receptor (absent) contrast
ax = axes[1,1]
cd200_lig = [col(g16,"CD200")[b16 & (g16.obs["tissue"]=="LN").values].mean()]
ax.bar(["CD200 (ligand, B LN)"], cd200_lig, color="#3b7dd8", edgecolor="k", lw=1.0)
ax.bar(["CD200R1 (receptor, max LN)"], [0.038], color="#d84b6a", edgecolor="k", lw=1.0)
ax.text(0, cd200_lig[0]+0.02, f"{cd200_lig[0]:.2f}", ha="center", fontsize=8)
ax.text(1, 0.038+0.01, "0.038", ha="center", fontsize=8)
style(ax, None, "expression (log1p)", "Ligand high, receptor absent"); panel(ax, "d")
fig.text(0.99, 0.01, "GSE161610 | CellPhoneDB, 1,000 permutations", ha="right", va="bottom", fontsize=6.5, color="gray")
plt.tight_layout()
save_all(fig, "Fig3_receptor_wiring"); plt.close()

# ============ Fig 4: Cross-cancer + disease gradient (2 panels) ============
fig, axes = plt.subplots(1, 2, figsize=(7.09, 2.58))
# a: 4-cancer ligand atlas
ax = axes[0]
ligands = ["CD200","CD70","TNFSF9","CD86","LGALS9","PVR"]
cll=[0.341,0.140,0.396,0.112,0.214,0.0003]; aml=[0.083,0.119,0.049,0.017,0.124,0.059]
cml=[0.000,0.009,0.074,0.000,0.024,0.031]; mm=[0.036,0.003,0.001,0.044,0.012,0.001]
x=np.arange(6); w=0.19
for i,(nm,vals,colr) in enumerate([("CLL",cll,"#d84b6a"),("AML",aml,"#e8912a"),("CML",cml,"#1f6fb2"),("MM",mm,"#17a2b8")]):
    ax.bar(x+(i-1.5)*w, vals, w, label=nm, color=colr, edgecolor="k", lw=1.0)
ax.set_xticks(x); ax.set_xticklabels(ligands, fontweight="bold"); ax.legend(frameon=False, fontsize=7)
style(ax, "Ligand", "mean log1p(CPM)", "Four-cancer ligand survey"); panel(ax, "a")
# b: CD200 disease-stage gradient
ax = axes[1]
stages=["Reactive\nLN","FL\nLN","CLL\nPB","CLL\nLN","Richter"]; vals=[0.065,0.027,0.132,0.341,0.367]
ax.bar(range(5), vals, color=["#9aa0a6","#17a2b8","#6cb2eb","#d84b6a","#7b2e3e"], edgecolor="k", lw=1.0)
for i,v in enumerate(vals): ax.text(i, v+0.01, f"{v:.2f}", ha="center", fontsize=7.5)
ax.set_xticks(range(5)); ax.set_xticklabels(stages, fontsize=7)
style(ax, "Disease stage", "CD200 (mean log1p)", "CD200 rises with progression"); panel(ax, "b")
fig.text(0.99, 0.01, "a: GSE161610/116256/117156, PMID 37919606 | b: GSE261917/161610/205704", ha="right", va="bottom", fontsize=6.5, color="gray")
plt.tight_layout()
save_all(fig, "Fig4_cross_cancer"); plt.close()

# ============ Fig 5: Survival + ibrutinib (4 panels) ============
m22762 = pd.read_csv(f"{OUT}/results/GSE22762_merged_151.tsv", sep="\t")
m39671 = pd.read_csv(f"{OUT}/results/GSE39671_merged.tsv", sep="\t")
_expr = pd.read_csv(f"{OUT}/data/GSE39671_expr.tsv", sep="\t")[["gsm","LPL"]]
m39671 = m39671.merge(_expr, on="gsm", how="left"); m39671["LPL_z"]=(m39671["LPL"]-m39671["LPL"].mean())/m39671["LPL"].std(ddof=0)
fig, axes = plt.subplots(1, 4, figsize=(7.09, 2.80))
def km(ax, m, gene, tcol, ecol, title):
    z = m[gene+"_z"].astype(float); hi = z>=z.median()
    t = m[tcol].astype(float); e = m[ecol].astype(int)
    lr = logrank_test(t[hi], t[~hi], e[hi], e[~hi])
    kf = KaplanMeierFitter()
    for grp,lbl,c,ls in [(True,"high","#1f6fb2","-"),(False,"low","#e8912a","--")]:
        kf.fit(t[hi==grp], e[hi==grp], label=f"{lbl} (n={(hi==grp).sum()})")
        kf.plot_survival_function(ax=ax, color=c, linestyle=ls, show_censors=True)
    ax.set_ylim(0,1.05); ax.set_title(f"{title}\nlog-rank p={lr.p_value:.3f}"); ax.legend(frameon=False, fontsize=6.5)
km(axes[0], m22762, "TNFSF9", "survival_time", "event", "GSE22762 TNFSF9 (OS)"); axes[0].set_xlabel("OS (days)"); axes[0].set_ylabel("survival"); panel(axes[0],"a")
km(axes[1], m39671, "TNFSF9", "ttft_years", "event", "GSE39671 TNFSF9 (TTFT)"); axes[1].set_xlabel("TTFT (years)"); panel(axes[1],"b")
km(axes[2], m39671, "LPL", "ttft_years", "event", "GSE39671 LPL (control)"); axes[2].set_xlabel("TTFT (years)"); panel(axes[2],"c")
# d: ibrutinib
ax = axes[3]
g11 = sc.read_h5ad(f"{OUT}/data/gse111014_processed.h5ad")
b11 = (g11.obs["celltype"]=="B").values
v = col(g11,"CD200"); tps=["d0","d30","d120","d150","d280"]
ys=[v[b11 & (g11.obs["timepoint"]==t).values].mean() for t in tps]
ax.plot(tps, ys, marker="o", color="#d84b6a", lw=2)
for t,y in zip(tps,ys): ax.text(t, y+0.012, f"{y:.2f}", ha="center", fontsize=6.5)
style(ax, "time", "CD200 (log1p)", "Ibrutinib (n=4)"); panel(ax,"d")
fig.text(0.99, 0.01, "a: GSE22762 | b,c: GSE39671 | d: GSE111014", ha="right", va="bottom", fontsize=6.5, color="gray")
plt.tight_layout()
save_all(fig, "Fig5_survival"); plt.close()

# ============ Fig 6: Druggability resource (1 panel) ============
fig, ax = plt.subplots(figsize=(3.46, 3.10))
dr = pd.read_csv(f"{OUT}/results/checkpoint_druggability.tsv", sep="\t")
dr = dr.sort_values("CLL", ascending=False)
# show top druggable checkpoints by CLL expression (raw values, not folds)
top = dr.head(8)
ax.barh(range(len(top))[::-1], top["CLL"], color="#3b7dd8", edgecolor="k", lw=1.0)
ax.set_yticks(range(len(top))[::-1]); ax.set_yticklabels([f"{r.ligand}→{r.receptor}" for _,r in top.iterrows()], fontsize=7.5)
for i,(_,r) in enumerate(top.iterrows()):
    ax.text(r["CLL"]+0.01, len(top)-1-i, f"{r.clinical_drug.split(' (')[0]}", va="center", fontsize=6.5, color="#555")
style(ax, "CLL expression (log1p)", None, "Druggable immune checkpoints (CLL LN)"); panel(ax, "a")
fig.text(0.99, 0.01, "GSE261917 (CLL vs reactive/FL)", ha="right", va="bottom", fontsize=6.5, color="gray")
plt.tight_layout()
save_all(fig, "Fig6_druggability"); plt.close()

print("saved 6 multi-panel figures (19 panels)")
