#!/usr/bin/env python3
"""6 figures (>=20 panels): innate immune compartment depletion in CLL."""
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "/home/zhh2/CLL-ScRNA-seqdatasets/NK_analysis"
plt.rcParams.update({"font.size": 8.5, "axes.labelsize": 9, "axes.labelweight": "bold",
                     "axes.titleweight": "bold", "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
                     "legend.fontsize": 7, "figure.dpi": 200, "pdf.fonttype": 42})
def add_panel(ax, L): ax.text(-0.16, 1.12, L, transform=ax.transAxes, fontsize=15, fontweight="bold", va="top", ha="left")
def style(ax, x=None, y=None, t=None):
    if x: ax.set_xlabel(x)
    if y: ax.set_ylabel(y)
    if t: ax.set_title(t)
    for s in ["top", "right"]: ax.spines[s].set_visible(False)

res = pd.read_csv(f"{OUT}/results/immune_composition.csv")
innate = ["NK_CD56dim","NK_CD56bright","NK_adaptive","Mono_classical","Mono_nonclassical","Mono_intermediate","cDC1","cDC2","pDC"]
CT_COLORS = {"B":"#3b7dd8","CD4_T":"#57a05a","CD8_T":"#e8912a","NK":"#9b5de5","Myeloid":"#d84b6a","pDC":"#17a2b8","Other":"#9aa0a6"}
spectrum = [("GSE213516","healthy","Healthy"), ("GSE161610","PB","CLL PB"),
            ("GSE161610","LN","CLL LN"), ("GSE205704","Richter","Richter")]

def frac(ds, g, ss):
    sub = res[(res.dataset==ds) & (res.group==g)]
    if len(sub)==0: return 0.0
    total = sub["total"].iloc[0]
    return sub[sub.subset.isin(ss)]["n"].sum()/total*100

def bar_spectrum(ax, ss, color, title):
    vals = [frac(ds, g, ss) for ds, g, _ in spectrum]
    ax.bar(range(4), vals, color=color, edgecolor="k", lw=0.5)
    for i, v in enumerate(vals): ax.text(i, v+max(vals)*0.03, f"{v:.1f}", ha="center", fontsize=7)
    ax.set_xticks(range(4)); ax.set_xticklabels([s[2] for s in spectrum], fontsize=7)
    ax.set_ylabel("% of cells"); ax.set_title(title)

# ============ Fig 1: atlas (4 panels) ============
def ct_composition(ds, g):
    sub = res[(res.dataset==ds) & (res.group==g)]
    if len(sub)==0: return []
    total = sub["total"].iloc[0]
    cts = ["B","CD4_T","CD8_T","NK","Myeloid","pDC"]
    return [sub[sub.subset.isin({"B":["B"],"CD4_T":["CD4_Treg","CD4_naive","CD4_EM"],"CD8_T":["CD8_naive","CD8_EM","CD8_EMRA"],"NK":["NK_CD56dim","NK_CD56bright","NK_adaptive"],"Myeloid":["Mono_classical","Mono_intermediate","Mono_nonclassical"],"pDC":["pDC"]}[ct])]["n"].sum()/total*100 if total else 0 for ct in cts]

fig, axes = plt.subplots(1, 4, figsize=(16, 4))
overview = [("Healthy (GSE213516)","GSE213516","healthy"),
            ("CLL (GSE161610)","GSE161610","PB"),
            ("CML (GSE241842)","GSE241842","CML")]
cts = ["B","CD4_T","CD8_T","NK","Myeloid","pDC"]
ctcols = ["#3b7dd8","#57a05a","#e8912a","#9b5de5","#d84b6a","#17a2b8"]
for i, (name, ds, g) in enumerate(overview):
    ax = axes[i]
    vals = ct_composition(ds, g)
    ax.bar(range(len(cts)), vals, color=ctcols, edgecolor="k", lw=0.5)
    ax.set_xticks(range(len(cts))); ax.set_xticklabels(cts, fontsize=7, rotation=30, ha="right")
    ax.set_ylabel("% of cells"); ax.set_title(name); add_panel(ax, "ABC"[i])
ax = axes[3]
bar_spectrum(ax, innate, "#1f3a5f", "Innate compartment")
add_panel(ax, "D")
plt.tight_layout(); plt.savefig(f"{OUT}/figures/IF1_overview.png", bbox_inches="tight")
plt.savefig(f"{OUT}/figures/IF1_overview.pdf", bbox_inches="tight"); plt.close()

# ============ Fig 2-4: subsets (each 4 panels) ============
def subset_fig(subsets, colors, titles, fname, Ls="ABCD"):
    fig, axes = plt.subplots(1, 4, figsize=(16, 4))
    for i, (ss, color, title) in enumerate(zip(subsets, colors, titles)):
        bar_spectrum(axes[i], [ss], color, title); add_panel(axes[i], Ls[i])
    plt.tight_layout(); plt.savefig(f"{OUT}/figures/{fname}.png", bbox_inches="tight")
    plt.savefig(f"{OUT}/figures/{fname}.pdf", bbox_inches="tight"); plt.close()

subset_fig([["NK_CD56dim"],["NK_CD56bright"],["NK_adaptive"],["NK_CD56dim","NK_CD56bright","NK_adaptive"]],
           ["#4b0082","#9b5de5","#c9a0e5","#6a3fa0"],
           ["NK CD56dim","NK CD56bright","NK adaptive","Total NK"], "IF2_NK")
subset_fig([["Mono_classical"],["Mono_intermediate"],["Mono_nonclassical"],["Mono_classical","Mono_intermediate","Mono_nonclassical"]],
           ["#d84b6a","#f4a7b9","#8b1a3a","#b03060"],
           ["Mono classical","Mono intermediate","Mono non-classical","Total monocytes"], "IF3_mono")
subset_fig([["cDC1"],["cDC2"],["pDC"],["cDC1","cDC2","pDC"]],
           ["#17a2b8","#6cc5d0","#0d6a75","#0f8a9a"],
           ["cDC1","cDC2","pDC","Total DC"], "IF4_DC")

# ============ Fig 5: tissue & therapy (4 panels) ============
fig, axes = plt.subplots(1, 4, figsize=(16, 4))
# A: PB vs LN
ax = axes[0]
for i, (ds, g, lab) in enumerate([("GSE161610","PB","PB"),("GSE161610","LN","LN")]):
    ax.bar(i, frac(ds, g, innate), color="#1f3a5f", edgecolor="k", lw=0.5)
    ax.text(i, frac(ds,g,innate)+0.05, f"{frac(ds,g,innate):.1f}", ha="center", fontsize=8)
ax.set_xticks(range(2)); ax.set_xticklabels(["PB","LN"]); style(ax, None, "% innate", "CLL PB vs LN"); add_panel(ax,"A")
# B: disease (Reactive/FL/CLL LN)
ax = axes[1]
for i, (ds, g, lab) in enumerate([("GSE261917","Reactive","Reactive"),("GSE261917","FL","FL"),("GSE261917","CLL","CLL")]):
    ax.bar(i, frac(ds, g, innate), color="#17a2b8", edgecolor="k", lw=0.5)
    ax.text(i, frac(ds,g,innate)+0.02, f"{frac(ds,g,innate):.1f}", ha="center", fontsize=8)
ax.set_xticks(range(3)); ax.set_xticklabels(["Reactive","FL","CLL"]); style(ax, None, "% innate", "LN disease"); add_panel(ax,"B")
# C: Richter restoration
ax = axes[2]
for i, (ds, g, lab) in enumerate([("GSE161610","PB","CLL PB"),("GSE161610","LN","CLL LN"),("GSE205704","Richter","Richter")]):
    ax.bar(i, frac(ds, g, innate), color=["#6cb2eb","#1f3a5f","#7b2e3e"][i], edgecolor="k", lw=0.5)
    ax.text(i, frac(ds,g,innate)+0.1, f"{frac(ds,g,innate):.1f}", ha="center", fontsize=8)
ax.set_xticks(range(3)); ax.set_xticklabels(["CLL PB","CLL LN","Richter"]); style(ax, None, "% innate", "Richter restoration"); add_panel(ax,"C")
# D: ibrutinib
ax = axes[3]
sub = res[res.dataset=="GSE111014"]
tps = ["d0","d30","d120","d150","d280"]
for tp in tps:
    g = sub[sub.group==tp]; total = g["total"].iloc[0] if len(g) else 0
    ax.bar(tps.index(tp), g[g.subset.isin(innate)]["n"].sum()/total*100 if total else 0, color="#d84b6a", edgecolor="k", lw=0.5)
ax.set_xticks(range(5)); ax.set_xticklabels(tps); style(ax, "time on ibrutinib", "% innate", "BTKi (GSE111014)"); add_panel(ax,"D")
plt.tight_layout(); plt.savefig(f"{OUT}/figures/IF5_tissue_therapy.png", bbox_inches="tight")
plt.savefig(f"{OUT}/figures/IF5_tissue_therapy.pdf", bbox_inches="tight"); plt.close()

# ============ Fig 6: cross-cancer (4 panels) ============
fig, axes = plt.subplots(1, 4, figsize=(16, 4))
cats = [("NK",["NK_CD56dim","NK_CD56bright","NK_adaptive"]),
        ("Mono",["Mono_classical","Mono_nonclassical","Mono_intermediate"]),
        ("DC",["cDC1","cDC2","pDC"])]
for i, (lab, ds, g) in enumerate([("Healthy","GSE213516","healthy"),("CLL PB","GSE161610","PB"),
                                  ("Richter","GSE205704","Richter"),("CML","GSE241842","CML")]):
    ax = axes[i]
    vals = [frac(ds, g, ss) for _, ss in cats]
    ax.bar(range(3), vals, color=["#9b5de5","#d84b6a","#17a2b8"], edgecolor="k", lw=0.5)
    for j, v in enumerate(vals): ax.text(j, v+max(vals)*0.03, f"{v:.1f}", ha="center", fontsize=7.5)
    ax.set_xticks(range(3)); ax.set_xticklabels(["NK","Mono","DC"])
    style(ax, None, "% of cells", lab); add_panel(ax, "ABCD"[i])
plt.tight_layout(); plt.savefig(f"{OUT}/figures/IF6_crosscancer.png", bbox_inches="tight")
plt.savefig(f"{OUT}/figures/IF6_crosscancer.pdf", bbox_inches="tight"); plt.close()
print("saved IF1-IF6 (23 panels)")
