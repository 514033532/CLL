#!/usr/bin/env python3
"""
Step 13 — Five main figures (reviewer-formatted, A-D/A-B panels).

Fig1 atlas+QC: UMAP celltype / UMAP tissue / proportions / marker dotplot
Fig2 discovery: CD200 feature / CD200 violin / paired per-patient / CD70+TNFSF9
Fig3 receptor+function: CD27/CD200R1 expr / CD27+effector / correlation / LN-PB exhaust
Fig4 disease: GSE261917 UMAP / disease bars / progression / (CD200+ B features)
Fig5 therapy+model: ibrutinib multi-ligand / dual-axis schematic
"""
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import warnings
warnings.filterwarnings("ignore")

OUT = "/home/zhh2/CLL-ScRNA-seqdatasets/NK_analysis"
plt.rcParams.update({"font.size": 9, "axes.labelsize": 10, "axes.labelweight": "bold",
                     "axes.titleweight": "bold", "axes.titlesize": 10,
                     "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 7.5,
                     "axes.linewidth": 0.8, "figure.dpi": 200, "pdf.fonttype": 42})

CT_COLORS = {"B": "#3b7dd8", "CD4_T": "#57a05a", "CD8_T": "#e8912a",
             "NK": "#9b5de5", "Myeloid": "#d84b6a", "Other": "#9aa0a6"}
TISSUE_COLORS = {"PB": "#c9d4e3", "LN": "#1f3a5f"}

def add_panel(ax, L):
    ax.text(-0.16, 1.12, L, transform=ax.transAxes, fontsize=16, fontweight="bold", va="top", ha="left")
def style(ax, x=None, y=None, t=None):
    if x: ax.set_xlabel(x)
    if y: ax.set_ylabel(y)
    if t: ax.set_title(t)
    for s in ["top", "right"]: ax.spines[s].set_visible(False)

g16 = sc.read_h5ad(f"{OUT}/data/gse161610_processed.h5ad")
g26 = sc.read_h5ad(f"{OUT}/data/gse261917_annotated.h5ad")
g11 = sc.read_h5ad(f"{OUT}/data/gse111014_processed.h5ad")
def col(d, g): return np.asarray(d.layers["log1p"][:, list(d.var_names).index(g)].toarray()).ravel()

# ================= Figure 1: atlas + QC (2x2) ==============================
fig, axes = plt.subplots(2, 2, figsize=(9, 8))
ax = axes[0,0]
for c in sorted(g16.obs["celltype"].unique()):
    m = g16.obs["celltype"]==c
    ax.scatter(g16.obsm["X_umap"][m,0], g16.obsm["X_umap"][m,1], s=1.5, alpha=0.5,
               color=CT_COLORS.get(c), label=c, rasterized=True)
style(ax, "UMAP1", "UMAP2", "Cell types")
ax.legend(markerscale=5, frameon=False, loc="center left", bbox_to_anchor=(1.02,0.5)); add_panel(ax,"A")

ax = axes[0,1]
for t in ["PB","LN"]:
    m = g16.obs["tissue"]==t
    ax.scatter(g16.obsm["X_umap"][m,0], g16.obsm["X_umap"][m,1], s=1.5, alpha=0.5,
               color=TISSUE_COLORS[t], label=t, rasterized=True)
style(ax, "UMAP1", "UMAP2", "Tissue"); ax.legend(frameon=False, loc="lower right"); add_panel(ax,"B")

ax = axes[1,0]
prop = pd.read_csv(f"{OUT}/results/gse161610_proportions.csv", index_col=0)
prop = prop[["B","CD4_T","CD8_T","Myeloid","NK"]]
prop.plot(kind="bar", stacked=True, ax=ax, color=[CT_COLORS[c] for c in prop.columns],
          edgecolor="k", lw=0.5, legend=False)
ax.set_ylabel("fraction"); ax.set_title("Cell composition"); ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
ax.set_xlabel("Tissue")
ax.legend(frameon=False, loc="center left", bbox_to_anchor=(1.02, 0.5)); add_panel(ax,"C")

ax = axes[1,1]
# marker dotplot (manual)
cts = ["B","CD4_T","CD8_T","NK","Myeloid"]
mks = ["CD79A","MS4A1","CD5","CD3D","CD4","CD8A","NKG7","GNLY","FCGR3A","CD14","LYZ"]
mks = [g for g in mks if g in g16.var_names]
frac = np.zeros((len(cts), len(mks))); mean = np.zeros((len(cts), len(mks)))
for i,ct in enumerate(cts):
    m = (g16.obs["celltype"]==ct).values
    for j,g in enumerate(mks):
        v = col(g16,g)[m]; frac[i,j]=(v>0).mean(); mean[i,j]=v.mean()
mean_norm = mean/mean.max()
for i in range(len(cts)):
    for j in range(len(mks)):
        ax.scatter(j, i, s=frac[i,j]*400, color=plt.cm.viridis(mean_norm[i,j]), edgecolor="k", lw=0.3)
ax.set_xticks(range(len(mks))); ax.set_xticklabels(mks, rotation=45, ha="right", fontsize=7)
ax.set_yticks(range(len(cts))); ax.set_yticklabels(cts)
sm = plt.cm.ScalarMappable(cmap="viridis", norm=plt.Normalize(0, 1))
sm.set_array([])
plt.colorbar(sm, ax=ax, shrink=0.6, label="scaled mean expr")
style(ax, "Marker genes", "Cell type", "Marker genes"); add_panel(ax,"D")
plt.tight_layout()
fig.text(0.99, 0.01, "GSE161610", ha="right", va="bottom", fontsize=7, color="gray")
plt.savefig(f"{OUT}/figures/F1_atlas.png", bbox_inches="tight")
plt.savefig(f"{OUT}/figures/F1_atlas.pdf", bbox_inches="tight"); plt.close()

# ================= Figure 2: discovery (2x2) ===============================
fig, axes = plt.subplots(2, 2, figsize=(9, 8))
ax = axes[0,0]
cd200 = col(g16, "CD200")
idx = np.argsort(cd200)
sc = ax.scatter(g16.obsm["X_umap"][idx,0], g16.obsm["X_umap"][idx,1], s=1.5, c=cd200[idx],
                cmap="Reds", alpha=0.7, rasterized=True)
plt.colorbar(sc, ax=ax, shrink=0.8, label="CD200")
style(ax, "UMAP1", "UMAP2", "CD200 feature"); add_panel(ax,"A")

ax = axes[0,1]
b16 = (g16.obs["celltype"]=="B").values
data=[cd200[b16 & (g16.obs["tissue"]=="PB").values], cd200[b16 & (g16.obs["tissue"]=="LN").values]]
vp=ax.violinplot(data, positions=[0,1], showmeans=True, showextrema=False)
for pc,c in zip(vp["bodies"], [TISSUE_COLORS["PB"],TISSUE_COLORS["LN"]]): pc.set_facecolor(c); pc.set_alpha(0.75); pc.set_edgecolor("k")
ax.set_xticks([0,1]); ax.set_xticklabels(["PB","LN"], fontweight="bold")
style(ax, "Tissue", "CD200 (log1p)", "CD200 PB vs LN (B cells)")
ax.legend([mpatches.Patch(color=TISSUE_COLORS["PB"]),mpatches.Patch(color=TISSUE_COLORS["LN"])],["PB","LN"],frameon=False,loc="upper left"); add_panel(ax,"B")

ax = axes[1,0]
# per-patient paired CD200
for p, colr in zip(sorted(g16.obs["patient"].unique()), plt.cm.tab10.colors):
    pb = cd200[b16 & (g16.obs["patient"]==p).values & (g16.obs["tissue"]=="PB").values].mean()
    ln = cd200[b16 & (g16.obs["patient"]==p).values & (g16.obs["tissue"]=="LN").values].mean()
    ax.plot([0,1],[pb,ln], marker="o", markersize=5, color=colr, label=p)
ax.set_xticks([0,1]); ax.set_xticklabels(["PB","LN"], fontweight="bold")
style(ax, "Tissue", "CD200 (mean)", "Paired PB→LN (5 patients)")
ax.legend(frameon=False, loc="center left", bbox_to_anchor=(1.02,0.5)); add_panel(ax,"C")

ax = axes[1,1]
for gi, (g, color) in enumerate([("CD70","#d84b6a"),("TNFSF9","#6cb2eb")]):
    off = gi * 0.35
    d=[col(g16,g)[b16 & (g16.obs["tissue"]=="PB").values], col(g16,g)[b16 & (g16.obs["tissue"]=="LN").values]]
    vp=ax.violinplot(d, positions=[0+off,1+off], showmeans=True, showextrema=False, widths=0.3)
    for pc in vp["bodies"]: pc.set_facecolor(color); pc.set_alpha(0.7); pc.set_edgecolor("k")
ax.set_xticks([0.17, 1.17]); ax.set_xticklabels(["PB","LN"], fontweight="bold")
style(ax, "Tissue", "log1p expr", "CD70 / TNFSF9 (B cells)")
ax.legend([mpatches.Patch(color="#d84b6a"), mpatches.Patch(color="#6cb2eb")], ["CD70","TNFSF9"], frameon=False, loc="upper left")
add_panel(ax,"D")
plt.tight_layout()
fig.text(0.99, 0.01, "GSE161610", ha="right", va="bottom", fontsize=7, color="gray")
plt.savefig(f"{OUT}/figures/F2_discovery.png", bbox_inches="tight")
plt.savefig(f"{OUT}/figures/F2_discovery.pdf", bbox_inches="tight"); plt.close()

# ================= Figure 3: receptor + function (2x2) =====================
fig, axes = plt.subplots(2, 2, figsize=(9, 8))
ax = axes[0,0]
# CD27 / CD200R1 expression across cell types
cts = ["B","CD4_T","CD8_T","NK","Myeloid"]
cd27 = [col(g16,"CD27")[(g16.obs["celltype"]==c).values].mean() for c in cts]
cd200r1 = [col(g16,"CD200R1")[(g16.obs["celltype"]==c).values].mean() for c in cts]
x=np.arange(5); w=0.38
ax.bar(x-w/2, cd27, w, label="CD27", color="#1f6fb2", edgecolor="k", lw=0.5)
ax.bar(x+w/2, cd200r1, w, label="CD200R1", color="#c0392b", edgecolor="k", lw=0.5)
ax.set_xticks(x); ax.set_xticklabels(cts)
style(ax, "Cell type", "mean log1p expr", "Receptors (CD27 / CD200R1)")
ax.legend(frameon=False, loc="center left", bbox_to_anchor=(1.02, 0.5)); add_panel(ax,"A")

ax = axes[0,1]
m=(g16.obs["celltype"]=="CD8_T").values
cd27pos=col(g16,"CD27")[m]>0
gzmb=col(g16,"GZMB")[m]; gnly=col(g16,"GNLY")[m]; prf1=col(g16,"PRF1")[m]
x=np.arange(2); w=0.25
for i,(nm,arr) in enumerate([("GZMB",gzmb),("GNLY",gnly),("PRF1",prf1)]):
    ax.bar(x-0.25+i*w, [arr[cd27pos].mean(), arr[~cd27pos].mean()], w, label=nm,
           color=plt.cm.tab10(i), edgecolor="k", lw=0.5)
ax.set_xticks(x); ax.set_xticklabels(["CD27+","CD27−"], fontweight="bold")
style(ax, "CD27 status", "mean log1p expr", "CD8 T cytotoxicity (CD27±)")
ax.legend(frameon=False, loc="upper left"); add_panel(ax,"B")

ax = axes[1,0]
cor = pd.read_csv(f"{OUT}/results/gse161610_donor_correlation.csv")
for t,cc in [("PB",TISSUE_COLORS["PB"]),("LN",TISSUE_COLORS["LN"])]:
    d=cor[cor["tissue"]==t]
    ax.scatter(d["B_CD200"], d["CD8_GZMB"], s=50, color=cc, edgecolor="k", label=t)
ax.set_xlabel("B-cell CD200 (mean)"); ax.set_ylabel("CD8 T GZMB (mean)")
ax.set_title("Donor-level correlation")
ax.text(0.05, 0.05, "ρ=−0.31, p=0.46 (n=8)", transform=ax.transAxes, fontsize=8)
ax.legend(frameon=False); add_panel(ax,"C")

ax = axes[1,1]
for ct,cc in [("CD8_T","#e8912a"),("NK","#9b5de5")]:
    m=(g16.obs["celltype"]==ct).values
    vals=[]
    for t in ["PB","LN"]:
        mm=m & (g16.obs["tissue"]==t).values
        vals.append(col(g16,"GZMB")[mm].mean())
    ax.plot([0,1],vals, marker="o", color=cc, label=ct)
ax.set_xticks([0,1]); ax.set_xticklabels(["PB","LN"], fontweight="bold")
style(ax, "Tissue", "GZMB (mean log1p)", "Effector cytotoxicity PB→LN")
ax.legend(frameon=False); add_panel(ax,"D")
plt.tight_layout()
fig.text(0.99, 0.01, "GSE161610", ha="right", va="bottom", fontsize=7, color="gray")
plt.savefig(f"{OUT}/figures/F3_function.png", bbox_inches="tight")
plt.savefig(f"{OUT}/figures/F3_function.pdf", bbox_inches="tight"); plt.close()

# ================= Figure 4: disease specificity (2x2) =====================
fig, axes = plt.subplots(2, 2, figsize=(9, 8))
ax = axes[0,0]
for c in sorted(g26.obs["celltype"].unique()):
    m=g26.obs["celltype"]==c
    ax.scatter(g26.obsm["X_umap"][m,0], g26.obsm["X_umap"][m,1], s=1.5, alpha=0.5,
               color=CT_COLORS.get(c,"#9aa0a6"), label=c, rasterized=True)
style(ax, "UMAP1","UMAP2","GSE261917 atlas")
ax.legend(markerscale=5, frameon=False, loc="center left", bbox_to_anchor=(1.02,0.5), fontsize=6.5); add_panel(ax,"A")

ax = axes[0,1]
g26g=list(g26.var_names)
def col26(g): return np.asarray(g26.layers["log1p"][:,g26g.index(g)].toarray()).ravel()
isb26=((col26("CD79A")>0.5)|(col26("MS4A1")>0.5))&(col26("CD3D")<0.5)
DIS=["Reactive","FL","CLL"]
cd200v=[col26("CD200")[isb26 & (g26.obs["disease"]==d).values].mean() for d in DIS]
cd70v=[col26("CD70")[isb26 & (g26.obs["disease"]==d).values].mean() for d in DIS]
x=np.arange(3); w=0.36
ax.bar(x-w/2, cd200v, w, label="CD200", color="#d84b6a", edgecolor="k", lw=0.5)
ax.bar(x+w/2, cd70v, w, label="CD70", color="#6cb2eb", edgecolor="k", lw=0.5)
ax.set_xticks(x); ax.set_xticklabels(DIS, fontweight="bold")
style(ax, "Disease", "mean log1p expr", "CD200/CD70 by disease (LN)")
ax.legend(frameon=False, loc="upper left"); add_panel(ax,"B")

ax = axes[1,0]
stages=["Reactive\nLN","FL\nLN","CLL\nPB","CLL\nLN","Richter"]
vals=[0.065,0.027,0.132,0.341,0.367]
cols=["#9aa0a6","#17a2b8","#6cb2eb","#d84b6a","#7b2e3e"]
ax.bar(range(5), vals, color=cols, alpha=0.85, edgecolor="k", lw=0.5)
for i,v in enumerate(vals): ax.text(i,v+0.012,f"{v:.2f}",ha="center",fontsize=8)
ax.set_xticks(range(5)); ax.set_xticklabels(stages, fontsize=7)
style(ax, "Disease stage", "CD200 (mean log1p)", "CD200 across progression"); add_panel(ax,"C")

ax = axes[1,1]
# CD200+ vs CD200- B cell features (LN)
bln=(g16.obs["celltype"]=="B").values & (g16.obs["tissue"]=="LN").values
pos=col(g16,"CD200")[bln]>0
feats=["MKI67","CD69","CD86","CXCR4","CD40","IL2RG"]
vals_pos=[col(g16,f)[bln][pos].mean() for f in feats]
vals_neg=[col(g16,f)[bln][~pos].mean() for f in feats]
x=np.arange(len(feats)); w=0.38
ax.bar(x-w/2, vals_pos, w, label="CD200+", color="#d84b6a", edgecolor="k", lw=0.5)
ax.bar(x+w/2, vals_neg, w, label="CD200−", color="#9aa0a6", edgecolor="k", lw=0.5)
ax.set_xticks(x); ax.set_xticklabels(feats, fontsize=8)
style(ax, "Gene", "mean log1p expr", "CD200+ vs CD200− B cells (LN)")
ax.legend(frameon=False, loc="upper left"); add_panel(ax,"D")
plt.tight_layout()
fig.text(0.99, 0.01, "A,B: GSE261917 | C: GSE261917+161610+205704 | D: GSE161610", ha="right", va="bottom", fontsize=7, color="gray")
plt.savefig(f"{OUT}/figures/F4_disease.png", bbox_inches="tight")
plt.savefig(f"{OUT}/figures/F4_disease.pdf", bbox_inches="tight"); plt.close()

# ================= Figure 5: therapy + model (1x2) ========================
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
ax = axes[0]
b11=(g11.obs["celltype"]=="B").values
tp=["d0","d30","d120","d150","d280"]
for g,cc in [("CD200","#d84b6a"),("CD70","#6cb2eb"),("CD86","#e8912a"),("CD40","#57a05a")]:
    v=col(g11,g); ys=[v[b11 & (g11.obs["timepoint"]==t).values].mean() for t in tp]
    ax.plot(range(5), ys, marker="o", markersize=4, color=cc, label=g)
ax.set_xticks(range(5)); ax.set_xticklabels(tp)
style(ax, "time on ibrutinib", "mean log1p expr", "Ligands decline on ibrutinib")
ax.legend(frameon=False, loc="center left", bbox_to_anchor=(1.02, 0.5), title="Ligand"); add_panel(ax,"A")

ax = axes[1]; ax.axis("off"); ax.set_xlim(0,10); ax.set_ylim(0,6)
def cell(x,y,r,label,color):
    ax.add_patch(plt.Circle((x,y),r,facecolor=color,edgecolor="k",lw=2,alpha=0.85,zorder=3))
    ax.text(x,y,label,ha="center",va="center",fontsize=8,fontweight="bold",zorder=4)
cell(5,3,1.2,"CLL B\ncell","#f4a7b9"); ax.text(5,4.4,"CD200  CD70",ha="center",fontsize=8,color="#a11",fontweight="bold")
cell(1.2,3,0.9,"Myeloid\ncell","#f9d77e"); ax.text(1.2,4.0,"CD200R1",ha="center",fontsize=8,color="#a11",fontweight="bold")
cell(8.8,3,0.9,"T/NK\ncell","#a9d6a5"); ax.text(8.8,4.0,"CD27",ha="center",fontsize=8,color="#a11",fontweight="bold")
ax.annotate("",xy=(2.1,3),xytext=(3.8,3.2),arrowprops=dict(arrowstyle="-|>",color="#c0392b",lw=2.5))
ax.text(2.9,3.85,"CD200→CD200R1\nsuppress myeloid",ha="center",fontsize=7.5,color="#c0392b",fontweight="bold")
ax.annotate("",xy=(7.9,3),xytext=(6.2,3.2),arrowprops=dict(arrowstyle="-|>",color="#1f6fb2",lw=2.5))
ax.text(7.1,3.85,"CD70→CD27\nT/NK exhaustion",ha="center",fontsize=7.5,color="#1f6fb2",fontweight="bold")
ax.text(5,5.6,"Model: dual-axis immune evasion",ha="center",fontsize=11,fontweight="bold")
add_panel(ax,"B")
plt.tight_layout()
fig.text(0.99, 0.01, "GSE111014", ha="right", va="bottom", fontsize=7, color="gray")
plt.savefig(f"{OUT}/figures/F5_model.png", bbox_inches="tight")
plt.savefig(f"{OUT}/figures/F5_model.pdf", bbox_inches="tight"); plt.close()

print("saved F1-F5 PNG")
print("DONE")
