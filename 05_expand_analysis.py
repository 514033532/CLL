#!/usr/bin/env python3
"""
Step 12 — Expand analyses for a 5-main-figure paper.

1. GSE261917 full clustering + cell-type annotation (B/T/NK/myeloid).
2. GSE161610: cell-type proportions PB vs LN.
3. GSE161610: donor-level correlation B-cell CD200 vs effector exhaustion.
4. GSE161610: CD200+ vs CD200- B-cell subset features.
5. GSE111014: ibrutinib dynamics of CD200/CD70/TNFSF9/CD86/IL2RG.
"""
import scanpy as sc
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings("ignore")

OUT = "/home/zhh2/CLL-ScRNA-seqdatasets/NK_analysis"
sc.settings.verbosity = 0

# ================= 1. GSE261917 full annotation ============================
print("== GSE261917 full clustering + annotation ==")
g26 = sc.read_h5ad(f"{OUT}/data/gse261917_processed.h5ad")
sc.pp.highly_variable_genes(g26, n_top_genes=2000)
g26h = g26[:, g26.var["highly_variable"]].copy()
sc.pp.scale(g26h, max_value=10)
sc.tl.pca(g26h, n_comps=40, svd_solver="randomized", random_state=0)
sc.pp.neighbors(g26h, n_pcs=25, n_neighbors=15)
sc.tl.umap(g26h, random_state=0)
sc.tl.leiden(g26h, resolution=0.6, random_state=0)
g26.obs["leiden"] = g26h.obs["leiden"].astype(str)
g26.obsm["X_umap"] = g26h.obsm["X_umap"]

genes = list(g26.var_names)
def col(d, g): return np.asarray(d.layers["log1p"][:, list(d.var_names).index(g)].toarray()).ravel()
# annotate clusters by markers
markers = ["CD79A","MS4A1","CD19","CD5","CD3D","CD3E","CD3G","CD8A","CD4",
           "NKG7","GNLY","FCGR3A","KLRF1","NCAM1","CD14","LYZ","FCN1","PPBP"]
mk = {}
for cl in g26.obs["leiden"].unique():
    m = g26.obs["leiden"] == cl
    mk[cl] = [col(g26, g)[m].mean() for g in markers]
mkdf = pd.DataFrame(mk, index=markers).T
ct = {}
for cl in mkdf.index:
    r = mkdf.loc[cl]
    if r["NKG7"] > 0.8 and r["CD3G"] < 0.5 and (r["FCGR3A"] > 0.5 or r["KLRF1"] > 0.5):
        ct[cl] = "NK"
    elif r["PPBP"] > 1.0: ct[cl] = "Platelet"
    elif r["CD14"] > 0.5 or r["FCN1"] > 0.5 or r["LYZ"] > 1.0: ct[cl] = "Myeloid"
    elif r["CD79A"] > 0.5 or r["MS4A1"] > 0.5 or r["CD19"] > 0.5: ct[cl] = "B"
    elif r["CD3G"] > 0.5: ct[cl] = "CD8_T" if r["CD8A"] >= r["CD4"] else "CD4_T"
    else: ct[cl] = "Other"
g26.obs["celltype"] = g26.obs["leiden"].map(ct).astype("category")
print(g26.obs["celltype"].value_counts())
g26.write_h5ad(f"{OUT}/data/gse261917_annotated.h5ad")

# ================= 2. GSE161610 proportions PB vs LN ========================
print("\n== GSE161610 cell proportions ==")
g16 = sc.read_h5ad(f"{OUT}/data/gse161610_processed.h5ad")
prop = g16.obs.groupby(["tissue", "celltype"], observed=True).size().unstack().fillna(0)
prop = prop.div(prop.sum(1), axis=0)
prop.to_csv(f"{OUT}/results/gse161610_proportions.csv")
print(prop.round(3))

# ================= 3. B CD200 vs effector exhaustion (donor-level) ==========
print("\n== donor-level correlation: B CD200 vs effector GZMB ==")
def gcol(g): return np.asarray(g16.layers["log1p"][:, list(g16.var_names).index(g)].toarray()).ravel()
rows = []
for s in g16.obs["sample"].unique():
    m = (g16.obs["sample"] == s).values
    b = m & (g16.obs["celltype"] == "B").values
    cd8 = m & (g16.obs["celltype"] == "CD8_T").values
    nk = m & (g16.obs["celltype"] == "NK").values
    rows.append(dict(sample=s, tissue=g16.obs["tissue"][m].iloc[0],
                     B_CD200=gcol("CD200")[b].mean(), B_CD70=gcol("CD70")[b].mean(),
                     CD8_GZMB=gcol("GZMB")[cd8].mean(), CD8_TIGIT=gcol("TIGIT")[cd8].mean(),
                     NK_GZMB=gcol("GZMB")[nk].mean()))
cor = pd.DataFrame(rows)
from scipy.stats import spearmanr
cor.to_csv(f"{OUT}/results/gse161610_donor_correlation.csv", index=False)
# correlate across samples (only those with >=10 CD8/NK cells)
sub = cor[cor["CD8_GZMB"] > 0]
if len(sub) >= 4:
    r1, p1 = spearmanr(sub["B_CD200"], sub["CD8_GZMB"])
    r2, p2 = spearmanr(sub["B_CD200"], sub["CD8_TIGIT"])
    print(f"  B_CD200 ~ CD8_GZMB: rho={r1:.3f} p={p1:.3f} (n={len(sub)})")
    print(f"  B_CD200 ~ CD8_TIGIT: rho={r2:.3f} p={p2:.3f}")
print(cor.round(3).to_string(index=False))

# ================= 4. CD200+ vs CD200- B subset features =====================
print("\n== CD200+ vs CD200- B cell features (LN) ==")
bln = (g16.obs["celltype"] == "B").values & (g16.obs["tissue"] == "LN").values
cd200_pos = gcol("CD200")[bln] > 0
feat = ["MKI67","CD69","CD38","CD86","CXCR4","HLA-DRA","CD40","TNFSF9","CD70","IL2RG"]
for f in feat:
    v = gcol(f)[bln]
    print(f"  {f:8s} CD200+={v[cd200_pos].mean():.3f}  CD200-={v[~cd200_pos].mean():.3f}")

# ================= 5. ibrutinib dynamics (multiple genes) ===================
print("\n== GSE111014 ibrutinib dynamics ==")
g11 = sc.read_h5ad(f"{OUT}/data/gse111014_processed.h5ad")
b11 = (g11.obs["celltype"] == "B").values
tp = ["d0","d30","d120","d150","d280"]
for g in ["CD200","CD70","TNFSF9","CD86","IL2RG","CD40"]:
    v = np.asarray(g11.layers["log1p"][:, list(g11.var_names).index(g)].toarray()).ravel()
    row = []
    for t in tp:
        m = b11 & (g11.obs["timepoint"] == t).values
        if m.sum() > 50: row.append(f"{t}:{v[m].mean():.3f}")
    print(f"  {g:8s} " + "  ".join(row))
print("\nDONE")
