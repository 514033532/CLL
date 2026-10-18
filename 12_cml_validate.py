#!/usr/bin/env python3
"""
Step 12 — Independent human CML scRNA validation (GSE241842, bone marrow).

Full QC/clustering/annotation, then check immune-evasion ligand expression on
CML tumor (myeloid) cells vs normal BM — to test the cross-cancer finding:
CLL = CD200/CD70/TNFSF9, CML = PVR (TIGIT ligand).
"""
import scanpy as sc
import pandas as pd
import numpy as np
import scipy.io as sio
import gzip, glob, os, warnings
warnings.filterwarnings("ignore")

EXT = "/home/zhh2/CML-ScRNA-seqdatasets/extracted/GSE241842_RAW"
OUT = "/home/zhh2/CLL-ScRNA-seqdatasets/NK_analysis"
sc.settings.verbosity = 0

# CML + Normal bone marrow samples (skip CD3-sorted)
samples = {
    "GSM7744314_CML1_BM": "CML", "GSM7744315_CML2_BM": "CML",
    "GSM7744316_CML3_BM": "CML", "GSM7744317_CML5_BM": "CML",
    "GSM7744318_CML6_BM": "CML",
    "GSM7744319_Normal2_BM": "Normal", "GSM7744320_Normal3_BM": "Normal",
}

def load(sid):
    mtx = sio.mmread(gzip.open(f"{EXT}/{sid}_matrix.mtx.gz", "rb")).tocsr().T
    fea = pd.read_csv(gzip.open(f"{EXT}/{sid}_features.tsv.gz", "rt"), sep="\t", header=None)
    sym = fea.iloc[:, 1].astype(str)
    seen, uniq = {}, []
    for s in sym:
        if s in seen: seen[s]+=1; uniq.append(f"{s}.{seen[s]}")
        else: seen[s]=0; uniq.append(s)
    bc = pd.read_csv(gzip.open(f"{EXT}/{sid}_barcodes.tsv.gz", "rt"), sep="\t", header=None)[0].astype(str)
    a = sc.AnnData(mtx, var=pd.DataFrame(index=uniq))
    a.var["gene"] = sym.values
    a.obs["sample"] = sid
    a.obs_names = [f"{sid}_{b}" for b in bc]
    return a

adatas = []
for sid, grp in samples.items():
    a = load(sid); a.obs["group"] = grp
    adatas.append(a)
adata = sc.concat(adatas, join="outer", index_unique=None)
adata.var_names_make_unique()
print(f"loaded {adata.n_obs} cells, {adata.n_vars} genes")

# QC
adata.var["mito"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mito"], inplace=True, log1p=False)
sc.pp.filter_cells(adata, min_genes=200)
adata = adata[adata.obs["pct_counts_mito"] < 20].copy()
sc.pp.filter_genes(adata, min_cells=3)
print(f"after QC: {adata.n_obs} cells")

# normalize + cluster
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers["log1p"] = adata.X.copy()
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
ah = adata[:, adata.var["highly_variable"]].copy()
sc.pp.scale(ah, max_value=10)
sc.tl.pca(ah, n_comps=40, svd_solver="randomized", random_state=0)
sc.pp.neighbors(ah, n_pcs=25, n_neighbors=15)
sc.tl.umap(ah, random_state=0)
sc.tl.leiden(ah, resolution=0.6, random_state=0)
adata.obs["leiden"] = ah.obs["leiden"].astype(str)
adata.obsm["X_umap"] = ah.obsm["X_umap"]

# annotate
genes = list(adata.var_names)
def col(g): return np.asarray(adata.layers["log1p"][:, genes.index(g)].toarray()).ravel() if g in genes else None
markers = ["MPO","ELANE","LYZ","CD33","CD14","FCGR3A","MS4A1","CD79A","CD19",
           "CD3D","CD3E","NKG7","GNLY","KLRD1","HBB","HBA1","MZB1","JCHAIN","PPBP"]
mk = {}
for cl in adata.obs["leiden"].unique():
    m = adata.obs["leiden"]==cl
    mk[cl] = [col(g)[m].mean() for g in markers]
mkdf = pd.DataFrame(mk, index=markers).T
ct = {}
for cl in mkdf.index:
    r = mkdf.loc[cl]
    if r["HBB"] > 0.5 or r["HBA1"] > 0.5: ct[cl] = "Erythroid"
    elif r["MPO"] > 0.5 or r["ELANE"] > 0.5: ct[cl] = "Myeloid/CML"
    elif r["CD14"] > 1.0: ct[cl] = "Monocyte"
    elif r["CD79A"] > 0.5 or r["MS4A1"] > 0.5: ct[cl] = "B"
    elif r["NKG7"] > 0.5 and r["CD3D"] < 0.5: ct[cl] = "NK"
    elif r["CD3D"] > 0.5: ct[cl] = "T"
    elif r["MZB1"] > 0.5: ct[cl] = "Plasma"
    elif r["PPBP"] > 1.0: ct[cl] = "Platelet"
    else: ct[cl] = "Other"
adata.obs["celltype"] = adata.obs["leiden"].map(ct).astype("category")
print("\ncelltype counts:")
print(adata.obs["celltype"].value_counts())

# ligand expression on CML tumor (Myeloid/CML) vs normal myeloid
LIG = ["CD200","CD70","TNFSF9","CD86","LGALS9","PVR"]
print("\n=== 配体表达：CML 肿瘤(髓系) vs 正常骨髓髓系 ===")
cml_myeloid = (adata.obs["celltype"]=="Myeloid/CML").values & (adata.obs["group"]=="CML").values
norm_myeloid = (adata.obs["celltype"]=="Myeloid/CML").values & (adata.obs["group"]=="Normal").values
rows=[]
for g in LIG:
    if col(g) is None: continue
    c = col(g)[cml_myeloid].mean() if cml_myeloid.sum()>0 else np.nan
    n = col(g)[norm_myeloid].mean() if norm_myeloid.sum()>0 else np.nan
    rows.append((g, c, n))
df = pd.DataFrame(rows, columns=["gene","CML_myeloid","Normal_myeloid"])
print(f"CML肿瘤细胞 n={cml_myeloid.sum()}, 正常髓系 n={norm_myeloid.sum()}")
print(df.round(3).to_string(index=False))
df.to_csv(f"{OUT}/results/gse241842_cml_ligand.tsv", sep="\t", index=False)

adata.X = adata.layers["log1p"].copy()
adata.write_h5ad(f"{OUT}/data/gse241842_cml.h5ad")
print("\nsaved gse241842_cml.h5ad + gse241842_cml_ligand.tsv")
