#!/usr/bin/env python3
"""Re-process GSE205704 (Richter) + GSE213516 (healthy) with full annotation."""
import scanpy as sc
import pandas as pd
import numpy as np
import scipy.io as sio
import gzip, glob, os, warnings
warnings.filterwarnings("ignore")

EXT = "/home/zhh2/CLL-ScRNA-seqdatasets/extracted"
OUT = "/home/zhh2/CLL-ScRNA-seqdatasets/NK_analysis/data"
sc.settings.verbosity = 0

def load_10x(dset, sid):
    d = f"{EXT}/{dset}"
    mtx = sio.mmread(gzip.open(glob.glob(f"{d}/{sid}_matrix.mtx.gz")[0], "rb")).tocsr().T
    fea = pd.read_csv(gzip.open(glob.glob(f"{d}/{sid}_features.tsv.gz")[0], "rt"), sep="\t", header=None)
    sym = fea.iloc[:, 1].astype(str)
    seen, uniq = {}, []
    for s in sym:
        if s in seen: seen[s]+=1; uniq.append(f"{s}.{seen[s]}")
        else: seen[s]=0; uniq.append(s)
    bc = pd.read_csv(gzip.open(glob.glob(f"{d}/{sid}_barcodes.tsv.gz")[0], "rt"), sep="\t", header=None)[0].astype(str)
    a = sc.AnnData(mtx, var=pd.DataFrame(index=uniq)); a.var["gene"]=sym.values
    a.obs["sample"]=sid; a.obs_names=[f"{sid}_{b}" for b in bc]
    return a

def process(dset, samples, label):
    ad = [load_10x(dset, s) for s in samples]
    adata = sc.concat(ad, join="outer", index_unique=None); adata.var_names_make_unique()
    adata.var["mito"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mito"], inplace=True, log1p=False)
    sc.pp.filter_cells(adata, min_genes=200)
    adata = adata[adata.obs["pct_counts_mito"] < 20].copy()
    sc.pp.normalize_total(adata, target_sum=1e4); sc.pp.log1p(adata)
    adata.layers["log1p"] = adata.X.copy()
    # 快速标记注释（不聚类，用标记基因）
    genes = list(adata.var_names)
    def col(g): return np.asarray(adata.layers["log1p"][:, genes.index(g)].toarray()).ravel() if g in genes else None
    cd3 = col("CD3D"); nkg7 = col("NKG7"); cd79a = col("CD79A"); ms4a1 = col("MS4A1")
    cd14 = col("CD14"); lyz = col("LYZ"); cd8a = col("CD8A"); cd4 = col("CD4"); mpo = col("MPO")
    ct = pd.Series("Other", index=adata.obs.index)
    ct[(cd79a > 0.5) | (ms4a1 > 0.5)] = "B"
    ct[(cd3 > 0.5) & (cd8a > 0.5) & (nkg7 < 0.5)] = "CD8_T"
    ct[(cd3 > 0.5) & (cd4 > 0.5)] = "CD4_T"
    ct[(nkg7 > 0.5) & (cd3 < 0.5)] = "NK"
    ct[(cd14 > 0.5) | (lyz > 1.0) | (mpo > 0.5)] = "Myeloid"
    adata.obs["celltype"] = ct.astype("category")
    adata.obs["group"] = label
    print(f"{dset} ({label}): {adata.n_obs} cells, {dict(ct.value_counts())}")
    adata.write_h5ad(f"{OUT}/{dset.lower()}_processed.h5ad")
    return adata

# GSE213516 healthy
hs = sorted(f[:-len("_matrix.mtx.gz")] for f in os.listdir(f"{EXT}/GSE213516_RAW") if f.endswith("_matrix.mtx.gz"))
process("GSE213516_RAW", hs, "healthy")
# GSE205704 Richter
rt = sorted(f[:-len("_matrix.mtx.gz")] for f in os.listdir(f"{EXT}/GSE205704_RAW") if f.endswith("_matrix.mtx.gz"))
process("GSE205704_RAW", rt, "Richter")
print("DONE")
