#!/usr/bin/env python3
"""
Composition analysis — full immune cell subset proportions across 6 cohorts.
Subsets: T(CD4 naive/EM/Treg, CD8 naive/EM/EMRA), NK(dim/bright/adaptive),
Mono(classical/nonclassical/intermediate), DC(cDC1/cDC2/pDC), B.
"""
import scanpy as sc
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings("ignore")

OUT = "/home/zhh2/CLL-ScRNA-seqdatasets/NK_analysis"
sc.settings.verbosity = 0

datasets = {
    "GSE161610": ("data/gse161610_processed.h5ad", "tissue"),
    "GSE111014": ("data/gse111014_processed.h5ad", "timepoint"),
    "GSE261917": ("data/gse261917_annotated.h5ad", "disease"),
    "GSE241842": ("data/gse241842_cml.h5ad", "group"),
    "GSE213516": ("data/gse213516_raw_processed.h5ad", "group"),
    "GSE205704": ("data/gse205704_raw_processed.h5ad", "group"),
}

def col(d, g):
    genes = list(d.var_names)
    return np.asarray(d.layers["log1p"][:, genes.index(g)].toarray()).ravel() if g in genes else None

rows = []
for name, (f, grpcol) in datasets.items():
    a = sc.read_h5ad(f"{OUT}/{f}")
    ct = a.obs["celltype"].astype(str)
    grp = a.obs[grpcol].astype(str)
    # ---- T subsets ----
    is_cd4 = ct.str.contains("CD4", regex=False)
    is_cd8 = ct.str.contains("CD8", regex=False)
    cd4 = col(a, "CD4"); cd8a = col(a, "CD8A"); foxp3 = col(a, "FOXP3")
    ccr7 = col(a, "CCR7"); il7r = col(a, "IL7R"); gzmb = col(a, "GZMB"); prf1 = col(a, "PRF1")
    treg = is_cd4 & (foxp3 > 0.5)
    cd4_naive = is_cd4 & (~treg) & (ccr7 > 0.5)
    cd4_em = is_cd4 & (~treg) & (ccr7 < 0.5)
    cd8_naive = is_cd8 & (ccr7 > 0.5) & (gzmb < 0.3)
    cd8_em = is_cd8 & (ccr7 < 0.5) & (gzmb > 0.3) & (prf1 > 0.3)
    cd8_emra = is_cd8 & (gzmb > 0.5) & (prf1 > 0.5)
    # ---- NK ----
    is_nk = ct.str.contains("NK", regex=False)
    fcgr3a = col(a, "FCGR3A"); ncam1 = col(a, "NCAM1"); klrg1 = col(a, "KLRG1")
    nk_dim = is_nk & (fcgr3a > 0.5) & (ncam1 < 0.3)
    nk_bright = is_nk & (ncam1 > 0.5) & (fcgr3a < 0.3)
    nk_adapt = is_nk & (klrg1 > 0.5)
    # ---- Mono ----
    is_myeloid = ct.str.contains("Myeloid", regex=False) | ct.str.contains("Mono", regex=False)
    cd14 = col(a, "CD14"); cd16m = col(a, "FCGR3A")
    mono_class = is_myeloid & (cd14 > 1.0) & (cd16m < 0.5)
    mono_nc = is_myeloid & (cd14 < 0.5) & (cd16m > 0.5)
    mono_int = is_myeloid & (cd14 > 0.5) & (cd16m > 0.5)
    # ---- DC ----
    cdc1 = is_myeloid & (col(a, "CLEC9A") > 0.5)
    cdc2 = is_myeloid & (col(a, "CD1C") > 0.5)
    pdc = is_myeloid & ((col(a, "CLEC4C") > 0.5) | (col(a, "LILRA4") > 0.5))
    # ---- B ----
    is_b = ct.str.contains("^B$", regex=True) | ct.str.contains("B cell", regex=False)

    subsets = {
        "CD4_Treg": treg, "CD4_naive": cd4_naive, "CD4_EM": cd4_em,
        "CD8_naive": cd8_naive, "CD8_EM": cd8_em, "CD8_EMRA": cd8_emra,
        "NK_CD56dim": nk_dim, "NK_CD56bright": nk_bright, "NK_adaptive": nk_adapt,
        "Mono_classical": mono_class, "Mono_nonclassical": mono_nc, "Mono_intermediate": mono_int,
        "cDC1": cdc1, "cDC2": cdc2, "pDC": pdc, "B": is_b,
    }
    for sub, mask in subsets.items():
        for g in grp.unique():
            rows.append((name, g, sub, int((mask.values & (grp == g).values).sum())))

res = pd.DataFrame(rows, columns=["dataset", "group", "subset", "n"])
# 比例 = 亚群数 / 该组总细胞数
total = res.groupby(["dataset", "group"])["n"].sum()
res["total"] = res.set_index(["dataset", "group"]).index.map(total)
res["frac"] = res["n"] / res["total"]
res.to_csv(f"{OUT}/results/immune_composition.csv", index=False)

# 汇总表
piv = res.pivot_table(index="subset", columns=["dataset", "group"], values="n", aggfunc="sum").fillna(0).astype(int)
print("=== 免疫亚群细胞数（跨数据集×分组）===")
print(piv.to_string())
print("\n已保存 immune_composition.csv")
