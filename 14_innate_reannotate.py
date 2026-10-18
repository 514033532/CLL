#!/usr/bin/env python3
"""
New pipeline — fine-grained re-annotation of the innate immune compartment
(NK subsets, monocyte subsets, dendritic cell subsets) across 6 CLL cohorts.

NK: CD56dim (FCGR3A+ NCAM1-), CD56bright (NCAM1+ FCGR3A-), adaptive (KLRG1+)
Mono: classical (CD14+CD16-/FCGR3A-), non-classical (CD14-CD16+), intermediate
DC: cDC1 (CLEC9A+), cDC2 (CD1C+), pDC (CLEC4C+/LILRA4+)
"""
import scanpy as sc
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings("ignore")

OUT = "/home/zhh2/CLL-ScRNA-seqdatasets/NK_analysis"
sc.settings.verbosity = 0

datasets = {
    "GSE161610": "data/gse161610_processed.h5ad",
    "GSE111014": "data/gse111014_processed.h5ad",
    "GSE261917": "data/gse261917_annotated.h5ad",
}

def col(d, g):
    genes = list(d.var_names)
    return np.asarray(d.layers["log1p"][:, genes.index(g)].toarray()).ravel() if g in genes else None

rows = []
for name, f in datasets.items():
    a = sc.read_h5ad(f"{OUT}/{f}")
    genes = list(a.var_names)
    ct = a.obs["celltype"].astype(str)
    # NK subsets
    is_nk = ct.str.contains("NK", regex=False)
    nkg7 = col(a, "NKG7"); fcgr3a = col(a, "FCGR3A"); ncam1 = col(a, "NCAM1"); klrg1 = col(a, "KLRG1")
    if is_nk.sum() > 0:
        nk_dim = is_nk & (fcgr3a > 0.5) & (ncam1 < 0.3)
        nk_bright = is_nk & (ncam1 > 0.5) & (fcgr3a < 0.3)
        nk_adapt = is_nk & (klrg1 > 0.5)
    else:
        nk_dim = nk_bright = nk_adapt = is_nk
    # monocyte subsets (within Myeloid)
    is_myeloid = ct.str.contains("Myeloid", regex=False) | ct.str.contains("Mono", regex=False)
    cd14 = col(a, "CD14"); cd16 = col(a, "FCGR3A"); lyz = col(a, "LYZ")
    if is_myeloid.sum() > 0:
        mono_class = is_myeloid & (cd14 > 1.0) & (cd16 < 0.5)
        mono_nc = is_myeloid & (cd14 < 0.5) & (cd16 > 0.5)
        mono_int = is_myeloid & (cd14 > 0.5) & (cd16 > 0.5)
    else:
        mono_class = mono_nc = mono_int = is_myeloid
    # DC subsets
    clec9a = col(a, "CLEC9A"); cd1c = col(a, "CD1C"); clec4c = col(a, "CLEC4C"); lilra4 = col(a, "LILRA4")
    if is_myeloid.sum() > 0:
        cdc1 = is_myeloid & (clec9a > 0.5)
        cdc2 = is_myeloid & (cd1c > 0.5)
        pdc = is_myeloid & ((clec4c > 0.5) | (lilra4 > 0.5))
    else:
        cdc1 = cdc2 = pdc = is_myeloid
    # group/disease
    grp = a.obs["tissue"] if "tissue" in a.obs else a.obs["disease"] if "disease" in a.obs else a.obs["group"] if "group" in a.obs else pd.Series("NA", index=a.obs.index)

    for mask, label in [
        (nk_dim, "NK_CD56dim"), (nk_bright, "NK_CD56bright"), (nk_adapt, "NK_adaptive"),
        (mono_class, "Mono_classical"), (mono_nc, "Mono_nonclassical"), (mono_int, "Mono_intermediate"),
        (cdc1, "cDC1"), (cdc2, "cDC2"), (pdc, "pDC")]:
        for g in grp.unique():
            n = (mask.values & (grp == g).values).sum()
            rows.append((name, g, label, n))
res = pd.DataFrame(rows, columns=["dataset", "group", "subset", "n_cells"])
print("=== 固有免疫亚群细胞数（跨数据集×分组）===")
print(res.pivot_table(index="subset", columns=["dataset","group"], values="n_cells", aggfunc="sum").fillna(0).astype(int).to_string())
res.to_csv(f"{OUT}/results/innate_subset_counts.csv", index=False)
print("\n已保存 innate_subset_counts.csv")
