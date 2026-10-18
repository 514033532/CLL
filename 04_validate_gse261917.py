#!/usr/bin/env python3
"""
Step 07 — GSE261917 independent validation of LN CD200/CD70 upregulation.

GSE261917 = FL/CLL proliferative LN microenvironment (3 CLL, 3 FL, 3 reactive LN).
Load CD45+ (immune) fractions, identify B cells, compare CD200/CD70/TNFSF9/CD86
across disease (Reactive vs FL vs CLL) to test CLL-specificity.
"""
import scanpy as sc
import pandas as pd
import numpy as np
import scipy.io as sio
import gzip, glob, os, warnings
warnings.filterwarnings("ignore")

OUT = "/home/zhh2/CLL-ScRNA-seqdatasets/NK_analysis"
EXT = "/home/zhh2/CLL-ScRNA-seqdatasets/extracted/GSE261917_RAW"
sc.settings.verbosity = 1

def load_10x(sample_id):
    d = EXT
    mtx = sio.mmread(gzip.open(glob.glob(f"{d}/{sample_id}_matrix.mtx.gz")[0], "rb")).tocsr().T
    fea = pd.read_csv(gzip.open(glob.glob(f"{d}/{sample_id}_features.tsv.gz")[0], "rt"),
                      sep="\t", header=None)
    sym = fea.iloc[:, 1].astype(str)
    seen, uniq = {}, []
    for s in sym:
        if s in seen:
            seen[s] += 1; uniq.append(f"{s}.{seen[s]}")
        else:
            seen[s] = 0; uniq.append(s)
    bc = pd.read_csv(gzip.open(glob.glob(f"{d}/{sample_id}_barcodes.tsv.gz")[0], "rt"),
                     sep="\t", header=None)[0].astype(str)
    a = sc.AnnData(mtx, var=pd.DataFrame(index=uniq))
    a.var["gene"] = sym.values
    a.obs["sample"] = sample_id
    a.obs_names = [f"{sample_id}_{b}" for b in bc]
    return a

# sample -> disease (use CD45+ and combined fractions; skip pure CD45- stroma)
samples = {
    "GSM8154135_CLL1_CD45+": "CLL", "GSM8154137_CLL2_CD45+": "CLL",
    "GSM8154139_CLL3_CD45+_and_CD45-": "CLL",
    "GSM8154130_FL1_CD45+": "FL", "GSM8154132_FL2_CD45+": "FL",
    "GSM8154134_FL3_CD45+_and_CD45-": "FL",
    "GSM8154127_Reactive2_CD45+": "Reactive",
    "GSM8154129_Reactive3_CD45+_and_CD45-": "Reactive",
}
# verify sample ids exist
existing = {}
for s in samples:
    if glob.glob(f"{EXT}/{s}_matrix.mtx.gz"):
        existing[s] = samples[s]
print(f"loading {len(existing)} samples:", existing)

adatas = []
for s, d in existing.items():
    a = load_10x(s); a.obs["disease"] = d
    adatas.append(a)
adata = sc.concat(adatas, join="outer", index_unique=None)
adata.var_names_make_unique()
adata.var["mito"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mito"], inplace=True, log1p=False)
sc.pp.filter_cells(adata, min_genes=200)
adata = adata[adata.obs["pct_counts_mito"] < 20].copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers["log1p"] = adata.X.copy()
print(f"after QC: {adata.n_obs} cells")

genes = list(adata.var_names)
def col(g): return np.asarray(adata.layers["log1p"][:, genes.index(g)].toarray()).ravel() if g in genes else None

# B cell identification (marker-based)
is_b = (col("CD79A") > 0.5) | (col("MS4A1") > 0.5) | (col("CD19") > 0.5)
is_t = (col("CD3D") > 0.5) | (col("CD3E") > 0.5)
is_b = is_b & (~is_t)
print(f"\nB cells: {is_b.sum()} / {adata.n_obs}")

# CD200/CD70/TNFSF9/CD86 on B cells, by disease
print("\n=== B 细胞 CD200/CD70/TNFSF9/CD86 表达, 按疾病 ===")
genes_oi = ["CD200","CD70","TNFSF9","CD86","CD79A","MS4A1","CD200R1"]
for g in genes_oi:
    if col(g) is None: continue
    row = []
    for d in ["Reactive","FL","CLL"]:
        m = is_b & (adata.obs["disease"]==d).values
        if m.sum() < 5: continue
        v = col(g)[m]
        row.append(f"{d}: frac={(v>0).mean():.2f}/mean={v.mean():.2f} (n={m.sum()})")
    print(f"  {g:9s} " + "  ".join(row))

# 统计检验: CLL vs Reactive (Mann-Whitney)
from scipy.stats import mannwhitneyu
print("\n=== Mann-Whitney (CLL vs Reactive B cells) ===")
for g in ["CD200","CD70","TNFSF9","CD86"]:
    if col(g) is None: continue
    v_c = col(g)[is_b & (adata.obs["disease"]=="CLL").values]
    v_r = col(g)[is_b & (adata.obs["disease"]=="Reactive").values]
    v_f = col(g)[is_b & (adata.obs["disease"]=="FL").values]
    u_cr, p_cr = mannwhitneyu(v_c, v_r)
    u_cf, p_cf = mannwhitneyu(v_c, v_f)
    print(f"  {g:8s}: CLL vs Reac p={p_cr:.2e}, CLL vs FL p={p_cf:.2e}  "
          f"(CLL mean={v_c.mean():.3f}, Reac={v_r.mean():.3f}, FL={v_f.mean():.3f})")

adata.X = adata.layers["log1p"].copy()
adata.write_h5ad(f"{OUT}/data/gse261917_processed.h5ad")
print("\nDONE")
