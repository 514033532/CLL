#!/usr/bin/env python3
"""
Step 03 — GSE161610 (PB/LN paired) independent validation.

Targets (per plan sections 5 & 7):
  1. Load 10 samples (5 patients x PB/LN), QC, cluster, annotate.
  2. CLL_NK module score on NK cells (AddModuleScore of the CNI 12-gene set —
     the exact "10 marker genes" are not enumerated in Liu et al., so the CNI
     12-gene core set is used as the CLL_NK proxy, documented here).
  3. Paired PB vs LN comparison of CLL_NK score (NK cells).
  4. HLA-E (malignant B) and KLRC1/NKG2A (NK) — LN vs PB, and their correlation.

Patient pairing: PB2642~LN2641, PB4174~LN4177, PB4522~LN4521, PB5208~LN5209,
PB5585~LN5586  ->  patient id = numeric part (2642,4174,4522,5208,5585).
"""
import scanpy as sc
import pandas as pd
import numpy as np
import scipy.io as sio
import gzip, io, glob, os, warnings
warnings.filterwarnings("ignore")

OUT = "/home/zhh2/CLL-ScRNA-seqdatasets/NK_analysis"
DST = f"{OUT}/data/GSE161610"
sc.settings.verbosity = 1

# ---- 1. load all samples ---------------------------------------------------
def load_sample(base, sample_id, tissue, patient):
    mtx = sio.mmread(f"{DST}/{base}_matrix.mtx.gz").tocsr().T      # cells x genes
    fea = pd.read_csv(f"{DST}/{base}_features.tsv.gz", sep="\t", header=None)
    sym = fea[1].astype(str)
    seen, uniq = {}, []
    for s in sym:
        if s in seen:
            seen[s] += 1; uniq.append(f"{s}.{seen[s]}")
        else:
            seen[s] = 0; uniq.append(s)
    bc = pd.read_csv(f"{DST}/{base}_barcodes.tsv.gz", sep="\t", header=None)[0].astype(str)
    a = sc.AnnData(mtx, var=pd.DataFrame(index=uniq))
    a.var["gene"] = sym.values
    a.obs["sample"] = sample_id
    a.obs["tissue"] = tissue
    a.obs["patient"] = patient
    a.obs_names = [f"{sample_id}_{b}" for b in bc]
    return a

pairs = {"2642": ("GSM4910769_PB2642","GSM4910774_LN2641"),
         "4174": ("GSM4910770_PB4174","GSM4910775_LN4177"),
         "4522": ("GSM4910771_PB4522","GSM4910776_LN4521"),
         "5208": ("GSM4910772_PB5208","GSM4910777_LN5209"),
         "5585": ("GSM4910773_PB5585","GSM4910778_LN5586")}

adatas = []
for p, (pb, ln) in pairs.items():
    adatas.append(load_sample(pb, f"PB{p}", "PB", p))
    adatas.append(load_sample(ln, f"LN{p}", "LN", p))
adata = sc.concat(adatas, join="outer", index_unique=None)
adata.var_names_make_unique()
print(f"loaded: {adata.n_obs} cells x {adata.n_vars} genes")

# ---- 2. QC -----------------------------------------------------------------
adata.var["mito"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mito"], inplace=True, log1p=False)
sc.pp.filter_cells(adata, min_genes=200)
adata = adata[adata.obs["pct_counts_mito"] < 20].copy()
sc.pp.filter_genes(adata, min_cells=3)
print(f"after QC: {adata.n_obs} cells")

# per-sample doublets
adata.obs["predicted_doublet"] = False
try:
    import scrublet as scr
    adata.obs["doublet_score"] = np.nan
    for s in adata.obs["sample"].unique():
        idx = adata.obs["sample"] == s
        if idx.sum() < 50:
            continue
        raw_s = adata[idx].X.tocsr().T
        d = scr.Scrublet(raw_s, expected_doublet_rate=0.06, random_state=0)
        ds, pred = d.scrub_doublets(verbose=False)
        adata.obs.loc[idx, "doublet_score"] = ds
        adata.obs.loc[idx, "predicted_doublet"] = pred.astype(bool)
    print("doublets:", int(adata.obs["predicted_doublet"].sum()))
except Exception as e:
    print("scrublet skip:", e)

# ---- 3. normalize + cluster ------------------------------------------------
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers["log1p"] = adata.X.copy()

sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3")
adata_hvg = adata[:, adata.var["highly_variable"]].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, n_comps=50, svd_solver="randomized", random_state=0)
sc.pp.neighbors(adata_hvg, n_pcs=30, n_neighbors=15)
sc.tl.umap(adata_hvg, random_state=0)
sc.tl.leiden(adata_hvg, resolution=0.8, random_state=0)
adata.obs["leiden"] = adata_hvg.obs["leiden"].astype(str)
adata.obsm["X_umap"] = adata_hvg.obsm["X_umap"]

# ---- 4. annotate clusters (CD3G-based NK rule, consistent w/ GSE111014) ----
markers = ["CD3D","CD3E","CD3G","CD8A","CD4","NKG7","GNLY","NCAM1","FCGR3A",
           "KLRD1","KLRF1","MS4A1","CD19","CD79A","CD5","CD14","LYZ","FCN1",
           "PPBP","LILRA4","CLEC4C","MZB1"]
markers = [g for g in markers if g in adata.var_names]
sub = adata[:, markers]
per_cl = {}
for cl in adata.obs["leiden"].unique():
    per_cl[cl] = np.asarray(sub[adata.obs["leiden"]==cl].layers["log1p"].mean(axis=0)).ravel()
mk = pd.DataFrame(per_cl, index=markers).T

celltype = {}
for cl in mk.index:
    r = mk.loc[cl]
    if (r["NKG7"] > 1.0 and r["CD3G"] < 0.5 and
        (r["FCGR3A"] > 0.5 or r["KLRF1"] > 0.5 or r["NCAM1"] > 0.15)):
        celltype[cl] = "NK"
    elif r["PPBP"] > 1.0:
        celltype[cl] = "Platelet"
    elif r["CLEC4C"] > 0.5:   # BDCA2 is pDC-specific (LILRA4 also marks B cells)
        celltype[cl] = "pDC"
    elif r["CD14"] > 0.5 or r["LYZ"] > 1.0 or r["FCN1"] > 0.5:
        celltype[cl] = "Myeloid"
    elif r["CD79A"] > 0.5 or r["MS4A1"] > 0.5 or r["CD19"] > 0.5 or r["CD5"] > 1.0:
        celltype[cl] = "B"
    elif r["CD3G"] > 0.5:
        celltype[cl] = "CD8_T" if r["CD8A"] >= r["CD4"] else "CD4_T"
    else:
        celltype[cl] = "Other"
adata.obs["celltype"] = adata.obs["leiden"].map(celltype).astype("category")
mk["celltype"] = mk.index.map(celltype)
mk.to_csv(f"{OUT}/results/gse161610_annotation_qc.csv")
print("\ncelltype counts:")
print(adata.obs["celltype"].value_counts())

# ---- 5. CLL_NK score (CNI 12 genes) + balance score on NK ------------------
CNI = ["SLAMF7","CST7","NKG7","CD63","EFHD2","DOK2","FCGR3A","TRBC1",
       "ASCL2","IL2RG","CD81","GZMB"]
FUNC = ["NKG7","ASCL2","GZMB","PRF1","GNLY","IL2RG"]
EXH  = ["TIGIT","CD96","LAG3","BTLA","HAVCR2","KLRC1"]
for gl, nm in [(CNI,"cll_nk_score"), (FUNC,"func_score"), (EXH,"exh_score")]:
    gl = [g for g in gl if g in adata.var_names]
    sc.tl.score_genes(adata, gl, score_name=nm, random_state=0)

nk = adata[adata.obs["celltype"] == "NK"].copy()
print(f"\nNK cells: {nk.n_obs}")
nk_counts = (nk.obs.groupby(["patient","tissue"], observed=True)
             .size().unstack().fillna(0).astype(int))
print(nk_counts)

# per-sample NK CLL_NK score + KLRC1
nk_sum = (nk.obs.groupby(["patient","tissue"], observed=True)
          .agg(cll_nk_mean=("cll_nk_score","mean"),
               nk_mean=("func_score","mean"),
               exh_mean=("exh_score","mean"),
               n=("cll_nk_score","size")).reset_index())
nk_sum.to_csv(f"{OUT}/results/gse161610_nk_cllnk.csv", index=False)
print("\nCLL_NK score per sample (NK cells):")
print(nk_sum.round(3).to_string(index=False))

# ---- 6. paired PB vs LN (CLL_NK score) -------------------------------------
wide = nk_sum.pivot(index="patient", columns="tissue", values="cll_nk_mean")
wide["diff_LN_minus_PB"] = wide["LN"] - wide["PB"]
wide.to_csv(f"{OUT}/results/gse161610_paired_cllnk.csv")
print("\npaired CLL_NK score (PB vs LN):")
print(wide.round(3).to_string())
if wide.dropna().shape[0] >= 3:
    from scipy.stats import wilcoxon, ttest_rel
    w = wilcoxon(wide["LN"], wide["PB"])
    t = ttest_rel(wide["LN"], wide["PB"])
    print(f"  paired Wilcoxon(LN vs PB): stat={w.statistic:.2f}, p={w.pvalue:.4f}")
    print(f"  paired t-test: t={t.statistic:.2f}, p={t.pvalue:.4f}")

# ---- 7. HLA-E (B) and KLRC1 (NK) — LN vs PB --------------------------------
b = adata[adata.obs["celltype"] == "B"].copy()
hla_e = b.obs.assign(hlae=b[:, "HLA-E"].X.toarray().ravel() if "HLA-E" in b.var_names else 0)
if "HLA-E" in b.var_names:
    b.obs["hla_e"] = np.asarray(b.layers["log1p"][:, list(b.var_names).index("HLA-E")].toarray()).ravel()
    hla_sum = (b.obs.groupby(["patient","tissue"], observed=True)["hla_e"]
               .mean().unstack())
    hla_sum["diff_LN_minus_PB"] = hla_sum["LN"] - hla_sum["PB"]
    hla_sum.to_csv(f"{OUT}/results/gse161610_hla_e.csv")
    print("\nHLA-E (B cells) PB vs LN:")
    print(hla_sum.round(3).to_string())

if "KLRC1" in nk.var_names:
    nk.obs["klrc1"] = np.asarray(nk.layers["log1p"][:, list(nk.var_names).index("KLRC1")].toarray()).ravel()
    klrc1_sum = (nk.obs.groupby(["patient","tissue"], observed=True)["klrc1"]
                 .mean().unstack())
    klrc1_sum["diff_LN_minus_PB"] = klrc1_sum["LN"] - klrc1_sum["PB"]
    klrc1_sum.to_csv(f"{OUT}/results/gse161610_klrc1.csv")
    print("\nKLRC1/NKG2A (NK cells) PB vs LN:")
    print(klrc1_sum.round(3).to_string())
    # NKG2A vs CLL_NK correlation (per-sample)
    both = nk_sum.merge(nk.obs.groupby(["patient","tissue"], observed=True)["klrc1"].mean().reset_index(),
                        on=["patient","tissue"])
    from scipy.stats import spearmanr
    rho, p = spearmanr(both["klrc1"], both["cll_nk_mean"])
    print(f"\nSpearman(KLRC1 ~ CLL_NK score): rho={rho:.3f}, p={p:.3f} (n={len(both)})")

# ---- 8. save ---------------------------------------------------------------
adata.X = adata.layers["log1p"].copy()
adata.write_h5ad(f"{OUT}/data/gse161610_processed.h5ad")
sc.pl.umap(adata, color=["celltype","tissue"], show=False, save="_gse161610.png")
print("\nDONE")
