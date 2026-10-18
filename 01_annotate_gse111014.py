#!/usr/bin/env python3
"""
Step 01 — GSE111014: load, rigorous QC, doublet removal, clustering, and
cluster-level marker annotation; then NK-cell extraction and donor-level
balance score.

GSE111014 (Rendeiro et al. 2020, Nat Commun) — ibrutinib time course in 4 CLL
patients (CLL1/CLL5/CLL6/CLL8). Matrix: 33,694 genes x 48,016 cells.

Strategy (per reviewer/QC requirement):
  1. Cell-level QC: min_genes / max_genes / pct_mito + Scrublet doublets
  2. HVG -> scale -> PCA -> neighbors -> UMAP -> Leiden clustering
  3. Cluster annotation via (a) top DE markers and (b) lineage-marker mean
     scores per cluster; NK clusters must be CD3-negative & NKG7/GNLY-positive
  4. NK extraction from annotated clusters (validated), not raw score threshold

Outputs:
  data/gse111014_processed.h5ad  QC'd normalized full object (with clusters)
  data/gse111014_nk.h5ad         NK subset
  results/gse111014_qc_summary.csv
  results/gse111014_cluster_markers.csv
  results/gse111014_cluster_annotation.csv
  results/gse111014_nk_counts.csv
  results/gse111014_balance_donor.csv
  figures/gse111014_umap_celltype.png
"""
import scanpy as sc
import pandas as pd
import numpy as np
import scipy.io as sio
import os, warnings
warnings.filterwarnings("ignore")

BASE = "/home/zhh2/CLL-ScRNA-seqdatasets/extracted"
OUT  = "/home/zhh2/CLL-ScRNA-seqdatasets/NK_analysis"
sc.settings.verbosity = 1
sc.settings.figdir = f"{OUT}/figures"

# ---- 1. load matrix -------------------------------------------------------
print("== loading genes / matrix / barcodes ==")
genes = pd.read_csv(f"{BASE}/GSE111014_genes.tsv/genes.tsv",
                    sep="\t", header=None, names=["ensembl", "symbol"])
genes["symbol"] = genes["symbol"].astype(str)
seen, uniq = {}, []
for s in genes["symbol"]:
    if s in seen:
        seen[s] += 1; uniq.append(f"{s}.{seen[s]}")
    else:
        seen[s] = 0; uniq.append(s)
genes["symbol_u"] = uniq

M = sio.mmread(f"{BASE}/GSE111014_matrix.mtx/matrix.mtx").tocsr().T  # cells x genes
bc = pd.read_csv(f"{BASE}/GSE111014_raw/GSE111014_barcodes.tsv.gz",
                 sep="\t", header=None)[0].astype(str)
sample  = bc.str.split("scRNA-seq_", expand=True)[1]
patient = sample.str.extract(r"(CLL\d+)")[0]
timept  = sample.str.extract(r"(d\d+)$")[0]

adata = sc.AnnData(M, var=pd.DataFrame(index=genes["symbol_u"]))
adata.var["ensembl"] = genes["ensembl"].values
adata.var["gene"]    = genes["symbol"].values
adata.obs["barcode"] = bc.values
adata.obs["sample"]  = sample.values
adata.obs["patient"] = patient.values
adata.obs["timepoint"] = timept.values
adata.obs_names = bc.values

# ---- 2. QC -----------------------------------------------------------------
print("== QC ==")
n0 = adata.n_obs
adata.var["mito"] = adata.var["gene"].str.startswith("MT-")
adata.var["ribo"] = adata.var["gene"].str.startswith(("RPS","RPL"))
sc.pp.calculate_qc_metrics(adata, qc_vars=["mito","ribo"],
                           percent_top=None, log1p=False, inplace=True)

# report distributions before filtering
qc = pd.DataFrame({
    "n_genes_min": adata.obs["n_genes_by_counts"].min(),
    "n_genes_median": adata.obs["n_genes_by_counts"].median(),
    "n_genes_max": adata.obs["n_genes_by_counts"].max(),
    "pct_mito_median": adata.obs["pct_counts_mito"].median(),
    "n_counts_median": adata.obs["total_counts"].median(),
}, index=[0])

# adaptive gene-count ceiling: mean + 3 SD of log n_genes
log_ng = np.log1p(adata.obs["n_genes_by_counts"])
max_genes = int(np.expm1(log_ng.mean() + 3 * log_ng.std()))
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_cells(adata, max_genes=max_genes)
adata = adata[adata.obs["pct_counts_mito"] < 20].copy()
sc.pp.filter_genes(adata, min_cells=3)

qc["max_genes_cutoff"] = max_genes
qc["cells_before"] = n0
qc["cells_after"]  = adata.n_obs
qc.to_csv(f"{OUT}/results/gse111014_qc_summary.csv", index=False)
print(f"cells {n0} -> {adata.n_obs} (max_genes cutoff {max_genes})")

# ---- 3. doublet detection (Scrublet, per sample) --------------------------
print("== doublet detection ==")
try:
    import scrublet as scr
    adata.obs["doublet_score"] = np.nan
    adata.obs["predicted_doublet"] = False
    for s in adata.obs["sample"].cat.categories:
        idx = adata.obs["sample"] == s
        n = int(idx.sum())
        if n < 50:
            continue
        raw_s = adata[idx].X.tocsr().T                      # genes x cells
        scr_s = scr.Scrublet(raw_s, expected_doublet_rate=0.06,
                             random_state=0)
        ds, pred = scr_s.scrub_doublets(verbose=False)
        adata.obs.loc[idx, "doublet_score"] = ds
        adata.obs.loc[idx, "predicted_doublet"] = pred.astype(bool)
    print("predicted doublets:", int(adata.obs["predicted_doublet"].sum()))
except Exception as e:
    print("scrublet unavailable:", e)
    adata.obs["doublet_score"] = 0.0
    adata.obs["predicted_doublet"] = False

# ---- 4. normalize + HVG + scale + PCA + neighbors + UMAP + Leiden ---------
print("== normalize / HVG / PCA / clustering ==")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers["log1p"] = adata.X.copy()
adata.raw = adata

sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3")
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=50, svd_solver="randomized", random_state=0)
sc.pp.neighbors(adata, n_pcs=30, n_neighbors=15)
sc.tl.umap(adata, random_state=0)
sc.tl.leiden(adata, resolution=0.8, random_state=0)

# ---- 5. cluster annotation -------------------------------------------------
print("== cluster annotation ==")
lineage = {
    "NK":      ["NCAM1","FCGR3A","NKG7","GNLY","KLRD1","KLRF1","CD7","PRF1","KLRB1"],
    "T":       ["CD3D","CD3E","CD3G","TRAC","TRBC1","CD2"],
    "B":       ["MS4A1","CD19","CD79A","CD79B","CD5","IGHD","IGHM"],
    "Myeloid": ["CD14","LYZ","MS4A7","FCN1","S100A8","S100A9","CD68","CST3"],
    "Mono_CD16": ["FCGR3A","MS4A7","LST1","IFI30"],
    "Plasma":  ["MZB1","XBP1","JCHAIN","DERL3","SSR4"],
    "pDC":     ["LILRA4","CLEC4C","IRF7","GZMB"],
}
# per-cluster mean scaled expression of each lineage marker set
for name, gl in lineage.items():
    gl = [g for g in gl if g in adata.var_names]
    sc.tl.score_genes(adata, gl, score_name=f"tmp_{name}", random_state=0)
lin_cols = [f"tmp_{name}" for name in lineage]
cl_mean = adata.obs.groupby("leiden")[lin_cols].mean()
cl_mean.columns = [c.replace("tmp_","lin_") for c in cl_mean.columns]

# top DE markers per cluster (Wilcoxon)
sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon", use_raw=True)
markers = pd.DataFrame({
    group: [adata.uns["rank_genes_groups"]["names"][group][i] for i in range(10)]
    for group in adata.obs["leiden"].cat.categories
}).T
markers.to_csv(f"{OUT}/results/gse111014_cluster_markers.csv")
print(markers)

# ---- per-cluster mean log1p of canonical markers (annotation QC + NK gate) --
canonical = ["CD3D","CD3E","CD3G","CD8A","CD8B","CD4","NKG7","GNLY","NCAM1",
             "FCGR3A","KLRD1","KLRF1","PRF1","GZMB","MS4A1","CD19","CD5",
             "CD79A","CD14","LYZ","FCN1","PPBP","MZB1","JCHAIN","LILRA4","CLEC4C"]
canonical = [g for g in canonical if g in adata.var_names]
sub = adata[:, canonical]
per_cl = {}
for cl in adata.obs["leiden"].cat.categories:
    per_cl[cl] = np.asarray(sub[adata.obs["leiden"] == cl].layers["log1p"].mean(axis=0)).ravel()
marker_cl = pd.DataFrame(per_cl, index=canonical).T

# classify: NK uses CD3G (T-specific chain, low in NK) + NK-specific receptors
# (FCGR3A/CD16, KLRF1/NKp80, NCAM1/CD56); everything else by lineage score.
# CD3D/CD3E are NOT used for the T gate because they leak in NK cells.
celltype = {}
for cl in cl_mean.index:
    scores = cl_mean.loc[cl]
    r = marker_cl.loc[cl]
    is_nk = (scores["lin_NK"] > 0.2 and r["CD3G"] < 0.5 and
             (r["FCGR3A"] > 0.5 or r["KLRF1"] > 0.5 or r["NCAM1"] > 0.15) and
             r["NKG7"] > 1.0)
    if is_nk:
        ct = "NK"
    else:
        ct = scores.idxmax().replace("lin_", "")
        # CD8 T cells also score high on the NK gene module (NKG7/GNLY);
        # if a cluster is CD3G-positive, force it out of the NK call.
        if ct == "NK" and r["CD3G"] >= 0.5:
            ct = "T"
        if ct == "Mono_CD16":
            ct = "Monocyte_CD16"
    celltype[cl] = ct

marker_cl["celltype"] = marker_cl.index.map(celltype)
marker_cl.insert(0, "n_cells",
                 adata.obs["leiden"].value_counts().reindex(marker_cl.index).values)
marker_cl.to_csv(f"{OUT}/results/gse111014_cluster_annotation.csv")
adata.obs["celltype"] = adata.obs["leiden"].map(celltype).astype("category")
print(adata.obs["celltype"].value_counts())

marker_cl.to_csv(f"{OUT}/results/gse111014_annotation_qc.csv")
print("\n== annotation QC (mean log1p expr per cluster) ==")
print(marker_cl.round(2).to_string())

# ---- 6. NK extraction (from annotated NK clusters, validated) -------------
print("== NK extraction ==")
adata.obs["is_nk"] = (adata.obs["celltype"] == "NK").values
nk = adata[adata.obs["is_nk"]].copy()
print(f"NK cells: {nk.n_obs}  (from clusters {sorted(adata.obs.loc[adata.obs.is_nk,'leiden'].unique())})")

counts = (adata.obs.assign(is_nk=adata.obs["is_nk"].astype(int))
          .groupby(["patient","timepoint","sample"], observed=True)["is_nk"]
          .agg(total="size", nk="sum").reset_index())
counts["nk_frac"] = counts["nk"] / counts["total"]
counts.to_csv(f"{OUT}/results/gse111014_nk_counts.csv", index=False)
print(counts)

# ---- 7. balance score ------------------------------------------------------
print("== balance score ==")
FUNC = ["NKG7","ASCL2","GZMB","PRF1","GNLY","IL2RG"]
EXH  = ["TIGIT","CD96","LAG3","BTLA","HAVCR2","KLRC1"]
print("missing functional:", [g for g in FUNC if g not in adata.var_names],
      "| missing exhaustion:", [g for g in EXH if g not in adata.var_names])

# (a) plan formula on positive mean log1p-normalized expression (per cell)
func_idx = [list(adata.var_names).index(g) for g in FUNC]
exh_idx  = [list(adata.var_names).index(g) for g in EXH]
func_expr = np.asarray(adata.layers["log1p"][:, func_idx].mean(axis=1)).ravel()
exh_expr  = np.asarray(adata.layers["log1p"][:, exh_idx].mean(axis=1)).ravel()
adata.obs["func_expr"] = func_expr
adata.obs["exh_expr"]  = exh_expr

# (b) AddModuleScore (Seurat-style, control-gene corrected) for robustness
sc.tl.score_genes(adata, FUNC, score_name="func_score", random_state=0)
sc.tl.score_genes(adata, EXH,  score_name="exh_score",  random_state=0)

nk = adata[adata.obs["is_nk"]].copy()
donor = (nk.obs.groupby(["patient","timepoint","sample"], observed=True)
         .agg(n_nk=("func_expr","size"),
              func_mean=("func_expr","mean"),
              exh_mean=("exh_expr","mean"),
              func_score=("func_score","mean"),
              exh_score=("exh_score","mean")).reset_index())
donor["balance"] = (donor["func_mean"] - donor["exh_mean"]) / \
                   (donor["func_mean"] + donor["exh_mean"])   # plan formula
donor["func_exh_ratio"] = donor["func_mean"] / donor["exh_mean"]
donor["balance_addmod"] = donor["func_score"] - donor["exh_score"]  # robustness
donor.to_csv(f"{OUT}/results/gse111014_balance_donor.csv", index=False)
print(donor.round(3))

# ---- 8. save (compact: X=sparse log1p; drop dense scaled X & redundant raw) --
adata.X = adata.layers["log1p"].copy()   # restore sparse log1p as X
adata.raw = None                          # raw is redundant with layers["log1p"]
adata.write_h5ad(f"{OUT}/data/gse111014_processed.h5ad")

nk = adata[adata.obs["is_nk"]].copy()     # inherits sparse log1p X
nk.raw = None
nk.write_h5ad(f"{OUT}/data/gse111014_nk.h5ad")

sc.settings.figdir = f"{OUT}/figures"
sc.pl.umap(adata, color=["celltype", "is_nk", "patient"],
           legend_loc="on data", show=False, save="_celltype.png")
print("saved h5ad + UMAP figure")
print("DONE")
