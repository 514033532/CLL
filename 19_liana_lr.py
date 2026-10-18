#!/usr/bin/env python3
"""
Step 19 — Rigorous ligand-receptor re-analysis (LIANA / CellPhoneDB permutation).
Question: is CD200->CD200R1 and CD70->CD27 a SIGNIFICANT B->effector interaction
in GSE161610, and does it differ LN vs PB?

Resource: OmniPath consensus (the only curated DB including CD200-CD200R1 AND CD70-CD27).
Method: CellPhoneDB (mean-based + 1000 label permutations for significance).
"""
import warnings; warnings.filterwarnings("ignore")
import scanpy as sc
import pandas as pd
import numpy as np
import liana as li
import liana.mt as mt

OUT = "/home/zhh2/CLL-ScRNA-seqdatasets/NK_analysis"
a = sc.read_h5ad(f"{OUT}/data/gse161610_processed.h5ad")

# keep main cell types only
keep = ["B", "CD4_T", "CD8_T", "NK", "Myeloid"]
a = a[a.obs["celltype"].isin(keep)].copy()
a.obs["celltype"] = a.obs["celltype"].astype(str)

TARGET_LIGANDS = {"CD200", "CD70", "TNFSF9", "CD86"}

def run(tissue_label, adata):
    mt.cellphonedb(adata, groupby="celltype", resource_name="consensus",
                   layer="log1p", use_raw=False, expr_prop=0.1, min_cells=3,
                   n_perms=1000, seed=1337, verbose=False, inplace=True)
    r = adata.uns["liana_res"].copy()
    r["tissue"] = tissue_label
    return r

# run on LN, PB, combined
res = []
for lab, sub in [("LN", a[a.obs["tissue"] == "LN"].copy()),
                 ("PB", a[a.obs["tissue"] == "PB"].copy()),
                 ("PB+LN", a.copy())]:
    r = run(lab, sub)
    res.append(r)
    print(f"{lab}: {len(r)} interaction rows")

full = pd.concat(res, ignore_index=True)

# inspect columns
print("\ncolumns:", list(full.columns))

# focus on target ligands (ligand_complex contains the gene symbol)
lig_col = "ligand_complex"
f = full[full[lig_col].astype(str).apply(lambda s: any(g in str(s) for g in TARGET_LIGANDS))]
f = f.copy()

# keep key columns
cols = [c for c in ["tissue", "source", "target", "ligand_complex", "receptor_complex",
                    "lr_means", "cellphone_pvals", "magnitude_rank", "specificity_rank"] if c in f.columns]
f = f[cols].sort_values(["tissue", "ligand_complex", "cellphone_pvals"])

f.to_csv(f"{OUT}/results/gse161610_liana_LR.tsv", sep="\t", index=False)
print("\n=== targeted LR interactions (p<0.1 shown) ===")
show = f[f["cellphone_pvals"] < 0.1] if "cellphone_pvals" in f.columns else f
print(show.to_string(index=False))

# also report receptor mean expression per cell type (LN) for CD27/CD200R1
print("\n=== receptor mean expr (LN, log1p) ===")
ln = a[a.obs["tissue"] == "LN"]
for rec in ["CD27", "CD200R1", "CD200R1L"]:
    if rec in ln.var_names:
        g = ln.var_names.get_loc(rec)
        v = np.asarray(ln.layers["log1p"][:, g].toarray()).ravel()
        row = {ct: round(v[ln.obs["celltype"] == ct].mean(), 3) for ct in keep}
        print(f"  {rec}: {row}")
print("\nDONE")
