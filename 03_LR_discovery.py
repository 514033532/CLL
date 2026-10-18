#!/usr/bin/env python3
"""
Step 06 — Cell-cell interaction differential analysis (PB vs LN), GSE161610.

Design: map ligand-receptor (L-R) communication between cell types, and find
which axes are rewired in the CLL lymph node vs peripheral blood.
Score = mean(log1p ligand in sender) x mean(log1p receptor in receiver).

Outputs: results/gse161610_LR_interactions.csv
"""
import scanpy as sc
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings("ignore")

OUT = "/home/zhh2/CLL-ScRNA-seqdatasets/NK_analysis"
adata = sc.read_h5ad(f"{OUT}/data/gse161610_processed.h5ad")
X = adata.layers["log1p"]
genes = list(adata.var_names)

def col(g): return np.asarray(X[:, genes.index(g)].toarray()).ravel() if g in genes else None

# ---- curated immune ligand-receptor database --------------------------------
# (ligand, receptor(s), category)
LR_DB = [
    # checkpoints / inhibitory
    ("LGALS9","HAVCR2","checkpoint"),
    ("LGALS9","CD44","checkpoint"),
    ("CD274","PDCD1","checkpoint"),
    ("PVR","TIGIT","checkpoint"),
    ("NECTIN2","TIGIT","checkpoint"),
    ("HLA-E","KLRC1","checkpoint"),
    ("HLA-G","LILRB1","checkpoint"),
    ("CD47","SIRPA","checkpoint"),
    ("CD80","CTLA4","checkpoint"),
    ("CD86","CTLA4","checkpoint"),
    ("IL10","IL10RA","suppressive_cytokine"),
    ("IL10","IL10RB","suppressive_cytokine"),
    ("TGFB1","TGFBR2","suppressive_cytokine"),
    ("VEGFA","FLT1","suppressive_cytokine"),
    ("MIF","CD74","pleiotropic"),
    # stimulatory / effector-supporting
    ("IL15","IL2RB","stimulatory_cytokine"),
    ("IL15","IL2RG","stimulatory_cytokine"),
    ("IL2","IL2RB","stimulatory_cytokine"),
    ("IL18","IL18R1","stimulatory_cytokine"),
    ("IFNG","IFNGR1","stimulatory_cytokine"),
    ("CD80","CD28","costim"),
    ("CD86","CD28","costim"),
    ("CD40LG","CD40","costim"),
    ("CD70","CD27","costim"),
    ("TNFSF13B","TNFRSF13C","B_survival"),   # BAFF-BAFFR
    ("TNFSF13B","TNFRSF13B","B_survival"),   # BAFF-TACI
    ("TNFSF13","TNFRSF17","B_survival"),     # APRIL-BCMA
    # chemokines (recruitment/homing)
    ("CXCL9","CXCR3","chemokine"),
    ("CXCL10","CXCR3","chemokine"),
    ("CXCL11","CXCR3","chemokine"),
    ("CXCL12","CXCR4","chemokine"),
    ("CXCL13","CXCR5","chemokine"),
    ("CXCL16","CXCR6","chemokine"),
    ("CCL19","CCR7","chemokine"),
    ("CCL21","CCR7","chemokine"),
    ("CCL5","CCR5","chemokine"),
    ("CCL5","CCR1","chemokine"),
    ("CCL3","CCR5","chemokine"),
    ("CCL4","CCR5","chemokine"),
    ("CCL2","CCR2","chemokine"),
    ("CCL8","CCR2","chemokine"),
    ("CCL17","CCR4","chemokine"),
    ("CCL22","CCR4","chemokine"),
    ("CX3CL1","CX3CR1","chemokine"),
]

# ---- per-tissue, per-celltype mean ligand/receptor expression --------------
celltypes = [c for c in adata.obs["celltype"].cat.categories if c != "Other"]
tissues = ["PB", "LN"]

# mean log1p expression matrix: (tissue, celltype) x gene
mean_expr = {}
frac_expr = {}
for t in tissues:
    sub = adata[adata.obs["tissue"] == t]
    for ct in celltypes:
        m = sub[sub.obs["celltype"] == ct]
        if m.n_obs < 5:
            continue
        mean_expr[(t, ct)] = np.asarray(m.X.mean(axis=0)).ravel()
        frac_expr[(t, ct)] = np.asarray((m.X > 0).mean(axis=0)).ravel()

def gene_mean(t, ct, g):
    if (t, ct) not in mean_expr or g not in genes:
        return 0.0
    return mean_expr[(t, ct)][genes.index(g)]

rows = []
for lig, rec, cat in LR_DB:
    if lig not in genes or rec not in genes:
        continue
    for t in tissues:
        # find sender with max ligand, receiver with max receptor
        senders = {ct: gene_mean(t, ct, lig) for ct in celltypes}
        recvs   = {ct: gene_mean(t, ct, rec) for ct in celltypes}
        # all pairwise interaction scores
        for s in celltypes:
            for r in celltypes:
                score = senders[s] * recvs[r]
                rows.append(dict(tissue=t, ligand=lig, receptor=rec, category=cat,
                                 sender=s, receiver=r, score=score))
LR = pd.DataFrame(rows)

# ---- dominant interaction per (tissue, LR pair) ----------------------------
dom = (LR.sort_values("score", ascending=False)
       .drop_duplicates(["tissue","ligand","receptor","category"])
       .reset_index(drop=True))
dom = dom.rename(columns={"sender":"dom_sender","receiver":"dom_receiver","score":"dom_score"})

# ---- LN vs PB differential -------------------------------------------------
piv = dom.pivot_table(index=["ligand","receptor","category"],
                      columns="tissue", values="dom_score", aggfunc="max")
piv["LN_minus_PB"] = piv["LN"] - piv["PB"]
piv["log2FC"] = np.log2((piv["LN"] + 1e-4) / (piv["PB"] + 1e-4))
piv = piv.sort_values("LN_minus_PB", ascending=False)

# attach dominant sender/receiver for LN
ln_dom = dom[dom["tissue"]=="LN"].set_index(["ligand","receptor"])[["dom_sender","dom_receiver"]]
pb_dom = dom[dom["tissue"]=="PB"].set_index(["ligand","receptor"])[["dom_sender","dom_receiver"]]
out = piv.join(ln_dom.rename(columns={"dom_sender":"LN_sender","dom_receiver":"LN_receiver"}))
out = out.join(pb_dom.rename(columns={"dom_sender":"PB_sender","dom_receiver":"PB_receiver"}))
out = out.reset_index()
out.to_csv(f"{OUT}/results/gse161610_LR_interactions.csv", index=False)

print("=== LN 上调的互作轴 (LN-PB > 0.5) ===")
print(out[out["LN_minus_PB"] > 0.5][["ligand","receptor","category","LN_sender","LN_receiver","PB","LN","LN_minus_PB"]].round(2).to_string(index=False))
print("\n=== LN 下调的互作轴 (LN-PB < -0.5) ===")
print(out[out["LN_minus_PB"] < -0.5][["ligand","receptor","category","LN_sender","LN_receiver","PB","LN","LN_minus_PB"]].round(2).to_string(index=False))
print("\nDONE")
