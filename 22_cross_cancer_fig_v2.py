#!/usr/bin/env python3
"""
F7_v2 — Cross-cancer immune-evasion ligand comparison (4 cancer types):
CLL B cells (GSE161610, LN, 10x) | AML malignant (GSE116256, 10x) |
CML tumor (PMID 37919606, 10x) | MM plasma (GSE117156, MARS-seq).
Mean log1p(CPM). Honest annotation of platform caveat.
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "/home/zhh2/CLL-ScRNA-seqdatasets/NK_analysis"
plt.rcParams.update({"font.size": 9, "axes.labelsize": 11, "axes.labelweight": "bold",
                     "axes.titleweight": "bold", "xtick.labelsize": 9, "ytick.labelsize": 9,
                     "legend.fontsize": 8, "figure.dpi": 200, "pdf.fonttype": 42})

LIGANDS = ["CD200", "CD70", "TNFSF9", "CD86", "LGALS9", "PVR"]
# mean log1p(CPM) — CLL/CML from cross_cancer_CLL_vs_CML.tsv, AML/MM from 21_cross_cancer_process.py
cll = [0.341, 0.140, 0.396, 0.112, 0.214, 0.0003]
aml = [0.083, 0.119, 0.049, 0.017, 0.124, 0.059]
cml = [0.000, 0.009, 0.074, 0.000, 0.024, 0.031]
mm  = [0.036, 0.003, 0.001, 0.044, 0.012, 0.001]

cancers = ["CLL", "AML", "CML", "MM"]
colors = {"CLL": "#d84b6a", "AML": "#e8912a", "CML": "#1f6fb2", "MM": "#17a2b8"}

fig, ax = plt.subplots(figsize=(9, 4.5))
x = np.arange(len(LIGANDS)); w = 0.19
for i, (name, vals) in enumerate(zip(cancers, [cll, aml, cml, mm])):
    off = (i - 1.5) * w
    ax.bar(x + off, vals, w, label=name, color=colors[name], edgecolor="k", lw=0.5)
    # value labels
    for j, v in enumerate(vals):
        ax.text(x[j] + off, v + 0.012, f"{v:.2f}", ha="center", fontsize=6.5, rotation=0)

ax.set_xticks(x); ax.set_xticklabels(LIGANDS, fontweight="bold")
ax.set_xlabel("Ligand")
ax.set_ylabel("mean log1p(CPM)")
ax.set_title("Immune-evasion ligands across 4 hematologic malignancies")
ax.legend(frameon=False, loc="upper right")
for s in ["top", "right"]: ax.spines[s].set_visible(False)

# data source + caveat
ax.text(0.02, -0.22,
        "CLL: GSE161610 (LN, 10x) | AML: GSE116256 (van Galen, 10x) | "
        "CML: PMID 37919606 (10x) | MM: GSE117156 (Ledergor, MARS-seq)\n"
        "Caveat: MM (MARS-seq) is shallower than 10x; MM absolute values may be underestimated.",
        transform=ax.transAxes, fontsize=6.5, color="gray")

plt.tight_layout()
plt.savefig(f"{OUT}/figures/F7_cross_cancer_v2.png", bbox_inches="tight")
plt.savefig(f"{OUT}/figures/F7_cross_cancer_v2.pdf", bbox_inches="tight")
plt.close()
print("saved F7_cross_cancer_v2.png/.pdf")

# save combined table
tab = pd.DataFrame({"CLL": cll, "AML": aml, "CML": cml, "MM": mm}, index=LIGANDS)
tab.round(4).to_csv(f"{OUT}/results/cross_cancer_4types.tsv", sep="\t")
print(tab.round(3).to_string())
