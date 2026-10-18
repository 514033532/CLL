#!/usr/bin/env python3
"""
Fig 7 — Cross-cancer comparison: CLL vs CML immune-evasion ligand program.
CLL B cells (GSE161610, LN) vs CML tumor cells (Krishnan 2024, PMID 37919606).
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "/home/zhh2/CLL-ScRNA-seqdatasets/NK_analysis"
plt.rcParams.update({"font.size": 9, "axes.labelsize": 11, "axes.labelweight": "bold",
                     "axes.titleweight": "bold", "xtick.labelsize": 9, "ytick.labelsize": 9,
                     "legend.fontsize": 9, "figure.dpi": 200, "pdf.fonttype": 42})

# CLL B cells (LN) vs CML tumor cells 配体表达
ligands = ["CD200", "CD70", "TNFSF9", "CD86", "LGALS9", "PVR"]
cll = [0.341, 0.140, 0.396, 0.112, 0.214, 0.000]
cml = [0.000, 0.009, 0.074, 0.000, 0.024, 0.031]

fig, ax = plt.subplots(figsize=(7, 4))
x = np.arange(len(ligands)); w = 0.38
b1 = ax.bar(x - w/2, cll, w, label="CLL B cells (LN)", color="#d84b6a", edgecolor="k", lw=0.5)
b2 = ax.bar(x + w/2, cml, w, label="CML tumor cells", color="#1f6fb2", edgecolor="k", lw=0.5)
# 标注倍数（CLL 显著更高的）
for i, (c1, c2, g) in enumerate(zip(cll, cml, ligands)):
    if c1 > c2 * 2 and c1 > 0.05:
        ax.text(i - w/2, c1 + 0.015, f"{c1:.2f}", ha="center", fontsize=7.5, color="#d84b6a")
    if c2 > c1 * 2 and c2 > 0.005:
        ax.text(i + w/2, c2 + 0.015, f"{c2:.3f}", ha="center", fontsize=7.5, color="#1f6fb2")
ax.set_xticks(x); ax.set_xticklabels(ligands, fontweight="bold")
ax.set_xlabel("Ligand")
ax.set_ylabel("mean log-normalized expression")
ax.set_title("Immune-evasion ligands: CLL vs CML")
ax.legend(frameon=False, loc="upper right")
for s in ["top", "right"]: ax.spines[s].set_visible(False)
# 数据来源标注
ax.text(0.98, 0.02, "CLL: GSE161610 | CML: PMID 37919606", transform=ax.transAxes,
        fontsize=7, color="gray", ha="right", va="bottom")

plt.tight_layout()
plt.savefig(f"{OUT}/figures/F7_cross_cancer.png", bbox_inches="tight")
plt.savefig(f"{OUT}/figures/F7_cross_cancer.pdf", bbox_inches="tight")
plt.close()
print("saved F7_cross_cancer.png/.pdf")
