#!/usr/bin/env python3
"""
Step 21 — Cross-cancer ligand expression: AML (GSE116256) + MM (GSE117156).
Compute mean log1p(CPM) of immune-evasion ligands on MALIGNANT/tumor cells,
same normalization as CLL/CML (normalize_total 1e4 + log1p).

AML: malignant = PredictionRefined=="malignant" (from .anno), UMI total from .anno TranscriptomeUMIs.
MM:  tumor = plasma cells (MZB1+/SDC1+/CD38+), UMI total = rowsum over full matrix.
"""
import gzip, glob, os, sys
import numpy as np
import pandas as pd

LIGANDS = ["CD200", "CD70", "TNFSF9", "CD86", "LGALS9", "PVR"]
MM_MARKERS = ["MZB1", "SDC1", "CD38", "XBP1"]

def read_gene_rows(path, genes, sep="\t"):
    """Read only requested gene rows from a gzipped genes-x-cells count matrix.
    Returns (cells, {gene: list_of_counts}). Handles both header formats:
    with a 'Gene' label column (AML) and without (MM)."""
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").split(sep)
        has_gene_col = header[0].lower() in ("gene", "symbol", "geneid", "gene_symbol", "genes")
        cells = header[1:] if has_gene_col else header[0:]
        found = {g: None for g in genes}
        for ln in f:
            parts = ln.rstrip("\n").split(sep)
            g = parts[0]
            if g in found and found[g] is None:
                found[g] = parts[1:]
            if all(v is not None for v in found.values()):
                break
    return cells, found

# ================= AML =================
def process_aml():
    aml_dir = "/home/zhh2/CLL-ScRNA-seqdatasets/AML_GSE116256/extracted"
    rows = []  # per malignant cell: ligand counts + total umi
    n_mal = 0
    for dem in sorted(glob.glob(f"{aml_dir}/*.dem.txt.gz")):
        anno = dem.replace(".dem.txt.gz", ".anno.txt.gz")
        if not os.path.exists(anno):
            continue
        ann = pd.read_csv(anno, sep="\t")
        ann = ann[ann["PredictionRefined"] == "malignant"]
        if len(ann) == 0:
            continue
        cells, gene_rows = read_gene_rows(dem, LIGANDS)
        # cell -> index mapping
        cell_idx = {c: i for i, c in enumerate(cells)}
        # total umi from anno (TranscriptomeUMIs)
        ann = ann.set_index("Cell")
        for g in LIGANDS:
            if gene_rows.get(g) is None:
                continue
            for c in ann.index:
                if c in cell_idx:
                    cnt = int(gene_rows[g][cell_idx[c]])
                    umi = ann.loc[c, "TranscriptomeUMIs"]
                    # log1p(CPM): counts/umi*1e4
                    cpm = cnt / max(umi, 1) * 1e4
                    rows.append((g, np.log1p(cpm)))
        n_mal += len(ann)
    df = pd.DataFrame(rows, columns=["gene", "expr"])
    out = df.groupby("gene")["expr"].mean()
    print(f"AML malignant cells: {n_mal}")
    return out

# ================= MM =================
def process_mm():
    mm_dir = "/home/zhh2/CLL-ScRNA-seqdatasets/MM_GSE117156/extracted"
    genes_wanted = LIGANDS + MM_MARKERS
    rows = []  # per plasma cell: ligand counts + total umi
    n_plasma = 0
    for f in sorted(glob.glob(f"{mm_dir}/*.txt.gz")):
        cells, gene_rows = read_gene_rows(f, genes_wanted)
        if "MZB1" not in gene_rows or gene_rows["MZB1"] is None:
            continue
        mzb1 = np.array([int(x) for x in gene_rows["MZB1"]])
        sdc1 = np.array([int(x) for x in gene_rows["SDC1"]]) if gene_rows.get("SDC1") else np.zeros_like(mzb1)
        cd38 = np.array([int(x) for x in gene_rows["CD38"]]) if gene_rows.get("CD38") else np.zeros_like(mzb1)
        # plasma cell: MZB1>=2 or (SDC1>=2 or CD38>=2)
        is_plasma = (mzb1 >= 2) | ((sdc1 >= 2) & (cd38 >= 2))
        # total UMI: full rowsum over the matrix (genes x cells)
        total_umi = np.zeros(len(cells), dtype=np.float64)
        with gzip.open(f, "rt") as fh:
            _header = fh.readline()  # skip
            for ln in fh:
                p = ln.rstrip("\n").split("\t")
                vals = p[1:]
                if len(vals) == len(total_umi):
                    for i, v in enumerate(vals):
                        if v and v != "0":
                            total_umi[i] += int(v)
        for g in LIGANDS:
            if gene_rows.get(g) is None:
                continue
            vals = np.array([int(x) for x in gene_rows[g]])
            for i in np.where(is_plasma)[0]:
                cpm = vals[i] / max(total_umi[i], 1) * 1e4
                rows.append((g, np.log1p(cpm)))
        n_plasma += int(is_plasma.sum())
    df = pd.DataFrame(rows, columns=["gene", "expr"])
    out = df.groupby("gene")["expr"].mean()
    print(f"MM plasma cells: {n_plasma}")
    return out

if __name__ == "__main__":
    print("=== processing AML ===")
    aml = process_aml()
    print("=== processing MM ===")
    mm = process_mm()
    res = pd.DataFrame({"AML": aml, "MM": mm}).reindex(LIGANDS)
    res.to_csv("/home/zhh2/CLL-ScRNA-seqdatasets/NK_analysis/results/cross_cancer_AML_MM.tsv", sep="\t")
    print("\n=== result (mean log1p CPM) ===")
    print(res.round(4).to_string())
