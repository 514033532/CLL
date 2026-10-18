#!/usr/bin/env python3
"""Inspect AML (GSE116256) .dem/.anno and MM (GSE117156) .txt formats."""
import tarfile, gzip, os, glob

def peek_anno(path):
    with gzip.open(path, "rt") as f:
        for i, ln in enumerate(f):
            if i < 3:
                print("    ", ln.rstrip())
            else:
                break

def peek_matrix(path, n_lines=3):
    with gzip.open(path, "rt") as f:
        for i, ln in enumerate(f):
            if i < n_lines:
                print("    ", ln.rstrip()[:200])
            else:
                break

# ---- AML ----
aml_dir = "/home/zhh2/CLL-ScRNA-seqdatasets/AML_GSE116256"
tar = f"{aml_dir}/GSE116256_RAW.tar"
if os.path.exists(tar):
    print("=== AML tar members (first 12) ===")
    with tarfile.open(tar) as t:
        names = t.getnames()
        print("  total members:", len(names))
        for n in names[:12]:
            print("  ", n)
        # extract first sample's .anno and .dem
        anno = [n for n in names if n.endswith(".anno.txt.gz")][:1]
        dem = [n for n in names if n.endswith(".dem.txt.gz")][:1]
        if anno:
            t.extract(anno[0], path=aml_dir + "/extracted")
            print("=== AML anno sample:", os.path.basename(anno[0]), "===")
            peek_anno(f"{aml_dir}/extracted/{anno[0]}")
        if dem:
            t.extract(dem[0], path=aml_dir + "/extracted")
            print("=== AML dem sample:", os.path.basename(dem[0]), "===")
            peek_matrix(f"{aml_dir}/extracted/{dem[0]}")
else:
    print("AML tar not yet downloaded")
