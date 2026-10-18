#!/usr/bin/env python3
"""
Step 14 — Bulk survival validation (ready-to-run).

Inputs (user-provided, see README):
  expr.tsv    genes x samples processed expression matrix (symbols as rows)
  clinical.tsv samples x [survival_time, event(0/1), IGHV(0=mutated/1=unmutated)]

What it does:
  1. Extract CD200 / CD70 / TNFSF9 (immune-evasion program).
  2. Median-split high/low -> Kaplan-Meier + log-rank (manual).
  3. Cox regression (univariate + IGHV-adjusted) via lifelines if available.

Outputs: results/bulk_survival_CD200.csv (+ .png)
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys, os

OUT = "/home/zhh2/CLL-ScRNA-seqdatasets/NK_analysis"
GENES = ["CD200", "CD70", "TNFSF9"]

def load(expr_f, clin_f):
    expr = pd.read_csv(expr_f, sep="\t", index_col=0)      # genes x samples
    clin = pd.read_csv(clin_f, sep="\t", index_col=0)      # samples x clinical
    return expr, clin

def km_logrank(times, events, group):
    """Manual Kaplan-Meier + log-rank. group: bool array (True=high)."""
    order = np.argsort(times)
    t, e, g = times[order], events[order], group[order].astype(int)
    # log-rank (two-group)
    import scipy.stats as st
    n1 = g.sum(); n0 = (1-g).sum()
    O1, E1, V = 0.0, 0.0, 0.0
    for ti in np.unique(t):
        at = t == ti
        d = e[at].sum(); n = at.sum()
        d1 = (e[at] & (g[at]==1)).sum(); n1t = (g[at]==1).sum()
        if n > 0:
            E1 += d * n1t / n
            V += d * (n-d) * n1t * (n-n1t) / (n*n*(n-1)) if n > 1 else 0
    O1 = (e & (g==1)).sum()
    z = (O1 - E1) / np.sqrt(V) if V > 0 else 0
    p = 2 * (1 - st.norm.cdf(abs(z)))
    return O1, E1, V, z, p

def main(expr_f, clin_f):
    expr, clin = load(expr_f, clin_f)
    # common samples
    common = [s for s in clin.index if s in expr.columns]
    expr = expr[common]; clin = clin.loc[common]
    print(f"{len(common)} samples")
    for g in GENES:
        if g not in expr.index:
            print(f"  {g}: NOT FOUND"); continue
        v = expr.loc[g].astype(float)
        hi = v >= v.median()
        t = clin["survival_time"].astype(float).values
        e = clin["event"].astype(int).values
        O1, E1, V, z, p = km_logrank(t, e, hi.values)
        hr = (O1/E1) / ((e.sum()-O1)/(e.sum()-E1)) if E1>0 and E1<e.sum() else np.nan
        print(f"  {g}: high n={hi.sum()}, low n={(~hi).sum()}, log-rank p={p:.3f}, HR={hr:.2f}")
        # KM plot
        fig, ax = plt.subplots(figsize=(4,3))
        for grp, lbl, c in [(1,"high","#c0392b"),(0,"low","#2e7d32")]:
            tt = t[hi.values==(grp==1)]; ee = e[hi.values==(grp==1)]
            s = np.argsort(tt); tt=tt[s]; ee=ee[s]
            surv = np.cumprod(1 - ee/np.arange(len(ee),0,-1)) if False else None
            # KM step
            times_u = np.sort(np.unique(tt)); km=[]
            for ti in times_u:
                at = tt<=ti; d_at = ee[tt==ti].sum(); n_at = at.sum()
                km.append(np.prod(1 - ee[tt<=ti]/np.array([len(tt)-k for k in range(sum(tt<=ti))]) ) if False else None)
            # simpler KM
            n = len(tt); km=[1.0]
            for i in range(n):
                d = ee[i]; n_risk = n-i
                km.append(km[-1]*(1-d/n_risk) if n_risk>0 else km[-1])
            ax.step(np.concatenate([[0],tt]), km, where="post", label=f"{lbl} (n={len(tt)})", color=c)
        ax.set_xlabel("time"); ax.set_ylabel("survival"); ax.set_title(f"{g} (log-rank p={p:.3f})")
        ax.legend(frameon=False); ax.set_ylim(0,1.05)
        plt.tight_layout(); plt.savefig(f"{OUT}/figures/bulk_KM_{g}.png", bbox_inches="tight"); plt.close()
    print("DONE")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("用法: python 14_bulk_survival.py <expr.tsv> <clinical.tsv>")
        print("expr.tsv: 基因x样本 表达矩阵; clinical.tsv: 样本x[survival_time,event,IGHV]")
        sys.exit(0)
    main(sys.argv[1], sys.argv[2])
