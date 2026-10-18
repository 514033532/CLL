#!/usr/bin/env python3
"""
Step 18e — GSE39671 survival with positive-control QC.
Core: TNFSF9/CD86/CD200/CD70. Controls: ZAP70, LPL, CD38 (established CLL markers).
If controls replicate (high -> shorter TTFT, p<0.05) then pipeline is sound and a
non-replication of the core genes is credible.
"""
import pandas as pd
import numpy as np
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test
import scipy.stats as st

OUT = "/home/zhh2/CLL-ScRNA-seqdatasets/NK_analysis"
expr = pd.read_csv(f"{OUT}/data/GSE39671_expr.tsv", sep="\t")
clin = pd.read_csv(f"{OUT}/data/GSE39671_clinical.tsv", sep="\t")
m = clin.merge(expr, on="gsm", how="inner")
t = m["ttft_years"].astype(float); e = m["event"].astype(int)

GENES = ["CD200", "CD70", "TNFSF9", "CD86", "ZAP70", "LPL", "CD38"]

print(f"{'gene':8s} {'role':10s} {'logrank_p':>10s} {'cox_HR':>7s} {'cox_p':>8s}  medTTFT(hi/lo)")
rows = []
for g in GENES:
    z = (m[g] - m[g].mean()) / m[g].std(ddof=0)
    hi = z >= z.median()
    lr = logrank_test(t[hi], t[~hi], e[hi], e[~hi])
    cph = CoxPHFitter()
    d = m[["ttft_years", "event"]].copy(); d["x"] = z.values
    cph.fit(d, "ttft_years", "event")
    hr = float(cph.hazard_ratios_["x"]); p = float(cph.summary.loc["x", "p"])
    role = "control" if g in ("ZAP70", "LPL", "CD38") else "core"
    print(f"{g:8s} {role:10s} {lr.p_value:10.4f} {hr:7.2f} {p:8.4f}  {t[hi].median():.2f}/{t[~hi].median():.2f}")

# --- ZAP70 individual probe 214032_at (canonical) from full matrix ---
print("\n--- ZAP70 individual probe check ---")
full = pd.read_csv(f"{OUT}/data/GSE39671_expr_full_probes.tsv", sep="\t", index_col=0)
if "214032_at" in full.index:
    zap = full.loc["214032_at"]
    zm = clin.merge(zap.rename("ZAP70_214032").reset_index().rename(columns={"index": "gsm"}), on="gsm")
    z = (zm["ZAP70_214032"] - zm["ZAP70_214032"].mean()) / zm["ZAP70_214032"].std(ddof=0)
    for split_name, hi in [("median", z >= z.median()),
                            ("top30%", z >= z.quantile(0.70))]:
        lr = logrank_test(t[hi], t[~hi], e[hi], e[~hi])
        print(f"ZAP70 214032_at {split_name}-split: logrank_p={lr.p_value:.4f}, "
              f"medTTFT hi/lo = {t[hi].median():.2f}/{t[~hi].median():.2f}")
    # spearman ZAP70 vs LPL (both should be IGHV-surrogates, correlated)
    lpl = m["LPL"]
    rho, p = st.spearmanr(zm["ZAP70_214032"], lpl)
    print(f"Spearman ZAP70~LPL: rho={rho:.3f}, p={p:.2e} (should be positive if both mark unmutated IGHV)")
