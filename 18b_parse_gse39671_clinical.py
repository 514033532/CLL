#!/usr/bin/env python3
"""
Parse GSE39671 series matrix (possibly truncated gz) -> clinical table.
Extracts: GSM accession, time-to-first-treatment (years), treatment event.
"""
import gzip
import re
import csv

SRC = "/home/zhh2/CLL-ScRNA-seqdatasets/GSE39671_RAW/GSE39671_series_matrix.txt.gz"
OUT = "/home/zhh2/CLL-ScRNA-seqdatasets/NK_analysis/data/GSE39671_clinical.tsv"

gsm = []
ttft = []
event = []
char_lines = []  # all characteristics lines, in order

# read lines, tolerate truncation of the gz stream
lines = []
try:
    with gzip.open(SRC, "rt", errors="replace") as fh:
        for ln in fh:
            lines.append(ln.rstrip("\n"))
except (EOFError, OSError):
    pass  # truncated file: header metadata is complete

def fields(ln):
    # split on tab, strip quotes
    return [x.strip('"') for x in ln.split("\t")]

for ln in lines:
    if ln.startswith("!Sample_geo_accession"):
        gsm = fields(ln)[1:]
    elif ln.startswith("!Sample_characteristics_ch1"):
        vals = fields(ln)[1:]
        if vals and vals[0].startswith("sampling time to first treatment"):
            ttft = vals
        elif vals and vals[0].startswith("treatment event"):
            event = vals

print(f"gsm n={len(gsm)}, ttft n={len(ttft)}, event n={len(event)}")

rows = []
for i in range(len(gsm)):
    t = None
    if i < len(ttft):
        m = re.search(r"([0-9.]+)", ttft[i])
        t = float(m.group(1)) if m else None
    e = None
    if i < len(event):
        m = re.search(r"([01])$", event[i])
        e = int(m.group(1)) if m else None
    rows.append((gsm[i], t, e))

with open(OUT, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["gsm", "ttft_years", "event"])
    for r in rows:
        w.writerow(r)

print(f"wrote {len(rows)} rows to {OUT}")
print("first 5:", rows[:5])
print("last 3:", rows[-3:])
print("events=1:", sum(1 for r in rows if r[2] == 1), "events=0:", sum(1 for r in rows if r[2] == 0))
print("ttft missing:", sum(1 for r in rows if r[1] is None))
