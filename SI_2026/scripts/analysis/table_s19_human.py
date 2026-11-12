"""Table S19 statistics for the human map binned as in Fig. S8, computed with the definitions of the table and checked
on the eight GenEra species of Dataset S14; also writes the human rows in the Dataset S14 format.
python table_s19_human.py <SI_2026/datasets> <figS8_human_binned.tsv> <out dir>"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

D, HUMAN, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
BINS = ["all_life", "eukaryotes", "kingdom_only", "family", "species"]


def wilson(k, n, z=1.959964):
    p = k / n
    centre = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return 100 * p, 100 * (centre - half), 100 * (centre + half)


def neff_gap(sub):
    both = sub[sub.age_bin.isin(["all_life", "species"])].dropna(subset=["log10_af2", "plddt_global"])
    if len(both) < 80:
        return None
    try:
        both = both.assign(q=pd.qcut(both.log10_af2, 4, duplicates="drop"))
    except ValueError:
        return None
    gaps = []
    for _, piece in both.groupby("q", observed=True):
        a = piece.loc[piece.age_bin.eq("all_life"), "plddt_global"]
        b = piece.loc[piece.age_bin.eq("species"), "plddt_global"]
        if len(a) >= 15 and len(b) >= 15:
            gaps.append(float(a.median() - b.median()))
    return float(np.mean(gaps)) if gaps else None


def row(sub):
    sub = sub[sub.age_bin.isin(BINS)].copy()
    sub["log10_af2"] = np.log10(sub["AF2-style_Neff"] + 1)
    n = len(sub)
    cells = {}
    for b in BINS:
        s = sub[sub.age_bin.eq(b)]
        cells[b] = f"{100 * len(s) / n:.1f}% / {s.plddt_global.median():.1f}"
    young = int(sub.age_bin.isin(["family", "species"]).sum())
    p, lo, hi = wilson(young, n)
    rho = spearmanr(sub["rank"], sub.plddt_global, nan_policy="omit").correlation
    g_all = neff_gap(sub)
    g_ord = neff_gap(sub[sub.frac_disordered < 0.2])
    fmt = lambda g: "–" if g is None else f"{g:.1f}"
    return {"n": f"{n:,}", **cells, "young": f"{p:.1f} ({lo:.1f}–{hi:.1f})", "rho": f"{rho:.2f}",
            "gap": f"{fmt(g_all)} ({fmt(g_ord)})"}


s14 = pd.read_csv(D / "Dataset_S14_published_gene_ages/published_gene_ages_afdb_v6.tsv", sep="\t")
print("check on Dataset S14 (Table S19 values: ARATH 26,914; 8.7 (8.3–9.0); −0.45; 11.3 (38.5)):")
for p in ("ARATH", "ORYSJ", "MOUSE"):
    print(" ", p, row(s14[s14.proteome.eq(p)]))

hum = pd.read_csv(HUMAN, sep="\t")
m = pd.read_parquet(D / "Dataset_S10_alignment_depth/merged_neff.parquet",
                    columns=["db_accession", "AF2-style_Neff", "tax_group7"])
hum = hum.merge(m, on="db_accession", how="left")
print("human rows:", len(hum), "| with AF2-style Neff:", int(hum["AF2-style_Neff"].notna().sum()))
print("HUMAN", row(hum))

hum = hum[s14.columns.intersection(hum.columns).tolist()]
for c in s14.columns:
    if c not in hum.columns:
        hum[c] = np.nan
hum = hum[s14.columns]
hum.to_csv(OUT / "s14_human_rows.tsv", sep="\t", index=False)
summ = []
for b in BINS:
    s = hum[hum.age_bin.eq(b)]
    summ.append({"proteome": "HUMAN", "species": "Human", "age_bin": b, "n": len(s),
                 "plddt_median": s.plddt_global.median(), "neff_median": s.neff.median(),
                 "af2style_median": s["AF2-style_Neff"].median(), "disorder_median": s.frac_disordered.median()})
pd.DataFrame(summ).to_csv(OUT / "s14_human_summary_rows.csv", index=False)
print("written", OUT / "s14_human_rows.tsv", "and s14_human_summary_rows.csv")
