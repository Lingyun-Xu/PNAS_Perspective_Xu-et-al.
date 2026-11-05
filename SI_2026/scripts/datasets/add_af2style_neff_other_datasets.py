"""Add AF2-style Neff to the other datasets that carry Neff (Datasets S2, S3, S13, S14 and S17, FigS20), from Dataset S10
by UniProt accession, next to the Neff column; add its median and quartiles to the summary tables of Datasets S13 and
S14. Before anything is written, the Neff of each table must equal the Neff of Dataset S10 for the same accession, and
the Neff columns of the summary tables must be reproduced from the per-protein table. All other columns stay as they are.
FigS19/neff_sets.tsv is written by plot_neff_sets.py.
Usage: python add_af2style_neff_other_datasets.py <SI_2026/datasets>"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

D = Path(sys.argv[1])
COL = "AF2-style_Neff"
s10 = pd.read_parquet(D / "Dataset_S10_alignment_depth/merged_neff.parquet", columns=["db_accession", "neff", COL])
assert s10.db_accession.is_unique
NEFF = s10.set_index("db_accession").neff
AF2 = s10.set_index("db_accession")[COL]


def with_column(t, acc, neff_col, name):
    """the table with AF2-style Neff right after its Neff column, once its Neff is found equal to Dataset S10"""
    assert COL not in t.columns, name
    ours, s = t[neff_col], t[acc].map(NEFF)
    both = ours.notna() & s.notna()
    diff = (ours[both] - s[both]).abs().max()
    only_here = int((ours.notna() & s.isna()).sum())
    assert diff < 1e-3 and only_here == 0, (name, diff, only_here)
    out = t.copy()
    out.insert(t.columns.get_loc(neff_col) + 1, COL, t[acc].map(AF2))
    print(f"{name}: {len(t):,} rows; Neff {int(ours.notna().sum()):,}, equal to Dataset S10 (largest difference {diff:.6f}); "
          f"AF2-style Neff {int(out[COL].notna().sum()):,}", flush=True)
    return out


# ---------------------------------------------------------------- Dataset S13
p = D / "Dataset_S13_homolog_class/panel_strata.parquet"
s13 = with_column(pd.read_parquet(p), "db_accession", "neff", "S13 panel_strata")
p2 = D / "Dataset_S13_homolog_class/panel_bin_summary.csv"
b13 = pd.read_csv(p2)
g = s13.groupby(["tax_group7", "panel_bin"])
sm = pd.DataFrame({"neff_median": g.neff.median(), "af2style_median": g[COL].median(),
                   "af2style_q25": g[COL].quantile(.25), "af2style_q75": g[COL].quantile(.75)}).reset_index()
m = b13.merge(sm, left_on=["group", "bin"], right_on=["tax_group7", "panel_bin"], how="left", suffixes=("", "_new"))
assert m.loc[m.n > 0, "tax_group7"].notna().all()                       # the summary keeps empty classes (n = 0)
assert (m.neff_median - m.neff_median_new).abs().max() < 0.06, (m.neff_median - m.neff_median_new).abs().max()
b13new = m.drop(columns=["tax_group7", "panel_bin", "neff_median_new"])
for c in ("af2style_median", "af2style_q25", "af2style_q75"):
    b13new[c] = b13new[c].round(1)
print("S13 summary: Neff medians reproduced;", len(b13new), "rows", flush=True)

# ---------------------------------------------------------------- Dataset S14
p3 = D / "Dataset_S14_published_gene_ages/published_gene_ages_afdb_v6.tsv"
s14 = with_column(pd.read_csv(p3, sep="\t", low_memory=False), "db_accession", "neff", "S14 published ages")
p4 = D / "Dataset_S14_published_gene_ages/suite_bin_summary.csv"
b14 = pd.read_csv(p4)
g = s14.groupby(["proteome", "age_bin"])
sm = pd.DataFrame({"neff_median_new": g.neff.median(), "af2style_median": g[COL].median()}).reset_index()
m = b14.merge(sm, on=["proteome", "age_bin"], how="left")
assert m.af2style_median.notna().all()
d = (m.neff_median - m.neff_median_new).abs().max()
print("S14 summary: largest difference of the Neff medians", d, flush=True)
assert d < 0.06
b14new = m.drop(columns=["neff_median_new"])
b14new.insert(b14new.columns.get_loc("neff_median") + 1, "af2style_median", b14new.pop("af2style_median").round(1))
p5 = D / "Dataset_S14_published_gene_ages/genera_balanced_sample_ages.csv"
s14b = with_column(pd.read_csv(p5), "gene", "neff", "S14 balanced sample")
p6 = D / "Dataset_S14_published_gene_ages/genera_balanced_sample_summary.csv"
b14b = pd.read_csv(p6)
g = s14b.groupby(["tax_group7", "age_bin"])
sm = pd.DataFrame({"neff_median_new": g.neff.median(), "af2style_median": g[COL].median()}).reset_index()
m = b14b.merge(sm, left_on=["group", "bin"], right_on=["tax_group7", "age_bin"], how="left")
assert m.loc[m.n > 0, "tax_group7"].notna().all()
d = (m.neff_median - m.neff_median_new).abs().max()
print("S14 balanced summary: largest difference of the Neff medians", d, flush=True)
assert d < 0.06
b14bnew = m.drop(columns=["tax_group7", "age_bin", "neff_median_new"])
b14bnew["af2style_median"] = b14bnew.af2style_median.round(1)

# ---------------------------------------------------------------- Dataset S17, FigS20
frames = {}
for tag in ("sel3800", "sel690"):
    q = D / f"Dataset_S17_source_data_SI_figures/FigS20/{tag}_genera_ages.csv"
    frames[q] = with_column(pd.read_csv(q), "db_accession", "neff", f"S17 FigS20 {tag}")

# ---------------------------------------------------------------- Datasets S2, S3 (xlsx)
xl = {}
for name in ("Dataset_S2_sel3800.xlsx", "Dataset_S3_sel690.xlsx"):
    q = D / name
    data = pd.read_excel(q, sheet_name="data")
    notes = pd.read_excel(q, sheet_name="README")
    t = with_column(data, "uniprot_accession", "neff_uniclust30", name)
    old = "Neff (OpenProteinSet UniClust30);"
    assert notes.README.str.contains(old, regex=False).sum() == 1
    notes["README"] = notes.README.str.replace(old, "Neff and AF2-style Neff (OpenProteinSet UniClust30 alignments; Section S4.1);", regex=False)
    xl[q] = (notes, t)

# ---------------------------------------------------------------- all checks passed: write
s13.to_parquet(p, index=False)
b13new.to_csv(p2, index=False)
s14.to_csv(p3, sep="\t", index=False)
b14new.to_csv(p4, index=False)
s14b.to_csv(p5, index=False)
b14bnew.to_csv(p6, index=False)
for q, t in frames.items():
    t.to_csv(q, index=False)
for q, (notes, t) in xl.items():
    with pd.ExcelWriter(q, engine="openpyxl") as xw:
        notes.to_excel(xw, sheet_name="README", index=False)
        t.to_excel(xw, sheet_name="data", index=False)
print("written", flush=True)
for f in [p, p2, p3, p4, p5, p6, *frames, *xl]:
    print(f"  {f.name}: {f.stat().st_size / 2 ** 20:.1f} MiB")
