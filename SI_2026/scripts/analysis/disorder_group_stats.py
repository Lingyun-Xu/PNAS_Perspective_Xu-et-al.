"""Disorder by group and plant proteome for the concise PNAS SI table on disorder (canonical entries of the 47 AFDB v6
proteomes, seven groups as in the other group tables).
  descriptive: n, median and IQR of the AIUPred disordered fraction, share of proteins with a disordered region of
               >= 30 residues, median global pLDDT
  disorder between groups: Kruskal-Wallis (epsilon^2); plants versus each group, Cliff's delta from Mann-Whitney U
               (Holm); among the four plant proteomes, Kruskal-Wallis and Dunn-Holm with Cliff's delta
  disorder and pLDDT: Spearman rho with 95% CI, plants versus each group (Fisher z, Holm), heterogeneity among the
               plant proteomes (Cochran Q); slope of global pLDDT on disordered fraction (least squares, per 10
               percentage points of disordered residues, 95% CI from heteroscedasticity-robust HC1 standard errors)
The tests use stats_vm.py (the functions of SI_2026/statistics); their results are checked against the installed
group_tests.tsv, among_plant_dunn.tsv and correlation_tests.tsv.
Usage: python disorder_group_stats.py <out_dir>"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "/path/to/project/template_esm/campaign/stats")
import stats_vm as sv  # noqa: E402

D = Path("/path/to/project/SI_2026/datasets")
ST = Path("/path/to/project/SI_2026/statistics")
out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)
GROUPS = ["Viridiplantae", "Archaea", "Bacteria", "Protist", "Fungi", "Animalia", "Parasitic worms"]
PL = {"Arabidopsis thaliana": "Arabidopsis", "Glycine max": "soybean", "Oryza sativa subsp. japonica": "rice",
      "Zea mays": "maize"}

dis = pd.read_parquet(D / "Dataset_S12_disorder/merged_per_protein.parquet",
                      columns=["db_accession", "organism_scientific", "plddt_global", "frac_disordered", "n_idr_30"])
dis = dis[~dis.db_accession.duplicated(keep=False)]                       # canonical: one model per accession
td = pd.read_parquet(D / "Dataset_S17_source_data_SI_figures/FigS06_S07/template_distances.parquet",
                     columns=["db_accession", "tax_group7"])
df = dis.merge(td, on="db_accession", how="inner")
df["tax_group7"] = sv.norm_group(df.tax_group7.astype(str))
df["plant"] = df.organism_scientific.where(df.organism_scientific.isin(list(sv.PLANTS)))
print("proteins:", len(df), df.tax_group7.value_counts().reindex(GROUPS).to_dict(), flush=True)


def slope(d):
    x, y = d.frac_disordered.to_numpy(float), d.plddt_global.to_numpy(float)
    xc = x - x.mean()
    b = (xc * (y - y.mean())).sum() / (xc ** 2).sum()
    e = y - y.mean() - b * xc
    se = np.sqrt((xc ** 2 * e ** 2).sum() / (xc ** 2).sum() ** 2 * len(x) / (len(x) - 2))     # HC1
    return 0.1 * b, 0.1 * (b - 1.96 * se), 0.1 * (b + 1.96 * se)


rows = []
subsets = [(g, df[df.tax_group7 == g]) for g in GROUPS] + [(PL[p], df[df.plant == p]) for p in PL]
for name, d in subsets:
    q1, q2, q3 = d.frac_disordered.quantile([0.25, 0.5, 0.75])
    s, lo, hi = slope(d)
    rows.append({"subset": name, "n": len(d), "median disordered fraction": q2, "q25": q1, "q75": q3,
                 "share with disordered region >= 30 residues": float((d.n_idr_30 > 0).mean()),
                 "median global pLDDT": d.plddt_global.median(),
                 "slope per 10% disordered residues": s, "slope ci_low": lo, "slope ci_high": hi})
desc = pd.DataFrame(rows)

# tests, with the functions of the installed statistics
g, a = sv.compare(df, "frac_disordered", "tax_group7", "Disordered fraction (AIUPred)")
c = sv.correlations(df, "frac_disordered", "plddt_global", "tax_group7", "pLDDT vs disordered fraction")
for name, new in (("group_tests.tsv", g), ("among_plant_dunn.tsv", a), ("correlation_tests.tsv", c)):
    old = pd.read_csv(ST / name, sep="\t")
    old = old[old.measure == new.measure.iloc[0]].reset_index(drop=True)
    num = [k for k in old.columns if k in new.columns and pd.api.types.is_numeric_dtype(old[k])]
    same = np.allclose(old[num].to_numpy(float), new.reset_index(drop=True)[num].to_numpy(float), equal_nan=True,
                       rtol=1e-9, atol=1e-12)
    print(f"{name}: reproduces the installed rows = {same}", flush=True)
desc.to_csv(out / "disorder_by_group_descriptive.tsv", sep="\t", index=False)
g.to_csv(out / "disorder_group_tests.tsv", sep="\t", index=False)
a.to_csv(out / "disorder_among_plant_dunn.tsv", sep="\t", index=False)
c.to_csv(out / "disorder_pLDDT_correlation_tests.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
print(desc.round(3).to_string(index=False))
print(c[["subset", "n", "rho", "ci_low", "ci_high", "p", "statistic", "p_holm"]].to_string(index=False))
print("DONE")
