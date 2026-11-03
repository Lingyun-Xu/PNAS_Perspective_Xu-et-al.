"""Checks for the expanded GenEra text: mouse species-level genes in the panel classes; Spearman p-values of the
published maps; sel3800/sel690 composition by crop species; installed gene-age test statistics."""
import numpy as np
import pandas as pd
from scipy import stats

D = "/path/to/project/SI_2026/datasets/"
s = pd.read_csv(D + "Dataset_S14_published_gene_ages/published_gene_ages_afdb_v6.tsv", sep="\t")
panel = pd.read_parquet(D + "Dataset_S13_homolog_class/panel_strata.parquet", columns=["db_accession", "panel_bin"])
m = s.merge(panel, on="db_accession", how="left")
print("mouse species-level genes by panel class:",
      m[(m.proteome == "MOUSE") & (m.age_bin == "species")].panel_bin.value_counts().to_dict())
BINS = ["all_life", "eukaryotes", "kingdom_only", "family", "species"]
for p, sub in s[s.age_bin.isin(BINS)].groupby("proteome", sort=False):
    r = stats.spearmanr(sub["rank"], sub.plddt_global, nan_policy="omit")
    print(p, "rho", round(r.statistic, 3), "p", r.pvalue, "n", int(sub[["rank", "plddt_global"]].dropna().shape[0]))

F = "/path/to/project/template_esm/campaign/figures/"
for tag, f in (("sel3800", F + "genera_ages_3800/sel3800_genera_ages.csv"), ("sel690", F + "genera_ages_690/sel690_genera_ages.csv")):
    g = pd.read_csv(f)
    aged = g[g.age_bin.isin(["all_life", "eukaryotes", "kingdom_only", "family_or_genus"])]
    print(tag, "rows", len(g), "aged", len(aged), "| species col:", [c for c in g.columns if "species" in c or "organism" in c][:3])
    sp = [c for c in g.columns if c in ("species", "organism_scientific", "organism")][0]
    for name in ["Arabidopsis thaliana", "Oryza sativa subsp. japonica", "Zea mays", "Glycine max"]:
        a = aged[aged[sp] == name]
        if len(a):
            print(f"  {name}: n={len(a)} family/genus {100 * (a.age_bin == 'family_or_genus').mean():.0f}% "
                  f"all life {100 * (a.age_bin == 'all_life').mean():.0f}%")
gt = pd.read_csv("/path/to/project/SI_2026/statistics/genera_tests.tsv", sep="\t")
pd.set_option("display.width", 250)
print(gt.to_string(index=False))
