"""Recompute the gene-age numbers of comprehensive SI Section S{B6}.2-3 from the datasets before they are used in the
concise PNAS SI: per-species class counts and medians, Spearman rho (rank vs pLDDT), the matched-Neff gaps (algorithm
of plot_published_suite.neff_gap), the agreement with the panel classes, and the balanced-sample summary."""
import numpy as np
import pandas as pd
from scipy import stats

D = "/path/to/project/SI_2026/datasets/"
s = pd.read_csv(D + "Dataset_S14_published_gene_ages/published_gene_ages_afdb_v6.tsv", sep="\t")
print("rows", len(s), "| sources", s.source.value_counts().to_dict(), "| species", s.proteome.nunique())
print("age_bin values", s.age_bin.value_counts(dropna=False).to_dict())
BINS = ["all_life", "eukaryotes", "kingdom_only", "family", "species"]
s["log10_neff"] = np.log10(s.neff.clip(lower=0) + 1)


def neff_gap(sub):
    both = sub[sub.age_bin.isin(["all_life", "species"])].dropna(subset=["log10_neff", "plddt_global"]).copy()
    if len(both) < 80:
        return None
    both["q"] = pd.qcut(both.log10_neff, 4, duplicates="drop")
    gaps = []
    for _, piece in both.groupby("q", observed=True):
        a = piece.loc[piece.age_bin.eq("all_life"), "plddt_global"]
        b = piece.loc[piece.age_bin.eq("species"), "plddt_global"]
        if len(a) >= 15 and len(b) >= 15:
            gaps.append(float(a.median() - b.median()))
    return round(float(np.mean(gaps)), 1) if gaps else None


rows = []
for p, sub in s.groupby("proteome", sort=False):
    use = sub[sub.age_bin.isin(BINS)]
    r = {"proteome": p, "matched": sub.db_accession.nunique(), "binned": len(use),
         "excluded (flag/TR<=30)": int((~sub.age_bin.isin(BINS)).sum())}
    for b in BINS:
        x = use[use.age_bin == b]
        r[b] = f"{len(x)}/{x.plddt_global.median():.1f}"
    r["pct species"] = round(100 * (use.age_bin == "species").mean(), 1)
    r["pct family"] = round(100 * (use.age_bin == "family").mean(), 1)
    rk = use.dropna(subset=["rank", "plddt_global"])
    r["rho rank"] = round(stats.spearmanr(rk["rank"], rk.plddt_global).statistic, 2)
    r["rho bin"] = round(stats.spearmanr(rk.age_bin.map({b: i for i, b in enumerate(BINS)}), rk.plddt_global).statistic, 2)
    r["gap matched Neff"] = neff_gap(use)
    r["gap ordered"] = neff_gap(use[use.frac_disordered < 0.2])
    sp = use[use.age_bin == "species"]
    r["species-level median Neff / disorder"] = f"{sp.neff.median():.1f} / {sp.frac_disordered.median():.2f}"
    rows.append(r)
pd.set_option("display.width", 260)
print(pd.DataFrame(rows).to_string(index=False))

panel = pd.read_parquet(D + "Dataset_S13_homolog_class/panel_strata.parquet")
print("panel columns:", list(panel.columns)[:12])
key = "db_accession" if "db_accession" in panel.columns else panel.columns[0]
bincol = [c for c in panel.columns if "bin" in c or "class" in c or "stratum" in c]
print("panel class column candidates:", bincol)
m = s[s.age_bin.isin(BINS)].merge(panel[[key] + bincol[:1]], left_on="db_accession", right_on=key, how="left")
pc = bincol[0]
for p, sub in m.groupby("proteome", sort=False):
    sp = sub[sub.age_bin == "species"]
    al = sub[sub.age_bin == "all_life"]
    print(p, "| species-level that are species-only in panel:", f"{(sp[pc] == 'species_only').mean():.3f}", f"(n={len(sp)})",
          "| all-life that are universal in panel:", f"{(al[pc] == 'universal').mean():.3f}")

g = pd.read_csv("/path/to/project/conservation/figures/sample_genera_ages.csv")
print("sample rows", len(g), "| proteomes", g.proteome.nunique(), "| age_bin", g.age_bin.value_counts(dropna=False).to_dict())
aged = g[g.age_bin.isin(["all_life", "eukaryotes", "kingdom_only", "family_or_genus"])]
print("aged", len(aged), aged.groupby("age_bin").plddt_global.agg(["size", "median"]).round(1).to_dict())
grp = aged.tax_group7.str.replace("Parasitic Worms", "Parasitic worms")
tab = pd.crosstab(grp.eq("Viridiplantae"), aged.age_bin)
chi2, pval, dof, _ = stats.chi2_contingency(tab)
print("plants vs others: chi2", round(chi2, 1), "df", dof, "p", round(pval, 4), "V", round(np.sqrt(chi2 / tab.values.sum() / (min(tab.shape) - 1)), 3))
print("per proteome counts (min/max):", g.groupby("proteome").size().min(), g.groupby("proteome").size().max())
