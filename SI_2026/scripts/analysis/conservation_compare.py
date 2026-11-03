"""Plants versus other groups and species in conservation and gene age (concise PNAS SI, Section S4.4).
A. Homolog classes within the 47 proteomes (Dataset S13): class composition per group and plant proteome, share
   restricted to the own kingdom, number of other same-kingdom proteomes in the panel, chi-square tests (plants vs
   each eukaryotic group; among plant proteomes), species-only share of every eukaryotic proteome.
B. Published GenEra maps (Dataset S14): young-gene shares with 95% Wilson intervals, chi-square tests (pooled plants vs
   pooled animals / fungi; each plant species vs each animal and fungal species), pLDDT within each age class (plants
   vs animals, plants vs fungi: Mann-Whitney U, Cliff's delta, Holm), and the plant mean pLDDT reweighted to the
   animal and fungal age composition.
C. Balanced GenEra sample (Dataset S14): plants vs each group, chi-square on the four age bins, young share with CI.
Usage: python conservation_compare.py <out_dir>"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

D = Path("/path/to/project/SI_2026/datasets")
out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)
pd.set_option("display.width", 250)


def chi2(tab):
    tab = tab.loc[:, tab.sum(axis=0) > 0]
    c, p, dof, _ = stats.chi2_contingency(tab)
    return c, dof, p, float(np.sqrt(c / tab.values.sum() / (min(tab.shape) - 1)))


def wilson(k, n, z=1.96):
    if n == 0:
        return np.nan, np.nan
    p = k / n
    den = 1 + z ** 2 / n
    c = (p + z ** 2 / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z ** 2 / (4 * n ** 2)) / den
    return c - h, c + h


def cliffs(a, b):
    u = stats.mannwhitneyu(a, b, alternative="two-sided")
    return 2 * u.statistic / (len(a) * len(b)) - 1, u.pvalue


def holm(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    adj = np.maximum.accumulate((len(p) - np.arange(len(p))) * p[o])
    res = np.empty_like(adj)
    res[o] = np.minimum(adj, 1)
    return res


# ---------------------------------------------------------------- A. homolog classes (panel)
P = pd.read_parquet(D / "Dataset_S13_homolog_class/panel_strata.parquet")
P["group"] = P.tax_group7.replace({"Parasitic Worms": "Parasitic worms"})
CL = ["universal", "eukaryote_wide", "kingdom_only", "family_or_genus", "species_only"]
PL = {"Arabidopsis thaliana": "Arabidopsis", "Glycine max": "soybean", "Oryza sativa subsp. japonica": "rice",
      "Zea mays": "maize"}
P["plant"] = P.organism_scientific.map(PL)
nprot = P.groupby("group").proteome.nunique()
same_kingdom = {"Viridiplantae": nprot["Viridiplantae"] - 1, "Animalia": nprot["Animalia"] + nprot["Parasitic worms"] - 1,
                "Parasitic worms": nprot["Animalia"] + nprot["Parasitic worms"] - 1, "Fungi": nprot["Fungi"] - 1}
rows = []
subsets = ([("Viridiplantae", P[P.group == "Viridiplantae"])] + [(n, P[P.plant == n]) for n in PL.values()]
           + [(g, P[P.group == g]) for g in ("Animalia", "Parasitic worms", "Fungi")])
for name, d in subsets:
    r = {"subset": name, "n": len(d), "other same-kingdom proteomes": same_kingdom.get(name, 3)}
    for c in CL:
        x = d[d.panel_bin == c]
        r[f"{c} n"], r[f"{c} %"] = len(x), 100 * len(x) / len(d)
        r[f"{c} median pLDDT"] = x.plddt_global.median() if len(x) else np.nan
    r["own kingdom only %"] = 100 * d.panel_bin.isin(["kingdom_only", "family_or_genus", "species_only"]).mean()
    rows.append(r)
A = pd.DataFrame(rows)
A.to_csv(out / "panel_class_composition.tsv", sep="\t", index=False)
print(A[["subset", "n", "other same-kingdom proteomes"] + [f"{c} %" for c in CL] + ["own kingdom only %"]].round(1).to_string(index=False))

tests = []
plants = P[P.group == "Viridiplantae"]
for g in ("Animalia", "Parasitic worms", "Fungi"):
    tab = pd.crosstab(P[P.group.isin(["Viridiplantae", g])].group, P[P.group.isin(["Viridiplantae", g])].panel_bin)[CL]
    c, dof, p, v = chi2(tab)
    o = P[P.group == g]
    tests.append({"data": "homolog classes (47 proteomes)", "comparison": f"plants vs {g}", "chi2": c, "df": dof, "p": p,
                  "V": v, "plants species-only %": 100 * (plants.panel_bin == "species_only").mean(),
                  "other species-only %": 100 * (o.panel_bin == "species_only").mean(),
                  "plants own-kingdom %": 100 * plants.panel_bin.isin(CL[2:]).mean(),
                  "other own-kingdom %": 100 * o.panel_bin.isin(CL[2:]).mean()})
    # the same, with the classes restricted to the own kingdom pooled (independent of the number of relatives)
    tab2 = pd.crosstab(P[P.group.isin(["Viridiplantae", g])].group,
                       P[P.group.isin(["Viridiplantae", g])].panel_bin.isin(CL[2:]))
    c2, dof2, p2, v2 = chi2(tab2)
    tests.append({"data": "homolog classes, own kingdom only (pooled) vs beyond", "comparison": f"plants vs {g}",
                  "chi2": c2, "df": dof2, "p": p2, "V": v2})
tab = pd.crosstab(plants.plant, plants.panel_bin)[CL]
tests.append(dict(zip(["data", "comparison", "chi2", "df", "p", "V"],
                      ["homolog classes (47 proteomes)", "among the four plant proteomes", *chi2(tab)])))
per_prot = (P[P.group.isin(["Viridiplantae", "Animalia", "Parasitic worms", "Fungi"])]
            .groupby(["group", "organism_scientific"]).panel_bin.apply(lambda s: 100 * (s == "species_only").mean())
            .rename("species-only %").reset_index())
per_prot.to_csv(out / "panel_species_only_per_proteome.tsv", sep="\t", index=False)
print(per_prot.groupby("group")["species-only %"].agg(["min", "median", "max", "size"]).round(1).to_string())

# ---------------------------------------------------------------- B. published GenEra maps (nr)
S = pd.read_csv(D / "Dataset_S14_published_gene_ages/published_gene_ages_afdb_v6.tsv", sep="\t")
BINS = ["all_life", "eukaryotes", "kingdom_only", "family", "species"]
S = S[S.age_bin.isin(BINS)].copy()
KIND = {"ARATH": "plant", "ORYSJ": "plant", "SOYBN": "plant", "MOUSE": "animal", "DROME": "animal", "CAEEL": "animal",
        "YEAST": "fungus", "SCHPO": "fungus"}
S["kind"] = S.proteome.map(KIND)
rows = []
for p, d in S.groupby("proteome", sort=False):
    n = len(d)
    k_sp, k_young = int((d.age_bin == "species").sum()), int(d.age_bin.isin(["family", "species"]).sum())
    k_own = int(d.age_bin.isin(["kingdom_only", "family", "species"]).sum())
    rows.append({"species": p, "n": n, "species-level %": 100 * k_sp / n,
                 "species-level CI": tuple(100 * x for x in wilson(k_sp, n)),
                 "young %": 100 * k_young / n, "young CI": tuple(100 * x for x in wilson(k_young, n)),
                 "own kingdom only %": 100 * k_own / n})
B = pd.DataFrame(rows)
B.to_csv(out / "published_maps_young_shares.tsv", sep="\t", index=False)
print(B.round(2).to_string(index=False))
for a, b in (("plant", "animal"), ("plant", "fungus")):
    d = S[S.kind.isin([a, b])]
    c, dof, p, v = chi2(pd.crosstab(d.kind, d.age_bin)[BINS])
    tests.append({"data": "published GenEra maps (nr)", "comparison": f"plants (3 species) vs {b}s (pooled)",
                  "chi2": c, "df": dof, "p": p, "V": v})
pair = []
for ps in ("ARATH", "ORYSJ", "SOYBN"):
    for os_ in ("MOUSE", "DROME", "CAEEL", "YEAST", "SCHPO"):
        d = S[S.proteome.isin([ps, os_])]
        c, dof, p, v = chi2(pd.crosstab(d.proteome, d.age_bin)[BINS])
        pair.append({"plant": ps, "other": os_, "chi2": c, "df": dof, "p": p, "V": v,
                     "young % plant": 100 * S[S.proteome == ps].age_bin.isin(["family", "species"]).mean(),
                     "young % other": 100 * S[S.proteome == os_].age_bin.isin(["family", "species"]).mean()})
pair = pd.DataFrame(pair)
pair["p_holm"] = holm(pair.p)
pair.to_csv(out / "published_maps_pairwise_composition.tsv", sep="\t", index=False)
print(pair.round(3).to_string(index=False))
for ps in ("ARATH", "ORYSJ", "SOYBN"):
    d = S[S.proteome.isin(["ARATH", "ORYSJ", "SOYBN"])]
c, dof, p, v = chi2(pd.crosstab(S[S.kind == "plant"].proteome, S[S.kind == "plant"].age_bin)[BINS])
tests.append({"data": "published GenEra maps (nr)", "comparison": "among the three plant species", "chi2": c, "df": dof,
              "p": p, "V": v})

wc = []
for b in BINS:
    pl = S[(S.kind == "plant") & (S.age_bin == b)].plddt_global.dropna()
    for other in ("animal", "fungus"):
        ot = S[(S.kind == other) & (S.age_bin == b)].plddt_global.dropna()
        if len(pl) >= 20 and len(ot) >= 20:
            dlt, p = cliffs(pl.to_numpy(), ot.to_numpy())
            wc.append({"age class": b, "comparison": f"plants vs {other}s", "n plants": len(pl), "n other": len(ot),
                       "median plants": pl.median(), "median other": ot.median(), "delta": dlt, "p": p})
wc = pd.DataFrame(wc)
wc["p_holm"] = holm(wc.p)
wc.to_csv(out / "published_maps_within_class_plddt.tsv", sep="\t", index=False)
print(wc.round(4).to_string(index=False))

rw = []
for other in ("animal", "fungus"):
    share = S[S.kind == other].age_bin.value_counts(normalize=True).reindex(BINS).fillna(0)
    plant_mean = S[S.kind == "plant"].groupby("age_bin").plddt_global.mean().reindex(BINS)
    rw.append({"reference": other, "plant mean": S[S.kind == "plant"].plddt_global.mean(),
               "plant mean, reweighted to reference age composition": float((plant_mean * share).sum()),
               "reference mean": S[S.kind == other].plddt_global.mean()})
print(pd.DataFrame(rw).round(2).to_string(index=False))
pd.DataFrame(rw).to_csv(out / "published_maps_reweighting.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- C. balanced GenEra sample (our runs)
G = pd.read_csv(D / "Dataset_S14_published_gene_ages/genera_balanced_sample_ages.csv")
G = G[G.age_bin.isin(["all_life", "eukaryotes", "kingdom_only", "family_or_genus"])].copy()
G["group"] = G.tax_group7.replace({"Parasitic Worms": "Parasitic worms"})
rows = []
for g, d in G.groupby("group"):
    k = int((d.age_bin == "family_or_genus").sum())
    rows.append({"group": g, "n": len(d), "family or genus %": 100 * k / len(d),
                 "CI": tuple(round(100 * x, 1) for x in wilson(k, len(d))),
                 "all life %": 100 * (d.age_bin == "all_life").mean()})
    if g != "Viridiplantae":
        dd = G[G.group.isin(["Viridiplantae", g])]
        c, dof, p, v = chi2(pd.crosstab(dd.group, dd.age_bin))
        tests.append({"data": "balanced GenEra sample (nr, 2026)", "comparison": f"plants vs {g}", "chi2": c, "df": dof,
                      "p": p, "V": v})
print(pd.DataFrame(rows).round(1).to_string(index=False))
pd.DataFrame(rows).to_csv(out / "balanced_sample_groups.tsv", sep="\t", index=False)
T = pd.DataFrame(tests)
bal = T.data.str.startswith("balanced")
T.loc[bal, "p_holm"] = holm(T.loc[bal, "p"])
T.to_csv(out / "conservation_composition_tests.tsv", sep="\t", index=False)
print(T.round(4).to_string(index=False))
