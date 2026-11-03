"""What determines the AFDB v6 global pLDDT of plant proteins? Joint model of alignment depth, template distance,
disorder, gene age, length and species, for plant proteins without their own PDB structure (the proteins the
AI-guided selection targets).

Proteins: canonical AFDB v6 entries of Arabidopsis, rice and soybean (the plant species with a published GenEra map;
none exists for maize) without their own PDB structure, with an AFDB MSA depth, AIUPred disorder and a gene age.
Factors (groups of model terms):
  MSA depth         log10 AFDB MSA depth (sequences in the AFDB alignment, the model input)
  template distance model 1: ESM-2 distance to the closest template used by the model (models with a mapped template)
                    model 2: ESM-2 distance to the nearest PDB sequence at model building (all proteins; the quantity
                    the selection reduces)
  disorder          AIUPred disordered fraction
  gene age          published GenEra age class (all life, eukaryotes, kingdom, family or genus, species)
  length            log10 sequence length
  species           Arabidopsis, rice, soybean
Influence: OLS with natural cubic splines (4 df) for the continuous factors and indicators for the classes.
  R2 alone          R2 of the factor on its own
  unique            drop in R2 when the factor is removed from the full model
  share             Shapley share of the full-model R2 (average gain in R2 over all orders of adding the factors)
Effect: linear model (log10 depth, distance, disorder, log10 length linear); change in pLDDT from the 25th to the 75th
  percentile of each continuous factor, and of each age class versus all life, with the other factors held equal.
Uncertainty: 200 bootstrap samples of proteins (95% percentile intervals).
Usage: python plant_plddt_factors.py <out_dir>"""
import itertools
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from patsy import dmatrix

D = Path("/path/to/project/SI_2026/datasets")
out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(20261101)
B = 200
SPECIES = {"Arabidopsis thaliana": "Arabidopsis", "Oryza sativa subsp. japonica": "rice", "Glycine max": "soybean"}
AGES = ["all_life", "eukaryotes", "kingdom_only", "family", "species"]

dis = pd.read_parquet(D / "Dataset_S12_disorder/merged_per_protein.parquet",
                      columns=["db_accession", "organism_scientific", "plddt_global", "frac_disordered",
                               "msa_raw_depth", "length"])
dis = dis[~dis.db_accession.duplicated(keep=False) & dis.organism_scientific.isin(list(SPECIES))]
td = pd.read_parquet(D / "Dataset_S17_source_data_SI_figures/FigS06_S07/template_distances.parquet",
                     columns=["db_accession", "in_pdb", "d_template_min", "d_nearest_asof"])
age = pd.read_csv(D / "Dataset_S14_published_gene_ages/published_gene_ages_afdb_v6.tsv", sep="\t",
                  usecols=["db_accession", "age_bin"]).drop_duplicates("db_accession")
df = dis.merge(td, on="db_accession", how="inner").merge(age, on="db_accession", how="inner")
df["species"] = df.organism_scientific.map(SPECIES)
df = df[~df.in_pdb.astype(bool) & df.age_bin.isin(AGES) & df.msa_raw_depth.notna() & df.frac_disordered.notna()].copy()
df["log_depth"] = np.log10(df.msa_raw_depth.clip(lower=1))
df["log_len"] = np.log10(df.length)
print("plant proteins without own structure, with depth, disorder and age:", len(df), df.species.value_counts().to_dict(),
      flush=True)


def blocks(d, tcol, spline=True):
    """design-matrix columns of each factor group"""
    s = (lambda v: np.asarray(dmatrix("cr(x, df=4) - 1", {"x": v}, return_type="matrix"))) if spline else \
        (lambda v: np.asarray(v, float)[:, None])
    return {"MSA depth": s(d.log_depth.to_numpy()),
            "template distance": s(d[tcol].to_numpy()),
            "disorder": s(d.frac_disordered.to_numpy()),
            "gene age": pd.get_dummies(pd.Categorical(d.age_bin, categories=AGES)).iloc[:, 1:].to_numpy(float),
            "length": s(d.log_len.to_numpy()),
            "species": pd.get_dummies(pd.Categorical(d.species, categories=["Arabidopsis", "rice", "soybean"])).iloc[:, 1:].to_numpy(float)}


def r2_table(blk, y):
    """R2 of every subset of factor groups, from the Gram matrix of the centred design"""
    names = list(blk)
    X = np.column_stack([blk[n] for n in names])
    X = X - X.mean(axis=0)
    yc = y - y.mean()
    G, g, syy = X.T @ X, X.T @ yc, yc @ yc
    cols, k = {}, 0
    for n in names:
        cols[n] = list(range(k, k + blk[n].shape[1]))
        k += blk[n].shape[1]
    out = {(): 0.0}
    for r in range(1, len(names) + 1):
        for sub in itertools.combinations(names, r):
            idx = sum((cols[n] for n in sub), [])
            beta = np.linalg.lstsq(G[np.ix_(idx, idx)], g[idx], rcond=None)[0]
            out[tuple(sorted(sub))] = float(g[idx] @ beta / syy)
    return names, out


def summarize(names, r2):
    full = r2[tuple(sorted(names))]
    m = len(names)
    res = {}
    for j in names:
        others = [n for n in names if n != j]
        sh = 0.0
        for r in range(m):
            for sub in itertools.combinations(others, r):
                w = math.factorial(r) * math.factorial(m - r - 1) / math.factorial(m)
                sh += w * (r2[tuple(sorted(sub + (j,)))] - r2[tuple(sorted(sub))])
        res[j] = {"R2 alone": r2[(j,)], "unique": full - r2[tuple(sorted(others))], "Shapley": sh,
                  "share %": 100 * sh / full}
    return full, res


def linear_effects(d, tcol):
    X = pd.DataFrame({"log_depth": d.log_depth, "dist": d[tcol], "disorder": d.frac_disordered,
                      "log_len": d.log_len})
    for a in AGES[1:]:
        X[a] = (d.age_bin == a).astype(float)
    for s_ in ("rice", "soybean"):
        X[s_] = (d.species == s_).astype(float)
    X.insert(0, "const", 1.0)
    beta = np.linalg.lstsq(X.to_numpy(float), d.plddt_global.to_numpy(float), rcond=None)[0]
    b = dict(zip(X.columns, beta))
    q = {c: d[c].quantile(0.75) - d[c].quantile(0.25) for c in ("log_depth", tcol, "frac_disordered", "log_len")}
    yhat = X.to_numpy(float) @ beta
    r2 = 1 - ((d.plddt_global - yhat) ** 2).sum() / ((d.plddt_global - d.plddt_global.mean()) ** 2).sum()
    return {"MSA depth (25th to 75th percentile)": b["log_depth"] * q["log_depth"],
            "template distance (25th to 75th percentile)": b["dist"] * q[tcol],
            "disorder (25th to 75th percentile)": b["disorder"] * q["frac_disordered"],
            "length (25th to 75th percentile)": b["log_len"] * q["log_len"],
            "age: eukaryotes vs all life": b["eukaryotes"], "age: kingdom vs all life": b["kingdom_only"],
            "age: family or genus vs all life": b["family"], "age: species vs all life": b["species"],
            "species: rice vs Arabidopsis": b["rice"], "species: soybean vs Arabidopsis": b["soybean"],
            "linear R2": r2}


def run(d, tcol, tag):
    y = d.plddt_global.to_numpy(float)
    names, r2 = r2_table(blocks(d, tcol), y)
    full, res = summarize(names, r2)
    eff = linear_effects(d, tcol)
    boots_s, boots_e = [], []
    for _ in range(B):
        i = rng.integers(0, len(d), len(d))
        db = d.iloc[i]
        nb, r2b = r2_table(blocks(db, tcol), db.plddt_global.to_numpy(float))
        fb, rb = summarize(nb, r2b)
        boots_s.append({k: v["share %"] for k, v in rb.items()} | {"full R2": fb} |
                       {f"unique {k}": v["unique"] for k, v in rb.items()})
        boots_e.append(linear_effects(db, tcol))
    bs, be = pd.DataFrame(boots_s), pd.DataFrame(boots_e)
    rows = []
    for k, v in res.items():
        rows.append({"model": tag, "factor": k, "R2 alone": v["R2 alone"], "unique": v["unique"],
                     "unique lo": bs[f"unique {k}"].quantile(0.025), "unique hi": bs[f"unique {k}"].quantile(0.975),
                     "share %": v["share %"], "share lo": bs[k].quantile(0.025), "share hi": bs[k].quantile(0.975)})
    shares = pd.DataFrame(rows)
    effects = pd.DataFrame([{"model": tag, "effect": k, "value": v, "lo": be[k].quantile(0.025),
                             "hi": be[k].quantile(0.975)} for k, v in eff.items()])
    print(f"\n== {tag}: n = {len(d)}; full R2 = {full:.3f} (95% CI {bs['full R2'].quantile(0.025):.3f}–"
          f"{bs['full R2'].quantile(0.975):.3f}); linear R2 = {eff['linear R2']:.3f}", flush=True)
    print(shares.round(3).to_string(index=False), flush=True)
    print(effects.round(2).to_string(index=False), flush=True)
    return shares, effects, full


iqr = {}
m1 = df[np.isfinite(df.d_template_min)]
m2 = df[np.isfinite(df.d_nearest_asof)]
for c in ("log_depth", "d_template_min", "frac_disordered", "log_len"):
    iqr[c] = (m1[c].quantile(0.25), m1[c].quantile(0.75))
iqr["d_nearest_asof"] = (m2.d_nearest_asof.quantile(0.25), m2.d_nearest_asof.quantile(0.75))
print("IQRs (model 1; nearest-PDB from model 2):", {k: tuple(round(x, 3) for x in v) for k, v in iqr.items()},
      "| MSA depth IQR in sequences:", tuple(round(10 ** x) for x in iqr["log_depth"]), flush=True)
s1, e1, f1 = run(m1, "d_template_min", "model 1: closest template used")
s2, e2, f2 = run(m2, "d_nearest_asof", "model 2: nearest PDB sequence at model building")

# per species (model 1, without the species factor)
per = []
for sp, d in m1.groupby("species"):
    blk = blocks(d, "d_template_min")
    blk.pop("species")
    names, r2 = r2_table(blk, d.plddt_global.to_numpy(float))
    full, res = summarize(names, r2)
    per.append({"species": sp, "n": len(d), "full R2": full} | {f"share % {k}": v["share %"] for k, v in res.items()})
per = pd.DataFrame(per)
print("\n== per species (model 1)\n" + per.round(1).to_string(index=False), flush=True)
pd.concat([s1, s2]).to_csv(out / "plant_plddt_factor_shares.tsv", sep="\t", index=False)
pd.concat([e1, e2]).to_csv(out / "plant_plddt_factor_effects.tsv", sep="\t", index=False)
per.to_csv(out / "plant_plddt_factor_shares_by_species.tsv", sep="\t", index=False)
models = []
for tag, d, tcol, f in (("model 1", m1, "d_template_min", f1), ("model 2", m2, "d_nearest_asof", f2)):
    q = lambda c: (d[c].quantile(0.25), d[c].quantile(0.75))
    models.append({"model": tag, "n": len(d), "full R2": f,
                   "n by species": d.species.value_counts().to_dict(),
                   "MSA depth IQR (sequences)": tuple(round(10 ** x) for x in q("log_depth")),
                   "template distance IQR": tuple(round(x, 2) for x in q(tcol)),
                   "disorder IQR": tuple(round(x, 2) for x in q("frac_disordered")),
                   "length IQR (residues)": tuple(round(10 ** x) for x in q("log_len"))})
pd.DataFrame(models).to_csv(out / "plant_plddt_factor_models.tsv", sep="\t", index=False)
print(pd.DataFrame(models).to_string(index=False))
print("DONE")
