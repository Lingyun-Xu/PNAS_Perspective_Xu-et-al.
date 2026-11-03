"""Why are plant proteins less confident than those of other groups? Decomposition of the difference in mean pLDDT
between another group and plants into the parts due to differences in template distance, alignment depth, disorder,
length (and gene age, where available), and a part not explained by these factors (Oaxaca-Blinder decomposition with
pooled coefficients; Oaxaca 1973, Blinder 1973).

For plants and one other group G: pLDDT ~ f(template distance) + f(log10 MSA depth) + f(disordered fraction)
+ f(log10 length) [+ gene-age classes] + plant indicator, fitted by ordinary least squares on the proteins of both groups
(f = natural cubic spline, 4 df, knots from the pooled proteins). Part due to factor k = sum over its terms of
coefficient x (mean term value in G - mean in plants); unexplained part = -(plant coefficient). The parts add up
exactly to the observed gap (mean pLDDT of G minus that of plants).
Outcomes: global pLDDT; ordered-residue pLDDT (mean over residues with AIUPred score <= 0.5).
Proteins: AFDB v6 entries of the 47 proteomes without their own PDB structure, with a template mapped to an embedding
(d_template_min), an AFDB MSA depth and AIUPred disorder. Gene-age runs: the eight species with published GenEra maps
(plants: Arabidopsis, rice, soybean; animals: mouse, fly, C. elegans; fungi: budding and fission yeast).
Uncertainty: 200 bootstrap samples of proteins (main runs). Checks (point estimates): coefficients from plants only or
from G only; no length in the model.
Usage: python plddt_gap_decomposition.py <out_dir>"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from patsy import dmatrix

D = Path("/path/to/project/SI_2026/datasets")
out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(20261102)
B = 200
AGES = ["all_life", "eukaryotes", "kingdom_only", "family", "species"]
FACTORS = [("template distance", "d_template_min"), ("alignment depth", "log_depth"), ("disorder", "frac_disordered"),
           ("length", "log_len")]

dis = pd.read_parquet(D / "Dataset_S12_disorder/merged_per_protein.parquet",
                      columns=["db_accession", "plddt_global", "plddt_ordered", "frac_disordered", "msa_raw_depth",
                               "length"])
dis = dis[~dis.db_accession.duplicated(keep=False)]
td = pd.read_parquet(D / "Dataset_S17_source_data_SI_figures/FigS06_S07/template_distances.parquet",
                     columns=["db_accession", "in_pdb", "d_template_min", "tax_group7"])
df = dis.merge(td, on="db_accession", how="inner")
df["group"] = df.tax_group7.replace({"Parasitic Worms": "Parasitic worms"})
df = df[~df.in_pdb.astype(bool) & np.isfinite(df.d_template_min) & df.msa_raw_depth.notna()
        & df.frac_disordered.notna()].copy()
df["log_depth"] = np.log10(df.msa_raw_depth.clip(lower=1))
df["log_len"] = np.log10(df.length)
age = pd.read_csv(D / "Dataset_S14_published_gene_ages/published_gene_ages_afdb_v6.tsv", sep="\t",
                  usecols=["db_accession", "proteome", "age_bin"]).drop_duplicates("db_accession")
sp = df.merge(age, on="db_accession", how="inner")
sp = sp[sp.age_bin.isin(AGES)].copy()
KIND = {"ARATH": "Viridiplantae", "ORYSJ": "Viridiplantae", "SOYBN": "Viridiplantae", "MOUSE": "Animalia",
        "DROME": "Animalia", "CAEEL": "Animalia", "YEAST": "Fungi", "SCHPO": "Fungi"}
sp["group"] = sp.proteome.map(KIND)


def design(d, factors, with_age):
    blk = {n: np.asarray(dmatrix("cr(x, df=4) - 1", {"x": d[c].to_numpy(float)}, return_type="matrix"))
           for n, c in factors}
    if with_age:
        blk["gene age"] = pd.get_dummies(pd.Categorical(d.age_bin, categories=AGES)).iloc[:, 1:].to_numpy(float)
    return blk


def decompose(d, outcome, factors, with_age, ref="pooled"):
    plant = (d.group == "Viridiplantae").to_numpy()
    blk = design(d, factors, with_age)
    names = list(blk)
    X = np.column_stack([blk[n] for n in names])
    y = d[outcome].to_numpy(float)
    if ref == "pooled":
        beta = np.linalg.lstsq(np.column_stack([np.ones(len(y)), X, plant.astype(float)]), y, rcond=None)[0]
        b, unexplained = beta[1:-1], -beta[-1]
    else:
        s = plant if ref == "plants" else ~plant
        b = np.linalg.lstsq(np.column_stack([np.ones(s.sum()), X[s]]), y[s], rcond=None)[0][1:]
        unexplained = None
    diff = X[~plant].mean(axis=0) - X[plant].mean(axis=0)
    res, k = {"gap": float(y[~plant].mean() - y[plant].mean())}, 0
    for n in names:
        w = blk[n].shape[1]
        res[n] = float(b[k:k + w] @ diff[k:k + w])
        k += w
    res["unexplained"] = float(unexplained) if unexplained is not None else res["gap"] - sum(
        v for n, v in res.items() if n != "gap")
    return res


rows, medians = [], []


def run(run_id, label, protein_set, d, outcome, with_age, boot=True):
    d = d[d[outcome].notna()]
    est = decompose(d, outcome, FACTORS, with_age)
    r = {"run": run_id, "comparison": label, "protein set": protein_set, "outcome": outcome, "with gene age": with_age,
         "n plants": int((d.group == "Viridiplantae").sum()), "n other": int((d.group != "Viridiplantae").sum()),
         "mean plants": float(d.loc[d.group == "Viridiplantae", outcome].mean()),
         "mean other": float(d.loc[d.group != "Viridiplantae", outcome].mean())}
    r.update(est)
    if boot:
        bt = []
        for _ in range(B):
            i = rng.integers(0, len(d), len(d))
            bt.append(decompose(d.iloc[i], outcome, FACTORS, with_age))
        bt = pd.DataFrame(bt)
        for k in est:
            r[f"{k} lo"], r[f"{k} hi"] = float(bt[k].quantile(0.025)), float(bt[k].quantile(0.975))
    for ref in ("plants", "other"):
        for k, v in decompose(d, outcome, FACTORS, with_age, ref).items():
            if k != "gap":
                r[f"{k} [{ref} coefficients]"] = v
    for k, v in decompose(d, outcome, FACTORS[:3], with_age).items():
        if k != "gap":
            r[f"{k} [no length]"] = v
    rows.append(r)
    for g, part in d.groupby("group"):
        medians.append({"run": run_id, "group": g, "n": len(part), "median pLDDT": part[outcome].median(),
                        "median template distance": part.d_template_min.median(),
                        "median MSA depth": part.msa_raw_depth.median(),
                        "median disordered fraction": part.frac_disordered.median(),
                        "median length": part.length.median()})
    print(run_id, label, {k: round(v, 2) for k, v in est.items()}, flush=True)


for g in ("Animalia", "Fungi"):
    pair = df[df.group.isin(["Viridiplantae", g])]
    run(f"global_{g}", f"{g} - plants", "47 proteomes", pair, "plddt_global", False)
    run(f"ordered_{g}", f"{g} - plants", "47 proteomes", pair, "plddt_ordered", False)
    run(f"age8_{g}", f"{g} - plants", "eight species", sp[sp.group.isin(["Viridiplantae", g])], "plddt_global", True)
for g in ("Parasitic worms", "Protist", "Bacteria", "Archaea"):
    run(f"global_{g}", f"{g} - plants", "47 proteomes", df[df.group.isin(["Viridiplantae", g])], "plddt_global",
        False)

R = pd.DataFrame(rows)
R.to_csv(out / "plddt_gap_decomposition.tsv", sep="\t", index=False)
M = pd.DataFrame(medians)
M.to_csv(out / "plddt_gap_factor_medians.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
show = ["run", "n plants", "n other", "gap", "template distance", "alignment depth", "disorder", "length", "gene age",
        "unexplained"]
print(R[[c for c in show if c in R.columns]].round(2).to_string(index=False))
print(R[["run"] + [c for c in R.columns if "[" in c]].round(2).T.to_string())
print(M.round(3).to_string(index=False))
