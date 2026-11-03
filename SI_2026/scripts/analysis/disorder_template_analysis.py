"""Disorder versus template distance as explanations of low pLDDT. Replaces the earlier regression, which used a template yes/no covariate
(93% of models have a template) and did not include template distance.

Proteins: canonical entries of the 47 AFDB v6 proteomes (Dataset S12, one model per accession) without their own PDB
structure, with a template mapped to an embedding and an AFDB MSA depth: 489,973 proteins. Template distance from
Dataset S17 (FigS06_S07/template_distances.parquet): ESM-2 distance to the closest template used (d_template_min).

1. Within each group and plant proteome: Spearman rho of global pLDDT with disordered fraction, template distance and
   MSA depth, and the partial Spearman rho of each with the other two held fixed (Pearson correlation of the rank
   residuals; t test with n - 4 degrees of freedom).
2. Mean pLDDT of each group minus that of plants: observed; and among proteins of the same class of the named measures
   (classes: deciles of each measure over all 489,973 proteins). Within each class the group mean minus the plant mean
   is taken and the class differences are averaged, weighted by the number of plant proteins per class; classes
   without proteins of both groups are left out (coverage reported). Ordered-residue pLDDT the same way.
3. Sensitivity: 20 classes per measure; medians (group proteins weighted to the plant mix of classes, weighted median
   minus plant median); protein length as a fourth matched measure (quintiles of all four measures).
Usage: python disorder_template_analysis.py <out_dir>"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

D = Path("/path/to/project/SI_2026/datasets")
out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)
PLANTS = {"Arabidopsis thaliana": "Arabidopsis", "Glycine max": "soybean", "Oryza sativa subsp. japonica": "rice",
          "Zea mays": "maize"}
GROUPS = ["Viridiplantae", "Archaea", "Bacteria", "Protist", "Fungi", "Animalia", "Parasitic worms"]
MEAS = [("frac_disordered", "disordered fraction"), ("d_template_min", "template distance"),
        ("log_depth", "MSA depth")]

dis = pd.read_parquet(D / "Dataset_S12_disorder/merged_per_protein.parquet",
                      columns=["db_accession", "organism_scientific", "plddt_global", "plddt_ordered",
                               "frac_disordered", "msa_raw_depth", "length"])
dis = dis[~dis.db_accession.duplicated(keep=False)]
td = pd.read_parquet(D / "Dataset_S17_source_data_SI_figures/FigS06_S07/template_distances.parquet",
                     columns=["db_accession", "tax_group7", "in_pdb", "d_template_min"])
df = dis.merge(td, on="db_accession", how="inner")
df["group"] = df.tax_group7.replace({"Parasitic Worms": "Parasitic worms"})
df = df[~df.in_pdb.astype(bool) & np.isfinite(df.d_template_min) & df.msa_raw_depth.notna()].copy()
df["log_depth"] = np.log10(df.msa_raw_depth.clip(lower=1))
df["plant"] = df.organism_scientific.map(PLANTS)
print("proteins:", len(df), df.group.value_counts().reindex(GROUPS).to_dict(), flush=True)


# ---------------------------------------------------------------- 1. correlations
def partial_rank_corr(d, x, controls):
    r = d[["plddt_global", x] + controls].rank().to_numpy(float)
    z = np.column_stack([np.ones(len(r)), r[:, 2:]])
    res = [r[:, k] - z @ np.linalg.lstsq(z, r[:, k], rcond=None)[0] for k in (0, 1)]
    pr = float(np.corrcoef(res[0], res[1])[0, 1])
    t = pr * np.sqrt((len(d) - 4) / (1 - pr ** 2))
    return pr, float(2 * stats.t.sf(abs(t), len(d) - 4))


rows = []
for name, d in [(g, df[df.group == g]) for g in GROUPS] + [(f"plant: {p}", df[df.plant == p]) for p in PLANTS.values()]:
    row = {"subset": name, "n": len(d), "n ordered-residue pLDDT": int(d.plddt_ordered.notna().sum())}
    for col, label in MEAS:
        row[f"rho {label}"] = stats.spearmanr(d.plddt_global, d[col]).statistic
        row[f"partial rho {label}"], row[f"partial p {label}"] = partial_rank_corr(d, col, [c for c, _ in MEAS if c != col])
    rows.append(row)
corr = pd.DataFrame(rows)
corr.to_csv(out / "disorder_vs_template_correlations.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
print(corr.round(3).to_string(index=False), flush=True)


# ---------------------------------------------------------------- 2. matched differences from plants
# classes over all 489,973 proteins (tied values share a class): deciles; 20 classes and quintiles for the sensitivity
for q, cols in ((10, [c for c, _ in MEAS]), (20, [c for c, _ in MEAS]), (5, [c for c, _ in MEAS] + ["length"])):
    for c in cols:
        df[f"q{q}_{c}"] = pd.qcut(df[c], q, labels=False, duplicates="drop")
        print(f"q{q}_{c}: {df[f'q{q}_{c}'].nunique()} classes", flush=True)


def classes(cols, q):
    return [f"q{q}_{c}" for c in cols]


def matched(ref, oth, y, cols, q=10, median=False):
    """group minus plants; returns (difference, share of plant proteins covered, share of group proteins covered)"""
    ref, oth = ref[ref[y].notna()], oth[oth[y].notna()]
    if not cols:
        f = np.median if median else np.mean
        return float(f(oth[y]) - f(ref[y])), 1.0, 1.0
    keys = classes(cols, q)
    nr = ref.groupby(keys).size().rename("nr")
    no = oth.groupby(keys).size().rename("no")
    both = pd.concat([nr, no], axis=1).dropna()
    rr = ref.join(both, on=keys, how="inner")
    oo = oth.join(both, on=keys, how="inner")
    w = (oo.nr / oo.no).to_numpy(float)                    # group proteins weighted to the plant mix of classes
    if median:
        v = oo[y].to_numpy(float)
        o = np.argsort(v)
        c = np.cumsum(w[o])
        diff = float(v[o][np.searchsorted(c, c[-1] / 2)] - rr[y].median())
    else:
        diff = float(np.average(oo[y], weights=w) - rr[y].mean())     # = plant-weighted mean of class differences
    return diff, len(rr) / len(ref), len(oo) / len(oth)


COMPARISONS = [("observed", [], "plddt_global"), ("same disorder", ["frac_disordered"], "plddt_global"),
               ("same template distance and MSA depth", ["d_template_min", "log_depth"], "plddt_global"),
               ("same template distance, MSA depth and disorder", ["d_template_min", "log_depth", "frac_disordered"],
                "plddt_global"),
               ("ordered-residue, observed", [], "plddt_ordered"),
               ("ordered-residue, same template distance and MSA depth", ["d_template_min", "log_depth"],
                "plddt_ordered")]
plants = df[df.group == "Viridiplantae"]
res, sens = [], []
for g in GROUPS[1:]:
    o = df[df.group == g]
    r = {"group": g, "n": len(o), "n ordered-residue pLDDT": int(o.plddt_ordered.notna().sum())}
    s = {"group": g}
    for tag, cols, y in COMPARISONS:
        diff, cp, co = matched(plants, o, y, cols)
        r[tag] = diff
        r[f"{tag}: coverage plants"], r[f"{tag}: coverage group"] = cp, co
        if cols:
            s[f"{tag}: 20 classes"] = matched(plants, o, y, cols, q=20)[0]
        s[f"{tag}: medians"] = matched(plants, o, y, cols, median=True)[0]
        if cols:
            s[f"{tag} + length (quintiles)"] = matched(plants, o, y, cols + ["length"], q=5)[0]
    for p in PLANTS.values():                              # each plant proteome against this group
        pp = plants[plants.plant == p]
        for tag, cols, y in COMPARISONS[:4]:
            r[f"{p}: {tag}"] = matched(pp, o, y, cols)[0]
    res.append(r)
    sens.append(s)
diffs = pd.DataFrame(res)
diffs.to_csv(out / "difference_from_plants_matched.tsv", sep="\t", index=False)
pd.DataFrame(sens).to_csv(out / "difference_from_plants_sensitivity.tsv", sep="\t", index=False)
print(diffs[["group", "n"] + [t for t, _, _ in COMPARISONS]].round(2).to_string(index=False), flush=True)
print(pd.DataFrame(sens).set_index("group").T.round(2).to_string(), flush=True)
print("DONE")
