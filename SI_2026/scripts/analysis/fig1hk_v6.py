"""Main-text Fig. 1h-k (Neff/len of plants vs human, yeast, E. coli; KDE + Jensen-Shannon distance) recomputed on the
revision alignments (Dataset S10), next to the original 34,921-row table with the same code; plus the assignment counts
of Dataset S10 (for the SI methods)."""
import numpy as np
import pandas as pd
from scipy.spatial.distance import jensenshannon
from scipy.stats import gaussian_kde, ks_2samp, spearmanr

B = "/path/to/project"
GRID = np.linspace(0, 120, 1000)


def jsd(a, b):
    pa, pb = gaussian_kde(a)(GRID), gaussian_kde(b)(GRID)
    return jensenshannon(pa / pa.sum(), pb / pb.sum(), base=2)


def report(label, groups):
    print(f"== {label}")
    for g, v in groups.items():
        print(f"  {g:6s} n={len(v):>7,}  median Neff/len {np.median(v):6.2f}  share >30 {100 * (v > 30).mean():5.1f}%")
    p = groups["PLANT"]
    for g in ("ECOLI", "HUMAN", "YEAST"):
        print(f"  plant vs {g:5s}: KS {ks_2samp(p, groups[g]).statistic:.3f}  JSD {jsd(p, groups[g]):.3f}")


# original table, original grouping (everything not human/E. coli/yeast = plant)
o = pd.read_csv(f"{B}/uniclust30_MSA/MSA_depth_with_plddt_34921.csv")
lab = o["organism_scientific"].map(lambda s: "HUMAN" if s.startswith("Homo") else "ECOLI" if s.startswith("Escherichia")
                                   else "YEAST" if s.startswith("Saccharomyces") else "PLANT")
report("original 34,921 rows (as in Fig. 1h-k)", {g: o.loc[lab == g, "Neff/len"].to_numpy() for g in ("PLANT", "HUMAN", "YEAST", "ECOLI")})

# revision: one row per protein of the same seven proteomes
m = pd.read_parquet(f"{B}/SI_2026/datasets/Dataset_S10_alignment_depth/merged_neff.parquet")
v6 = pd.read_parquet(f"{B}/SI_2026/datasets/Dataset_S6_AFDB_v6_47_proteomes/AF2_48_proteomes_v6.parquet", columns=["db_accession", "proteome_tar"])
v6["code"] = v6["proteome_tar"].str.extract(r"^UP\d+_\d+_([A-Z0-9]+)_v6")[0]
code = {"ARATH": "PLANT", "ORYSJ": "PLANT", "MAIZE": "PLANT", "SOYBN": "PLANT", "HUMAN": "HUMAN", "YEAST": "YEAST", "ECOLI": "ECOLI"}
v7 = v6[v6["code"].isin(code)].drop_duplicates("db_accession")
r = m.merge(v7[["db_accession", "code"]], on="db_accession")
r = r[np.isfinite(r["neff_per_len"])]
r["grp"] = r["code"].map(code)
report("revision, one row per protein, same seven proteomes", {g: r.loc[r.grp == g, "neff_per_len"].to_numpy() for g in ("PLANT", "HUMAN", "YEAST", "ECOLI")})
print("  rho(Neff/len, pLDDT) seven proteomes: %.3f" % spearmanr(r.neff_per_len, r.plddt_global)[0])

# Dataset S10 assignment counts (for the SI methods)
print("\n== Dataset S10 assignment")
print("proteins", len(m), "| with Neff", int(np.isfinite(m["neff"]).sum()))
print(m["how"].value_counts(dropna=False).to_dict())
pl = m[m["tax_group7"] == "Viridiplantae"]
print("plants", len(pl), pl["how"].value_counts(dropna=False).to_dict())
has = m[np.isfinite(m["neff"])]
print("distinct alignments used: all", has["rep_id"].nunique(), "| plants", has.loc[has.tax_group7 == "Viridiplantae", "rep_id"].nunique())
print("source:", m["source"].value_counts().to_dict())
mm = m[m["how"] == "mmseqs"]
print("mmseqs hits: median identity %.2f, median query coverage %.2f" % (mm["pident"].median(), mm["qcov"].median()))
