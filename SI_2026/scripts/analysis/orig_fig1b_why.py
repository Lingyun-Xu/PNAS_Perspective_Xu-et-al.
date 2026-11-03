"""Why the original Fig. 1b correlation was weak (rho = 0.20): composition of the 34,921-row table and the same
measure on the revision alignments (Dataset S10), restricted step by step."""
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

B = "/path/to/project"
o = pd.read_csv(f"{B}/uniclust30_MSA/MSA_depth_with_plddt_34921.csv")
o["kind"] = np.where(o["Source"].astype(str).str.match(r"^[0-9][a-z0-9]{3}_"), "pdb_chain", "uniclust30")


def rho(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    k = np.isfinite(x) & np.isfinite(y)
    return f"rho={spearmanr(x[k], y[k])[0]:.3f} n={int(k.sum()):,}"


print("== original table")
print("rows", len(o), o["kind"].value_counts().to_dict())
print("alignments", o["Source"].nunique(), o.groupby("kind")["Source"].nunique().to_dict())
print("accessions", o["db_accession"].nunique(), o.groupby("kind")["db_accession"].nunique().to_dict())
print("alignment-accession pairs", len(o[["Source", "db_accession"]].drop_duplicates()))
print(o.groupby("organism_scientific").agg(rows=("plddt", "size"), proteins=("db_accession", "nunique"),
                                           alignments=("Source", "nunique")).to_string())
for f, col in (("pdb_MSA.csv", None), ("uniclust30_MSA.csv", None)):
    lst = pd.read_csv(f"{B}/uniclust30_MSA/{f}")
    print(f, len(lst), list(lst.columns)[:3], lst.iloc[:2, 0].tolist())
print("Neff/len vs pLDDT (v4), all rows:", rho(o["Neff/len"], o["plddt"]))
print("Neff vs pLDDT (v4), all rows:", rho(o["Neff"], o["plddt"]))
print("median Neff", round(o["Neff"].median()), "median pLDDT", o["plddt"].median(), "mean pLDDT", round(o["plddt"].mean(), 2))

a = o.groupby("db_accession").agg(nl=("Neff/len", "median"), neff=("Neff", "median"), pl4=("plddt", "median"),
                                  org=("organism_scientific", "first"))
print("one row per protein: Neff/len", rho(a.nl, a.pl4), "| Neff", rho(a.neff, a.pl4))

# ---------------------------------------------------------------- revision alignments, same seven proteomes
m = pd.read_parquet(f"{B}/SI_2026/datasets/Dataset_S10_alignment_depth/merged_neff.parquet",
                    columns=["db_accession", "length", "source", "how", "plddt_global", "neff", "neff_per_len",
                             "query_length", "msa_raw_depth"])
v6 = pd.read_parquet(f"{B}/SI_2026/datasets/Dataset_S6_AFDB_v6_47_proteomes/AF2_48_proteomes_v6.parquet",
                     columns=["db_accession", "proteome_tar"])
v6["code"] = v6["proteome_tar"].str.extract(r"^UP\d+_\d+_([A-Z0-9]+)_v6")[0]
SEVEN = {"HUMAN", "YEAST", "ECOLI", "ARATH", "ORYSJ", "MAIZE", "SOYBN"}
v7 = v6[v6["code"].isin(SEVEN)].drop_duplicates("db_accession")[["db_accession", "code"]]
m7 = m.merge(v7, on="db_accession")
m7 = m7[np.isfinite(m7["neff"])]
m7["in_orig"] = m7["db_accession"].isin(a.index)
print("\n== revision alignments (Dataset S10), seven proteomes")
print("proteins with Neff", len(m7), "of which in original set", int(m7["in_orig"].sum()))
for lab, d in (("all seven proteomes", m7), ("original proteins", m7[m7.in_orig]), ("other proteins", m7[~m7.in_orig])):
    print(f"{lab:22s} Neff/len {rho(d.neff_per_len, d.plddt_global)} | log Neff {rho(d.neff, d.plddt_global)}"
          f" | median Neff {d.neff.median():.0f} | Neff>=1000 {100 * (d.neff >= 1000).mean():.1f}%"
          f" | median pLDDT {d.plddt_global.median():.1f} | pLDDT>=70 {100 * (d.plddt_global >= 70).mean():.1f}%"
          f" | pLDDT IQR {d.plddt_global.quantile(.25):.1f}-{d.plddt_global.quantile(.75):.1f}")
for code in ["HUMAN", "ECOLI", "YEAST", "ARATH", "ORYSJ", "MAIZE", "SOYBN"]:
    d = m7[m7.code == code]
    print(f"  {code}: all {rho(d.neff_per_len, d.plddt_global)} | original {rho(d[d.in_orig].neff_per_len, d[d.in_orig].plddt_global)}")

# original Neff/len against v6 pLDDT (pLDDT version check), per protein
j = a.join(m7.set_index("db_accession")[["plddt_global", "length", "neff", "neff_per_len"]], how="inner", rsuffix="_rev")
print("\noriginal Neff/len vs v6 pLDDT (per protein):", rho(j.nl, j.plddt_global))
print("revision Neff/len vs v6 pLDDT, same proteins:", rho(j.neff_per_len, j.plddt_global))

# ---------------------------------------------------------------- PDB chain coverage of the full-length protein
p = o[o.kind == "pdb_chain"].drop_duplicates(["Source", "db_accession"]).copy()
p["qlen"] = p["Neff"] / p["Neff/len"]
p = p.merge(m7[["db_accession", "length", "plddt_global"]], on="db_accession")
p["cov"] = (p["qlen"] / p["length"]).clip(upper=1)
print("\n== PDB-chain alignments: chain length / protein length")
print("pairs", len(p), "median coverage", round(p["cov"].median(), 3), "share < 0.5:", round(100 * (p["cov"] < 0.5).mean(), 1), "%")
for lo, hi in ((0, 0.5), (0.5, 0.9), (0.9, 1.01)):
    d = p[(p["cov"] >= lo) & (p["cov"] < hi)]
    print(f"  coverage {lo}-{hi}: Neff/len vs v6 pLDDT {rho(d['Neff/len'], d.plddt_global)} | Neff {rho(d['Neff'], d.plddt_global)}")
print("chain length median", round(p["qlen"].median()), "| chains < 50 residues", round(100 * (p["qlen"] < 50).mean(), 1), "%",
      "| protein length median", int(p["length"].median()))
print("chains per protein: median", p.groupby("db_accession").size().median(), "max", p.groupby("db_accession").size().max())
best = p.sort_values("cov", ascending=False).drop_duplicates("db_accession")
print("one chain per protein (best coverage): Neff/len", rho(best["Neff/len"], best.plddt_global),
      "| Neff", rho(best["Neff"], best.plddt_global))
j["neff_per_protein_len"] = j["neff"] / j["length"]
print("per protein, original Neff / full protein length vs v6 pLDDT:", rho(j["neff_per_protein_len"], j.plddt_global))
print("per protein, original Neff vs v6 pLDDT:", rho(j["neff"], j.plddt_global))
print("human rows per protein: median", o[o.organism_scientific == "Homo sapiens"].groupby("db_accession").size().median(),
      "| rows that repeat an alignment", int(len(o) - len(o[["Source", "db_accession"]].drop_duplicates())))
dup = o.groupby(["Source", "db_accession"]).agg(n=("plddt", "size"), spread=("plddt", lambda s: s.max() - s.min()))
print("alignments repeated:", int((dup.n > 1).sum()), "| median pLDDT spread among repeats", round(dup[dup.n > 1].spread.median(), 1))
