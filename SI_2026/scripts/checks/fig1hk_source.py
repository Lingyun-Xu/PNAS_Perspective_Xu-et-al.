import numpy as np, pandas as pd
from scipy.spatial.distance import jensenshannon
from scipy.stats import gaussian_kde, ks_2samp
B = "/path/to/project"
GRID = np.linspace(0, 120, 1000)
def jsd(a, b):
    pa, pb = gaussian_kde(a)(GRID), gaussian_kde(b)(GRID)
    return jensenshannon(pa / pa.sum(), pb / pb.sum(), base=2)
o = pd.read_csv(f"{B}/uniclust30_MSA/MSA_depth_with_plddt_34921.csv")
o["kind"] = np.where(o["Source"].astype(str).str.match(r"^[0-9][a-z0-9]{3}_"), "pdb", "uc30")
o["grp"] = o["organism_scientific"].map(lambda s: "HUMAN" if s.startswith("Homo") else "ECOLI" if s.startswith("Escherichia") else "YEAST" if s.startswith("Saccharomyces") else "PLANT")
print(o.groupby(["grp", "kind"]).agg(rows=("Neff", "size"), med_neff=("Neff", "median"), med_nl=("Neff/len", "median")).round(2).to_string())
pp = o[(o.grp == "PLANT") & (o.kind == "pdb")]["Neff/len"].to_numpy()
for g in ("ECOLI", "HUMAN", "YEAST"):
    x = o[o.grp == g]["Neff/len"].to_numpy()
    print(f"plant PDB-chain rows only vs {g}: KS {ks_2samp(pp, x).statistic:.3f} JSD {jsd(pp, x):.3f}")
t = o[o.db_accession == "Q8WZ42"]
print("titin Q8WZ42: rows", len(t), "| alignments", t.Source.nunique(), "| Neff/len range %.1f-%.1f" % (t["Neff/len"].min(), t["Neff/len"].max()),
      "| pLDDT range %.1f-%.1f" % (t.plddt.min(), t.plddt.max()), "| chain lengths", sorted((t.Neff / t["Neff/len"]).round().astype(int).unique())[:8])
