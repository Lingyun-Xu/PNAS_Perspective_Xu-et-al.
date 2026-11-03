import numpy as np, pandas as pd
from scipy.stats import spearmanr
B = "/path/to/project"
o = pd.read_csv(f"{B}/uniclust30_MSA/MSA_depth_with_plddt_34921.csv")
o["kind"] = np.where(o["Source"].astype(str).str.match(r"^[0-9][a-z0-9]{3}_"), "pdb_chain", "uniclust30")
a = o.groupby("db_accession").agg(nl=("Neff/len", "median"), neff=("Neff", "median"), pl=("plddt", "median"))
print("per protein v4: Neff/len %.4f  Neff %.4f" % (spearmanr(a.nl, a.pl)[0], spearmanr(a.neff, a.pl)[0]))
pc = o[o.kind == "pdb_chain"]
n = pc.groupby("db_accession")["Source"].nunique()
print("PDB chains per protein: median", n.median(), "max", n.max(), "proteins", len(n), "| proteins with >1 chain", int((n > 1).sum()))
dup = o.groupby(["Source", "db_accession", "organism_scientific"]).size().rename("n").reset_index()
r = dup[dup.n > 1]
print("repeated alignments", len(r), "proteins", r.db_accession.nunique(), r.organism_scientific.value_counts().to_dict(),
      "extra rows", int((r.n - 1).sum()), "max repeats", int(r.n.max()))
m = pd.read_parquet(f"{B}/SI_2026/datasets/Dataset_S10_alignment_depth/merged_neff.parquet", columns=["db_accession", "length"])
p = pc.groupby(["Source", "db_accession"]).agg(neff=("Neff", "first"), nl=("Neff/len", "first"), pl=("plddt", "median")).reset_index()
p["qlen"] = p.neff / p.nl
p = p.merge(m.drop_duplicates("db_accession"), on="db_accession")
p["cov"] = (p.qlen / p.length).clip(upper=1)
print("pairs", len(p), "median cov %.3f" % p["cov"].median(), "cov<0.5 %.1f%%" % (100 * (p["cov"] < .5).mean()), "cov>=0.9 %.1f%%" % (100 * (p["cov"] >= .9).mean()))
for lo, hi in ((0, .5), (.5, .9), (.9, 1.01)):
    d = p[(p["cov"] >= lo) & (p["cov"] < hi)]
    print("cov %.1f-%.1f: n=%d Neff/len vs v4 pLDDT %.3f" % (lo, hi, len(d), spearmanr(d.nl, d.pl)[0]))
