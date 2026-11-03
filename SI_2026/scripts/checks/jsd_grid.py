import numpy as np, pandas as pd
from scipy.spatial.distance import jensenshannon
from scipy.stats import gaussian_kde
B = "/path/to/project"
m = pd.read_parquet(f"{B}/SI_2026/datasets/Dataset_S10_alignment_depth/merged_neff.parquet", columns=["db_accession", "neff_per_len"])
v6 = pd.read_parquet(f"{B}/SI_2026/datasets/Dataset_S6_AFDB_v6_47_proteomes/AF2_48_proteomes_v6.parquet", columns=["db_accession", "proteome_tar"])
v6["code"] = v6["proteome_tar"].str.extract(r"^UP\d+_\d+_([A-Z0-9]+)_v6")[0]
code = {"ARATH": "PLANT", "ORYSJ": "PLANT", "MAIZE": "PLANT", "SOYBN": "PLANT", "HUMAN": "HUMAN", "YEAST": "YEAST", "ECOLI": "ECOLI"}
r = m.merge(v6[v6.code.isin(code)].drop_duplicates("db_accession"), on="db_accession")
r = r[np.isfinite(r.neff_per_len)]
r["grp"] = r.code.map(code)
k = {g: gaussian_kde(r.loc[r.grp == g, "neff_per_len"].to_numpy()) for g in ("PLANT", "HUMAN", "YEAST", "ECOLI")}
for hi in (30, 60, 120):
    grid = np.linspace(0, hi, 1000)
    d = {g: k[g](grid) for g in k}
    print(f"grid 0-{hi}: " + "  ".join(f"plant-{g} {jensenshannon(d['PLANT'] / d['PLANT'].sum(), d[g] / d[g].sum(), base=2):.3f}" for g in ("ECOLI", "HUMAN", "YEAST")))
