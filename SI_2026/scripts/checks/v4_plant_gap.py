import numpy as np, pandas as pd
d = pd.read_parquet("/path/to/project/SI_2026/datasets/AFDB_v4_extracted/AF2_48_model_oganisms.parquet", columns=["AAseq_one_letter_code", "organism_scientific"])
d = d[d.organism_scientific.str.contains("Arabidopsis|Glycine max|Oryza sativa|Zea mays")].copy()
d["seq"] = d.AAseq_one_letter_code.str.lstrip(";").str.replace(r"\s+", "", regex=True)
u = d.drop_duplicates("seq")
L = u.seq.str.len()
print("unique", len(u), "| non-standard letters:", int(u.seq.str.contains("[^ACDEFGHIKLMNPQRSTVWY]").sum()))
for cut in (1000, 1022, 1024, 1500, 2000, 2500, 2700, 3000):
    print(f"  length <= {cut}: {(L <= cut).sum():,}")
target = 161460
s = np.sort(L.to_numpy())
print("length cutoff giving 161,460:", s[target - 1], "->", s[target] if target < len(s) else None)
