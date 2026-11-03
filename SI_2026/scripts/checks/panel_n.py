import pandas as pd
p = pd.read_parquet("/path/to/project/SI_2026/datasets/Dataset_S13_homolog_class/panel_strata.parquet")
print(len(p), list(p.columns))
print(p.tax_group7.value_counts().to_dict())
for g in ("Viridiplantae", "Animalia", "Fungi"):
    s = p[p.tax_group7 == g]
    print(g, len(s), (s.panel_bin.value_counts(normalize=True) * 100).round(1).to_dict())
