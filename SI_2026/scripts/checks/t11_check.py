import numpy as np, pandas as pd
from scipy.stats import spearmanr
m = pd.read_parquet("/path/to/project/SI_2026/datasets/Dataset_S10_alignment_depth/merged_neff.parquet",
                    columns=["tax_group7", "msa_raw_depth", "neff", "plddt_global"])
order = ["All", "Viridiplantae", "Fungi", "Animalia", "Protist", "Parasitic Worms", "Bacteria", "Archaea"]
for g in order:
    d = m if g == "All" else m[m.tax_group7 == g]
    a, b = d[d.msa_raw_depth.notna()], d[np.isfinite(d.neff)]
    print(f"{g:16s} n_depth {len(a):>7,}  n_neff {len(b):>7,}  median depth {a.msa_raw_depth.median():>8,.1f}"
          f"  rho depth {spearmanr(a.msa_raw_depth, a.plddt_global)[0]:.3f}  rho log10depth {spearmanr(np.log10(a.msa_raw_depth), a.plddt_global)[0]:.3f}"
          f"  rho neff {spearmanr(b.neff, b.plddt_global)[0]:.3f}  depth range {a.msa_raw_depth.min():,.0f}-{a.msa_raw_depth.max():,.0f}")
