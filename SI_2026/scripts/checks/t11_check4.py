import numpy as np, pandas as pd
from scipy.stats import spearmanr
m = pd.read_parquet("/path/to/project/SI_2026/datasets/Dataset_S10_alignment_depth/merged_neff.parquet", columns=["tax_group7", "neff", "neff_per_len", "plddt_global"])
p = m[(m.tax_group7 == "Viridiplantae") & np.isfinite(m.neff)]
print("rho Neff %.6f | rho Neff/len %.6f" % (spearmanr(p.neff, p.plddt_global)[0], spearmanr(p.neff_per_len, p.plddt_global)[0]))
