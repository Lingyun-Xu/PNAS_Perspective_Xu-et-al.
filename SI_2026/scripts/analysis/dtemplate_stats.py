#!/usr/bin/env python3
"""Plant versus other groups, and among the four plant proteomes, for the ESM-2 distance to templates
(concise PNAS SI Figs. S6-S7 and Tables S16-S18; full SI Section S{B7}.3). Same tests as stats_vm.py (imported):
  distributions: Kruskal-Wallis (epsilon^2); plants vs each group, Cliff's delta (Holm); among plants, Dunn-Holm
  correlation of global pLDDT with each distance per group and plant proteome (Fisher z vs plants, Holm; Cochran Q)
Measures (AFDB v6 models without their own PDB structure, template_distances.parquet = Dataset S17 FigS06_S07):
  d_template_min  distance to the closest template used by the model (906,513 models with a mapped template)
  d_nearest_asof  distance to the nearest PDB sequence available when the model was built (951,630 models)
  d_nearest_now   distance to the nearest PDB sequence in the September 2026 PDB (951,630 models; added 1 Nov 2026)
Usage: python dtemplate_stats.py <template_distances.parquet> <out_dir>
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import stats_vm as sv  # noqa: E402

MEASURES = (("d_template_min", "ESM-2 distance to closest template used (no own PDB structure)"),
            ("d_nearest_asof", "ESM-2 distance to nearest PDB sequence at model building (no own PDB structure)"),
            ("d_nearest_now", "ESM-2 distance to nearest PDB sequence, September 2026 (no own PDB structure)"))
cols = [c for c, _ in MEASURES]

src, out = Path(sys.argv[1]), Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)
df = pd.read_parquet(src, columns=["db_accession", "plddt_global", "tax_group7", "in_pdb", "organism_scientific"] + cols)
df["tax_group7"] = sv.norm_group(df.tax_group7.astype(str))
df["plant"] = df.organism_scientific.where(df.organism_scientific.isin(list(sv.PLANTS)))
nopdb = df[~df.in_pdb.astype(bool)].copy()
print("models", len(df), "without own PDB structure", len(nopdb),
      "with d_template_min", int(np.isfinite(nopdb.d_template_min).sum()), flush=True)
print(nopdb.groupby("tax_group7")[cols].median().round(3).to_string(), flush=True)

tests, among, corr = [], [], []
for col, label in MEASURES:
    g, a = sv.compare(nopdb, col, "tax_group7", label)
    tests.append(g)
    if a is not None:
        among.append(a)
    corr.append(sv.correlations(nopdb, col, "plddt_global", "tax_group7",
                                f"global pLDDT vs {col} (no own PDB structure)"))
corr = pd.concat(corr, ignore_index=True)
pd.concat(tests, ignore_index=True).to_csv(out / "dtemplate_group_tests.tsv", sep="\t", index=False)
pd.concat(among, ignore_index=True).to_csv(out / "dtemplate_among_plant_dunn.tsv", sep="\t", index=False)
corr.to_csv(out / "dtemplate_correlation_tests.tsv", sep="\t", index=False)
pd.set_option("display.width", 220)
print(pd.concat(tests, ignore_index=True).round(4).to_string(index=False), flush=True)
print(corr.round(4).to_string(index=False), flush=True)
print("DTEMPLATE_STATS_DONE", flush=True)
