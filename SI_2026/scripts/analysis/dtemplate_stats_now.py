"""Tables S16-S18 with the nearest-PDB distance in the September 2026 PDB (d_nearest_now, Fig. S7 D-F) added, with the
tests of dtemplate_stats.py (stats_vm.py functions). Also recomputes the ESMFold rows of the main statistics files on
all 92,385 proteins (Dataset S11). Writes to <out_dir>; the caller compares with SI_2026/statistics and installs."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "/path/to/project/template_esm/campaign/stats")
import stats_vm as sv  # noqa: E402

D = Path("/path/to/project/SI_2026/datasets")
ST = Path("/path/to/project/SI_2026/statistics")
out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------- template distances (Tables S16-S18)
df = pd.read_parquet(D / "Dataset_S17_source_data_SI_figures/FigS06_S07/template_distances.parquet",
                     columns=["db_accession", "plddt_global", "tax_group7", "in_pdb", "organism_scientific",
                              "d_template_min", "d_nearest_asof", "d_nearest_now"])
df["tax_group7"] = sv.norm_group(df.tax_group7.astype(str))
df["plant"] = df.organism_scientific.where(df.organism_scientific.isin(list(sv.PLANTS)))
nopdb = df[~df.in_pdb.astype(bool)].copy()
print("models", len(df), "| without own PDB structure", len(nopdb),
      "| with d_nearest_now", int(np.isfinite(nopdb.d_nearest_now).sum()), flush=True)
MEAS = [("d_template_min", "ESM-2 distance to closest template used (no own PDB structure)",
         "global pLDDT vs d_template_min (no own PDB structure)"),
        ("d_nearest_asof", "ESM-2 distance to nearest PDB sequence at model building (no own PDB structure)",
         "global pLDDT vs d_nearest_asof (no own PDB structure)"),
        ("d_nearest_now", "ESM-2 distance to nearest PDB sequence, September 2026 (no own PDB structure)",
         "global pLDDT vs d_nearest_now (no own PDB structure)")]
tests, among, corr = [], [], []
for col, label, clabel in MEAS:
    g, a = sv.compare(nopdb, col, "tax_group7", label)
    tests.append(g)
    among.append(a)
    corr.append(sv.correlations(nopdb, col, "plddt_global", "tax_group7", clabel))
    print("done", col, flush=True)
pd.concat(tests, ignore_index=True).to_csv(out / "dtemplate_group_tests.tsv", sep="\t", index=False)
pd.concat(among, ignore_index=True).to_csv(out / "dtemplate_among_plant_dunn.tsv", sep="\t", index=False)
pd.concat(corr, ignore_index=True).to_csv(out / "dtemplate_correlation_tests.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- ESMFold rows on all 92,385 proteins
es = pd.read_csv(D / "Dataset_S11_ESMFold_confidence/esmfold_vs_afdb_v6_per_protein.tsv", sep="\t")
es["tax_group7"] = sv.norm_group(es.tax_group7)
code = es.proteome_tar.str.extract(r"^UP\d+_\d+_([A-Z0-9]+)_v6")[0]
es["plant"] = code.map({"ARATH": "Arabidopsis thaliana", "SOYBN": "Glycine max",
                        "ORYSJ": "Oryza sativa subsp. japonica", "MAIZE": "Zea mays"})
es["plddt_global"] = es.afdb_v6_plddt_global
es["af2_minus_esmfold"] = es.plddt_global - es.esmfold_plddt
eg, ea = [], []
for col, label in (("esmfold_plddt", "ESMFold pLDDT (ESM Atlas)"), ("af2_minus_esmfold", "AF2 - ESMFold pLDDT")):
    g, a = sv.compare(es, col, "tax_group7", label)
    eg.append(g)
    ea.append(a)
ec = sv.correlations(es, "plddt_global", "esmfold_plddt", "tax_group7", "AF2 vs ESMFold pLDDT")
for name, new_rows in (("group_tests.tsv", pd.concat(eg)), ("among_plant_dunn.tsv", pd.concat(ea)),
                       ("correlation_tests.tsv", ec)):
    old = pd.read_csv(ST / name, sep="\t")
    keep = old[~old.measure.isin(set(new_rows.measure))]
    merged = pd.concat([keep, new_rows[old.columns.intersection(new_rows.columns)]], ignore_index=True)
    merged.to_csv(out / name, sep="\t", index=False)
    print(f"{name}: replaced {len(old) - len(keep)} ESMFold rows with {len(new_rows)}", flush=True)
print("DONE")
