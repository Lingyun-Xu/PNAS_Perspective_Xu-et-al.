import sys
import pandas as pd
sys.path.insert(0, "/path/to/project/template_esm/campaign/stats")
import stats_vm as sv
d = pd.read_csv("/path/to/project/SI_2026/datasets/Dataset_S11_ESMFold_confidence/esmfold_vs_afdb_v6_per_protein.tsv", sep="\t")
d["tax_group7"] = sv.norm_group(d.tax_group7)
code = d.proteome_tar.str.extract(r"^UP\d+_\d+_([A-Z0-9]+)_v6")[0]
d["plant"] = code.map({"ARATH": "Arabidopsis thaliana", "SOYBN": "Glycine max", "ORYSJ": "Oryza sativa subsp. japonica", "MAIZE": "Zea mays"})
d["plddt_global"] = d.afdb_v6_plddt_global
t = sv.correlations(d, "plddt_global", "esmfold_plddt", "tax_group7", "AF2 vs ESMFold pLDDT")
pd.set_option("display.width", 250)
print(t.to_string())
