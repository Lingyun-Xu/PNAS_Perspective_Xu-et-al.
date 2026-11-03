import pandas as pd
d = pd.read_parquet("/path/to/project/SI_2026/datasets/Dataset_S16_source_data_main_text_figures/Fig2e/MD_PDB_entries_with_superkingdom.parquet")
print(len(d), list(d.columns))
idc = [c for c in d.columns if "pdb" in c.lower() and "id" in c.lower()] or [d.columns[0]]
print("PDB id column", idc[0], "unique", d[idc[0]].nunique(), "| organisms", d.organism_name.nunique() if "organism_name" in d else "-")
