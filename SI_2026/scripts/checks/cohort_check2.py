import pandas as pd
c = pd.read_csv("/path/to/project/template_esm/campaign/figures/search_density_change_hhprob50/foldseek_cohort.tsv", sep="\t")
for ai in ("plus3800", "plus690"):
    cols = ["d_minus3825", "d_pdbonly", f"d_{ai}"]
    print(ai, "vs pdb3825: hits in", cols, "->", int(c[cols].notna().all(1).sum()))
