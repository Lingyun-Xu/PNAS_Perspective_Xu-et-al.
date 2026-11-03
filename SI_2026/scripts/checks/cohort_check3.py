import pandas as pd
c = pd.read_csv("/path/to/project/template_esm/campaign/figures/search_density_change_hhprob50/foldseek_cohort.tsv", sep="\t")
d = c[["d_pdbonly", "d_minus3825", "d_plus3800", "d_plus690"]]
print("rows", len(c), "| no hit in any search", int(d.isna().all(1).sum()), "| hit in >= 1", int(d.notna().any(1).sum()), "| all four", int(d.notna().all(1).sum()))
