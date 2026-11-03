import pandas as pd
c = pd.read_csv("/path/to/project/template_esm/campaign/figures/search_density_change_hhprob50/foldseek_cohort.tsv", sep="\t")
base = {"plus3800": "pdbonly", "plus690": "pdbonly", "pdbonly": "minus3825"}
print("foldseek rows", len(c), "| all four non-null:", int(c[["d_pdbonly", "d_minus3825", "d_plus3800", "d_plus690"]].notna().all(1).sum()))
for after, before in base.items():
    print(f"  {before} -> {after}: both hits {int(c[[f'd_{before}', f'd_{after}']].notna().all(1).sum()):,}")
h = pd.read_csv("/path/to/project/template_esm/campaign/figures/search_density_change_hhprob50/hhsearch_cohort.tsv", sep="\t")
print("hhsearch rows", len(h), "| all four non-null:", int(h[[x for x in h.columns if x.startswith("d_") and x != "d_homolog"]].notna().all(1).sum()))
