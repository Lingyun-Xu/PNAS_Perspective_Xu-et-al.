import pandas as pd
p = "/path/to/project/SI_2026/datasets/Dataset_S16_source_data_main_text_figures/Fig2/Fig2_plotted_counts.tsv"
c = pd.read_csv(p, sep="\t", keep_default_na=False)
c.loc[(c.panel == "2b") & (c.group.isin(["nan", ""])), "group"] = "hand-assigned, mapping not saved"
e = c[c.panel == "2e"].copy()
others = int(e.loc[e.group.isin(["Others", "Nematoda"]), "count"].sum())
c = c[~((c.panel == "2e") & c.group.isin(["Others", "Nematoda"]))]
c = pd.concat([c, pd.DataFrame([{"panel": "2e", "unit": e.unit.iloc[0], "group": "Others (incl. Nematoda 36)", "count": others}])])
fig = {"Bacteria": 29.6, "H. sapiens": 27.6, "Other Mammalia": 10.3, "Viruses": 9.1, "Others": 6.9, "Fungi": 6.4,
       "Viridiplantae": 4.6, "Archaea": 2.7, "Other Vertebrata": 1.6, "Insecta": 1.3}
c = pd.concat([c, pd.DataFrame([{"panel": "2b", "unit": "percentage printed on the published pie (of 176,322 sequences)", "group": g, "count": v} for g, v in fig.items()])])
c.to_csv(p, sep="\t", index=False)
print(c.groupby("panel").size().to_dict()); print(c[c.panel.isin(["2b", "2e"])].to_string(index=False)[:2500])
