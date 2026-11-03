import pandas as pd
p = "/path/to/project/SI_2026/datasets/Dataset_S16_source_data_main_text_figures/Fig2/Fig2_plotted_counts.tsv"
c = pd.read_csv(p, sep="\t", keep_default_na=False)
c["group"] = c.group.replace({"others": "Others", "other mammalia": "Other Mammalia"})
c = c.groupby(["panel", "unit", "group"], sort=False, as_index=False)["count"].sum()
c["value"] = c["count"].map(lambda v: f"{int(v)}" if float(v).is_integer() and "percentage" not in "" else f"{v:g}")
c.loc[c.unit.str.startswith("percentage"), "value"] = c.loc[c.unit.str.startswith("percentage"), "count"].map(lambda v: f"{v:.1f}")
c[["panel", "unit", "group", "value"]].to_csv(p, sep="\t", index=False)
print(c[c.panel == "2b"][["group", "value"]].to_string(index=False))
