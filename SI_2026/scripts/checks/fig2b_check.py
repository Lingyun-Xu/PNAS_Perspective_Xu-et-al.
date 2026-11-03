import pandas as pd
c = pd.read_csv("/path/to/project/SI_2026/datasets/Dataset_S16_source_data_main_text_figures/Fig2/Fig2_plotted_counts.tsv", sep="\t")
b = c[c.panel == "2b"].set_index("group")["count"]
fig = {"Bacteria": 29.6, "Homo": 27.6, "Other Mammalia": 10.3, "Viruses": 9.1, "Others": 6.9, "Fungi": 6.4, "Viridiplantae": 4.6, "Archaea": 2.7, "Other Vertebrata": 1.6, "Insecta": 1.3}
b2 = b.copy()
for k in ("Synthetic", "Nematoda"):
    if k in b2: b2["Others"] = b2.get("Others", 0) + b2.pop(k)
na = int(b2.pop("nan")) if "nan" in b2 else 0
N = 176322
print(f"assigned {int(b2.sum()):,} | unassigned {na:,} | total {int(b2.sum()) + na:,}")
for g, pct in fig.items():
    n = int(b2.get(g, 0)); need = pct / 100 * N
    print(f"  {g:16s} assigned {n:>7,} ({100*n/N:5.2f}% of 176,322) | figure {pct:4.1f}% -> {need:>9,.0f} | missing {need-n:>8,.0f}")
