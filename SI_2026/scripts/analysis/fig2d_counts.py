"""Main-text Fig. 2d (MaveDB score sets by group): rebuild the plotted counts from the organism table and the group
lists of Fig2_Pie_chart.ipynb cell 58 (list literals only, read with ast), and compare with the percentages printed by
cell 62. Writes nothing unless --write <out.tsv> is given."""
import ast
import json
import sys

import pandas as pd

NB = "/path/to/home/GEO_dataset/Fig2_Pie_chart.ipynb"
nb = json.load(open(NB))
tree = ast.parse("".join(nb["cells"][58]["source"]))
lists = {}
for node in tree.body:
    if isinstance(node, ast.Assign) and isinstance(node.value, ast.List):
        lists[node.targets[0].id] = set(ast.literal_eval(node.value))
LABEL = {"Other_Mammalia": "Other Mammalia", "Other_Vertebrata": "Other Vertebrata", "synthetic": "Synthetic"}
ORDERS = {
    "commented loop": ["Bacteria", "Fungi", "Insecta", "Other_Mammalia", "Homo", "Nematoda", "Viridiplantae",
                       "Other_Vertebrata", "Archaea", "Viruses", "Others", "synthetic"],
    "meltome loop": ["Bacteria", "Fungi", "Others", "Other_Mammalia", "Homo", "Nematoda", "Insecta", "Viridiplantae",
                     "Other_Vertebrata", "Archaea", "Viruses", "synthetic"],
}
mave = pd.read_csv("/path/to/home/GEO_dataset/mavedb-organisms_species_and_numbers.csv")
results = {}
for name, order in ORDERS.items():
    def group(s):
        for k in order:
            if s in lists[k]:
                return LABEL.get(k, k)
        return None
    m = mave.assign(group=mave["Species"].map(group))
    c = m.groupby("group")["Associated Score Sets"].sum().sort_values(ascending=False)
    results[name] = (m, c)
    print(f"== {name}: unassigned organisms {int(m.group.isna().sum())}")
    print("   " + "; ".join(f"{g} {n} ({100 * n / c.sum():.2f}%)" for g, n in c.items()))
printed = "".join("".join(o.get("text", "")) for o in nb["cells"][62].get("outputs", []))
print("== printed by cell 62:\n" + printed[:700])
if len(sys.argv) > 2 and sys.argv[1] == "--write":
    m, c = results["commented loop"]
    m.rename(columns={"group": "pie_group"}).to_csv(sys.argv[2], sep="\t", index=False)
    print("wrote", sys.argv[2])
