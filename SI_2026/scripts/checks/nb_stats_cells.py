import json, re
nb = json.load(open("/path/to/home/uniclust30_MSA/Refer_back_to folder_names_and_calculate_MSA_depth.ipynb"))
keys = re.compile(r"jensenshannon|gaussian_kde|kdeplot|kruskal|posthoc|dunn|ks_2samp|spearmanr|pearsonr|linregress|OLS|> *30|linspace|xlim", re.I)
for i, c in enumerate(nb["cells"]):
    if c["cell_type"] != "code":
        continue
    src = "".join(c["source"])
    if not keys.search(src):
        continue
    outs = []
    for o in c.get("outputs", []):
        if "text" in o:
            outs.append("".join(o["text"]))
        elif "data" in o and "text/plain" in o["data"]:
            outs.append("".join(o["data"]["text/plain"]))
    print(f"===== cell {i}\n{src[:2500]}\n--- output:\n{''.join(outs)[:1200]}")
