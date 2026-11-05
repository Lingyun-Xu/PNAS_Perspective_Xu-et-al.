"""Pfam accession and InterPro entry name of every Pfam family shown in the annotation figure (si_pnas_local.py
annotation, panel D: the most frequent best Pfam families of sel3800, sel690 and the full-length pdb3825 proteins),
with the local hmmscan description where one exists.
Writes SI_2026/pnas_figures/parts/annotation_pfam_names.tsv"""
import json
import urllib.request
from pathlib import Path

import pandas as pd

from si_pnas_local import ANN_SETS, PC, top_families

B = Path(r"/path/to/project")
C = B / "template_esm" / "campaign"
pc = pd.read_csv(PC, sep="\t", low_memory=False)
names = []
for s in ANN_SETS:
    shown = top_families(pc.loc[pc.set.eq(s), "pfam_best_family"].dropna().value_counts())
    print(s, len(shown), "families, at least", shown.min(), "proteins")
    names += [n for n in shown.index if n not in names]
acc, desc = {}, {}
for f in (C / "figures" / "pfam_3800" / "selection_3800_pfam.tsv", C / "figures" / "pfam_690" / "selection_690_pfam.tsv"):
    t = pd.read_csv(f, sep="\t", dtype=str).dropna(subset=["best_name"])
    for r in t.itertuples():
        acc.setdefault(r.best_name, r.best_pfam.split(".")[0])
        desc.setdefault(r.best_name, r.best_desc)
for f in ("pdb3825_pfam_fulllength.tsv", "pdb3825_pfam_constructs.tsv"):
    t = pd.read_csv(C / "pdb3825" / "results" / f, sep="\t", dtype=str).dropna(subset=["pfam_best_family"])
    for r in t.itertuples():
        acc.setdefault(r.pfam_best_family, r.pfam_best_accession.split(".")[0])
rows = []
for n in names:
    with urllib.request.urlopen(f"https://www.ebi.ac.uk/interpro/api/entry/pfam/{acc[n]}", timeout=60) as r:
        meta = json.load(r)["metadata"]
    rows.append({"pfam_name": n, "pfam_accession": acc[n], "interpro_name": meta["name"]["name"],
                 "interpro_short": meta["name"].get("short"), "hmmscan_description": desc.get(n, "")})
out = pd.DataFrame(rows)
out.to_csv(B / "SI_2026" / "pnas_figures" / "parts" / "annotation_pfam_names.tsv", sep="\t", index=False)
print(out[["pfam_name", "pfam_accession", "interpro_name"]].to_string(index=False))
