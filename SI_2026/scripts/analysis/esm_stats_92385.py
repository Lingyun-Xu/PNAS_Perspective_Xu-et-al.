"""ESMFold rows of Tables S24/S25 recomputed on all 92,385 proteins of Dataset S11 (stats_vm.py read only
scores_precise.jsonl, 92,322 proteins, omitting the 63 Atlas scores already stored at full precision). Same functions
as stats_vm.py: Kruskal-Wallis with epsilon^2 = H/(n-1); plants vs each group, Cliff's delta from Mann-Whitney U,
Holm; among the four plant proteomes, Kruskal-Wallis, epsilon^2, Dunn (tie-corrected, Holm) with Cliff's delta."""
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, "/path/to/project/template_esm/campaign/stats")
import stats_vm as sv  # noqa: E402

d = pd.read_csv("/path/to/project/SI_2026/datasets/Dataset_S11_ESMFold_confidence/esmfold_vs_afdb_v6_per_protein.tsv", sep="\t")
d["tax_group7"] = sv.norm_group(d.tax_group7)
code = d.proteome_tar.str.extract(r"^UP\d+_\d+_([A-Z0-9]+)_v6")[0]
d["plant"] = code.map({"ARATH": "Arabidopsis thaliana", "SOYBN": "Glycine max", "ORYSJ": "Oryza sativa subsp. japonica",
                       "MAIZE": "Zea mays"})
d["esmfold"] = d.esmfold_plddt
d["af2_minus"] = d.afdb_v6_plddt_global - d.esmfold_plddt
print("rows", len(d), "| with ESMFold", int(d.esmfold.notna().sum()), "| groups", d.tax_group7.value_counts().to_dict())
print("plant proteome n:", d.plant.value_counts().to_dict())
for col, label in (("esmfold", "ESMFold pLDDT"), ("af2_minus", "AF2 global - ESMFold")):
    g, a = sv.compare(d, col, "tax_group7", label)
    kw = g[g.comparison.eq("Kruskal-Wallis across groups")].iloc[0]
    vs = g[g.comparison.str.startswith("Viridiplantae vs")]
    print(f"\n== {label}: n {int(kw.n):,} eps2 {kw.effect:.3f} plant median {vs.median_plants.iloc[0]:.2f}")
    for _, r in vs.iterrows():
        print(f"   {r.comparison:34s} delta {r.effect:+.3f} p_holm {r.p_holm:.3g}")
    kp = g[g.comparison.eq("Kruskal-Wallis among plant proteomes")].iloc[0]
    print(f"   plants: n {int(kp.n):,} eps2 {kp.effect:.3f} p {kp.p:.3g} medians "
          + ", ".join(f"{k} {kp['median_' + k]:.2f}" for k in ("Arabidopsis", "soybean", "rice", "maize")))
    top = a.reindex(a.cliffs_delta_a_vs_b.abs().sort_values(ascending=False).index).iloc[0]
    print(f"   largest pairwise delta: {top.a} vs {top.b}: {top.cliffs_delta_a_vs_b:+.3f} (Dunn-Holm p {top.p_holm:.3g})")
