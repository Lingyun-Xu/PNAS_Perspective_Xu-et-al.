#!/usr/bin/env python3
"""Tables S11-S14, S25 and S26 of the concise PNAS SI with AF2-style Neff added: the tests of stats_vm.py on
log10(AF2-style_Neff + 1) from Dataset S10, by group and plant proteome (organisms from Dataset S17; proteome entries, not
the Swiss-Prot entries of the same species). The new rows go after those of Neff in group_tests.tsv, among_plant_dunn.tsv and correlation_tests.tsv; existing rows are not changed.
First recomputes the Neff rows with the same code and stops unless they reproduce the existing ones.
Also writes af2style_neff_by_group.tsv (n, medians and Spearman rho with global pLDDT, all proteins and by group).
Usage: python af2style_neff_stats.py <SI_2026/datasets> <SI_2026/statistics>"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, "/path/to/project/template_esm/campaign/stats")
import stats_vm as sv  # noqa: E402

D, ST = Path(sys.argv[1]), Path(sys.argv[2])
LABEL, CORR = "log10 AF2-style Neff", "pLDDT vs log10 AF2-style Neff"
AFTER, CORR_AFTER = "log10 Neff (OpenProteinSet)", "pLDDT vs log10 Neff"
m = pd.read_parquet(D / "Dataset_S10_alignment_depth/merged_neff.parquet",
                    columns=["db_accession", "tax_group7", "plddt_global", "source", "neff", "AF2-style_Neff"])
org = pd.read_parquet(D / "Dataset_S17_source_data_SI_figures/FigS06_S07/template_distances.parquet",
                      columns=["db_accession", "organism_scientific"]).drop_duplicates("db_accession")
m["tax_group7"] = sv.norm_group(m.tax_group7.astype(str))
m["organism"] = m.db_accession.map(org.set_index("db_accession").organism_scientific)
m["plant"] = m.organism.where(m.organism.isin(list(sv.PLANTS)) & m.source.eq("proteome"))   # the four proteomes only
m["log10_neff"] = np.log10(m.neff + 1)
m["log10_af2"] = np.log10(m["AF2-style_Neff"] + 1)


def insert_after(old, new, key, after):
    """new rows directly after the block of `after`; earlier rows of the new measure are replaced"""
    old = old[old[key].ne(new[key].iat[0])]
    last = old.index[old[key].eq(after)].max()
    pos = old.index.get_loc(last) + 1
    return pd.concat([old.iloc[:pos], new.reindex(columns=old.columns), old.iloc[pos:]], ignore_index=True)


# check: the Neff rows from the same code and the same file must reproduce the existing ones
gt, dn, ct = (pd.read_csv(ST / n, sep="\t") for n in ("group_tests.tsv", "among_plant_dunn.tsv", "correlation_tests.tsv"))
g0, a0 = sv.compare(m, "log10_neff", "tax_group7", AFTER)
c0 = sv.correlations(m, "log10_neff", "plddt_global", "tax_group7", CORR_AFTER)
ref = gt[gt.measure.eq(AFTER)].set_index("comparison")
worst = max(abs(float(ref.loc[r.comparison, "effect"]) - r.effect) for r in g0.itertuples())
ref_c = ct[ct.measure.eq(CORR_AFTER) & ct.rho.notna()].set_index("subset")
c1 = c0[c0.rho.notna()].set_index("subset")
worst_c = float((ref_c.rho - c1.rho).abs().max())
print(f"Neff rows recomputed: largest difference in effect {worst:.2e}, in rho {worst_c:.2e}; n equal: {(ref_c.n == c1.n).all()}", flush=True)
assert worst < 1e-6 and worst_c < 1e-6, "existing Neff rows not reproduced"

g, a = sv.compare(m, "log10_af2", "tax_group7", LABEL)
c = sv.correlations(m, "log10_af2", "plddt_global", "tax_group7", CORR)
insert_after(gt, g, "measure", AFTER).to_csv(ST / "group_tests.tsv", sep="\t", index=False)
insert_after(dn, a.assign(measure=LABEL), "measure", AFTER).to_csv(ST / "among_plant_dunn.tsv", sep="\t", index=False)
insert_after(ct, c, "measure", CORR_AFTER).to_csv(ST / "correlation_tests.tsv", sep="\t", index=False)

w = m[m["AF2-style_Neff"].notna()]
rows = []
for name, s in [("All", w)] + [(k, w[w.tax_group7.eq(k)]) for k in sv.GROUPS] + \
        [(f"plant: {v}", w[w.plant.eq(k)]) for k, v in sv.PLANTS.items()]:
    rows.append({"group": name, "n": len(s), "median AF2-style Neff": s["AF2-style_Neff"].median(), "median Neff": s.neff.median(),
                 "Spearman rho, pLDDT vs AF2-style Neff": spearmanr(s["AF2-style_Neff"], s.plddt_global).correlation,
                 "Spearman rho, pLDDT vs Neff": spearmanr(s.neff, s.plddt_global).correlation,
                 "Spearman rho, AF2-style Neff vs Neff": spearmanr(s["AF2-style_Neff"], s.neff).correlation,
                 "share AF2-style Neff < 10 (%)": (s["AF2-style_Neff"] < 10).mean() * 100})
pd.DataFrame(rows).to_csv(ST / "af2style_neff_by_group.tsv", sep="\t", index=False)
pd.set_option("display.width", 200)
print(pd.DataFrame(rows).round(3).to_string(index=False), flush=True)
print(g[["comparison", "n", "effect", "median_plants", "median_other", "p_holm"]].round(3).to_string(index=False), flush=True)
print(a[["a", "b", "median_a", "median_b", "cliffs_delta_a_vs_b", "p_holm"]].round(3).to_string(index=False), flush=True)
print(c.round(3).to_string(index=False), flush=True)
print("DONE", flush=True)
