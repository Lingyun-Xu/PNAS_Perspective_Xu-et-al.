#!/usr/bin/env python3
"""Neff/length (AF2-style Neff / query length) statistics for the SI tables in which the whole-alignment Neff and
Neff/length rows were replaced (author request, 13 Nov 2026): Tables S9, S11, S22 and S23 (group and among-plant tests)
and Table S10 (Spearman rho with global pLDDT by group and plant proteome, Fisher z vs plants, Cochran Q).
Same code and inputs as af2style_neff_stats.py (stats_vm.compare / stats_vm.correlations on Dataset S10; plant
proteomes = proteome entries of the four plant species). First recomputes the existing AF2-style Neff rows (and,
with the Dataset S10 of before 13 Nov 2026, the Neff per residue rows of the whole-alignment Neff, whose columns were
then removed) and stops unless they reproduce group_tests.tsv,
among_plant_dunn.tsv and correlation_tests.tsv. Existing statistics files are not changed; the new rows are written to
<out>/af2style_neff_length_{group_tests,among_plant_dunn,correlation_tests}.tsv.
Usage: python af2style_neff_length_stats.py <SI_2026/datasets> <SI_2026/statistics> <out dir>"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.stats import spearmanr

sys.path.insert(0, "/path/to/project/template_esm/campaign/stats")
import stats_vm as sv  # noqa: E402

D, ST, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
OUT.mkdir(parents=True, exist_ok=True)
LABEL, CORR = "Neff/length (AF2-style Neff / query length)", "pLDDT vs Neff/length (AF2-style)"
S10 = D / "Dataset_S10_alignment_depth/merged_neff.parquet"
WHOLE = [c for c in ("neff", "neff_per_len") if c in pq.read_schema(S10).names]   # removed from Dataset S10 on 13 Nov 2026
m = pd.read_parquet(S10, columns=["db_accession", "tax_group7", "plddt_global", "source", "AF2-style_Neff",
                                  "query_length"] + WHOLE)
org = pd.read_parquet(D / "Dataset_S17_source_data_SI_figures/FigS06_S07/template_distances.parquet",
                      columns=["db_accession", "organism_scientific"]).drop_duplicates("db_accession")
m["tax_group7"] = sv.norm_group(m.tax_group7.astype(str))
m["organism"] = m.db_accession.map(org.set_index("db_accession").organism_scientific)
m["plant"] = m.organism.where(m.organism.isin(list(sv.PLANTS)) & m.source.eq("proteome"))   # the four proteomes only
m["log10_af2"] = np.log10(m["AF2-style_Neff"] + 1)
m["af2_per_len"] = m["AF2-style_Neff"] / m.query_length
print("proteins with Neff/length:", int(m.af2_per_len.notna().sum()), "| query_length missing where AF2-style Neff present:",
      int((m["AF2-style_Neff"].notna() & m.query_length.isna()).sum()), flush=True)

# check: the existing rows must be reproduced by the same code on the same file
gt, dn, ct = (pd.read_csv(ST / n, sep="\t") for n in ("group_tests.tsv", "among_plant_dunn.tsv", "correlation_tests.tsv"))
checks = ([("neff_per_len", "Neff per residue")] if "neff_per_len" in WHOLE else []) + [("log10_af2", "log10 AF2-style Neff")]
for col, label in checks:
    g0, a0 = sv.compare(m, col, "tax_group7", label)
    ref = gt[gt.measure.eq(label)].set_index("comparison")
    worst = max(abs(float(ref.loc[r.comparison, "effect"]) - r.effect) for r in g0.itertuples())
    refa = dn[dn.measure.eq(label)].set_index(["a", "b"])
    worst_a = max(abs(float(refa.loc[(r.a, r.b), "cliffs_delta_a_vs_b"]) - r.cliffs_delta_a_vs_b) for r in a0.itertuples())
    print(f"{label}: largest difference from the existing rows, group effects {worst:.2e}, among-plant delta {worst_a:.2e}",
          flush=True)
    assert worst < 1e-6 and worst_a < 1e-6, f"existing rows of {label} not reproduced"
c0 = sv.correlations(m, "log10_af2", "plddt_global", "tax_group7", "pLDDT vs log10 AF2-style Neff")
ref_c = ct[ct.measure.eq("pLDDT vs log10 AF2-style Neff") & ct.rho.notna()].set_index("subset")
worst_c = float((ref_c.rho - c0[c0.rho.notna()].set_index("subset").rho).abs().max())
print(f"AF2-style Neff correlations: largest rho difference {worst_c:.2e}", flush=True)
assert worst_c < 1e-6, "existing AF2-style Neff correlations not reproduced"

g, a = sv.compare(m, "af2_per_len", "tax_group7", LABEL)
a = a.assign(measure=LABEL)
c = sv.correlations(m, "af2_per_len", "plddt_global", "tax_group7", CORR)
g.to_csv(OUT / "af2style_neff_length_group_tests.tsv", sep="\t", index=False)
a.to_csv(OUT / "af2style_neff_length_among_plant_dunn.tsv", sep="\t", index=False)
c.to_csv(OUT / "af2style_neff_length_correlation_tests.tsv", sep="\t", index=False)
w = m[m.af2_per_len.notna()]
print("all proteins: n", len(w), "| Spearman rho with global pLDDT", round(spearmanr(w.af2_per_len, w.plddt_global).correlation, 4),
      flush=True)
pd.set_option("display.width", 220)
print(g[["comparison", "n", "effect", "median_plants", "median_other", "p", "p_holm"]].to_string(index=False), flush=True)
print(g.filter(like="median_").iloc[[1]].to_string(index=False), flush=True)
print(a[["a", "b", "median_a", "median_b", "dunn_z", "cliffs_delta_a_vs_b", "p_holm"]].to_string(index=False), flush=True)
print(c.to_string(index=False), flush=True)
print("DONE", flush=True)
