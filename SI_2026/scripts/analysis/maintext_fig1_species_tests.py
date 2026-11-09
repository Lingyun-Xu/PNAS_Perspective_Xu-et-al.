"""Main-text Fig. 1 (AFDB v6), four plant proteomes vs E. coli, S. cerevisiae and H. sapiens: Kruskal-Wallis test
across the four species and Dunn's test (tie-corrected, Holm over the three comparisons with plants) for the
closest-template distance (panels b, d-g), AF2-style Neff per residue (panels c, h-k) and global pLDDT.

Inputs are the SI datasets: S10 merged_neff.parquet (global pLDDT, AF2-style Neff, query length), S17
FigS06_S07/template_distances.parquet (closest-template distance) and S12 merged_per_protein.parquet (proteome of
each accession). The selection of proteins is that of updated_maintext_figure/scripts/maintext_statistics.py.

Usage: python maintext_fig1_species_tests.py <SI_2026/datasets> <out.tsv>
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

D, OUT = Path(sys.argv[1]), Path(sys.argv[2])
SPECIES = {"PLANT": ("_ARATH", "_ORYSJ", "_MAIZE", "_SOYBN"), "ECOLI": ("_ECOLI",), "YEAST": ("_YEAST",),
           "HUMAN": ("_HUMAN",)}
ORDER = ["PLANT", "ECOLI", "YEAST", "HUMAN"]


def holm(p):
    p = np.asarray(p, float)
    adj, running = np.empty_like(p), 0.0
    for rank, i in enumerate(np.argsort(p)):
        running = max(running, (len(p) - rank) * p[i])
        adj[i] = min(1.0, running)
    return adj


def dunn_vs_first(groups):
    """Dunn's test of the first group against each other group, ranks over all groups, tie-corrected"""
    vals = np.concatenate(groups)
    ranks = stats.rankdata(vals)
    n = len(vals)
    _, counts = np.unique(vals, return_counts=True)
    tie = (counts ** 3 - counts).sum() / (12 * (n - 1))
    bounds = np.cumsum([0] + [len(g) for g in groups])
    mean_rank = [ranks[bounds[k]:bounds[k + 1]].mean() for k in range(len(groups))]
    out = []
    for k in range(1, len(groups)):
        se = np.sqrt((n * (n + 1) / 12 - tie) * (1 / len(groups[0]) + 1 / len(groups[k])))
        z = (mean_rank[0] - mean_rank[k]) / se
        out.append((z, 2 * stats.norm.sf(abs(z))))
    return out


m = pd.read_parquet(D / "Dataset_S10_alignment_depth/merged_neff.parquet",
                    columns=["db_accession", "plddt_global", "AF2-style_Neff", "query_length"])
t = pd.read_parquet(D / "Dataset_S17_source_data_SI_figures/FigS06_S07/template_distances.parquet",
                    columns=["db_accession", "n_templates", "d_template_min"])
t = t[(t.n_templates.fillna(0) > 0) & np.isfinite(t.d_template_min)]
p = pd.read_parquet(D / "Dataset_S12_disorder/merged_per_protein.parquet", columns=["db_accession", "proteome"])
p = p.dropna(subset=["proteome"]).drop_duplicates("db_accession")
lab = pd.Series(np.nan, index=p.index, dtype=object)
for sp, suffixes in SPECIES.items():
    lab[p.proteome.str.endswith(suffixes)] = sp
species = p.assign(species=lab).dropna(subset=["species"]).set_index("db_accession").species
m["species"] = m.db_accession.map(species)
t["species"] = t.db_accession.map(species)
ma = m.dropna(subset=["species", "AF2-style_Neff"]).copy()
ma["af2_per_len"] = ma["AF2-style_Neff"] / ma.query_length

measures = {
    "closest-template distance (b, d-g)": {s: t.loc[t.species == s, "d_template_min"].to_numpy(float) for s in ORDER},
    "AF2-style Neff per residue (c, h-k)": {s: ma.loc[ma.species == s, "af2_per_len"].to_numpy(float) for s in ORDER},
    "global pLDDT, all proteins of the four species": {s: m.loc[m.species == s, "plddt_global"].to_numpy(float)
                                                       for s in ORDER},
}
rows = []
for name, by in measures.items():
    groups = [by[s] for s in ORDER]
    h, pk = stats.kruskal(*groups)
    n = sum(len(g) for g in groups)
    dunn = dunn_vs_first(groups)
    padj = holm([d[1] for d in dunn])
    row = {"measure": name, "n": n, **{f"n {s}": len(by[s]) for s in ORDER}, "Kruskal-Wallis H (3 df)": h,
           "p": pk, "epsilon^2": h / (n - 1)}
    for (z, _), pa, s in zip(dunn, padj, ORDER[1:]):
        row[f"Dunn z, plants vs {s}"] = z
        row[f"Dunn p Holm, plants vs {s}"] = pa
    rows.append(row)
res = pd.DataFrame(rows)
OUT.parent.mkdir(parents=True, exist_ok=True)
res.to_csv(OUT, sep="\t", index=False)
pd.set_option("display.width", 250)
print(res.T.to_string())
