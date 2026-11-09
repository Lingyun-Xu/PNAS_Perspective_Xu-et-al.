#!/usr/bin/env python3
"""Why Neff, AF2-style Neff and Neff per residue have two peaks in the revision (SI Section S4.1), and why the original
Fig. 1h-k had one. KDEs as in the figures (Silverman bandwidth, 1,000 points). Peaks: local maxima of the curve, padded
so that a maximum at the left edge counts, with a prominence of at least 3% of the highest point (1% for Neff per
residue, whose upper peak is lower and broader); the value in parentheses is the height relative to the highest peak.

Tables written to <out>:
  by_assignment.tsv     peaks by how the protein received its alignment (and without the MMseqs2-assigned proteins)
  by_alignment.tsv      one row per alignment instead of one row per protein
  rows_per_alignment.tsv, rows_vs_depth.tsv   rows per alignment around 1,000, for all 341,470 alignments
  saturation.tsv, saturation_correlation.tsv   AF2-style Neff and Neff against the AFDB MSA depth of the same proteins
  homolog_class.tsv     homolog class (Dataset S13) of the proteins below 100, 100-800 and at least 800
  length_bands.tsv      peaks by query length: the upper peak of AF2-style Neff stays, that of Neff per residue moves
  original_data.tsv     the data of the original Fig. 1h-k (Dataset S16), PDB-chain folders and UniClust30 alignments
  same_proteins.tsv     the proteins of the original analysis with the revision alignments, against all proteins
  same_chains.tsv, same_chains_rows.tsv   the three OpenProteinSet alignments of the same PDB chains (from two_peaks_server.py)
Figure: two_peaks_diagnostics.png
Usage: python neff_two_peaks.py <SI_2026/datasets> <folder with the two_peaks_server.py tables> <out>"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
from scipy.stats import gaussian_kde, spearmanr

D, SERVER, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
OUT.mkdir(parents=True, exist_ok=True)
COL = "AF2-style_Neff"
SPECIES = {"UP000006548_3702_ARATH": "PLANT", "UP000007305_4577_MAIZE": "PLANT", "UP000008827_3847_SOYBN": "PLANT",
           "UP000059680_39947_ORYSJ": "PLANT", "UP000000625_83333_ECOLI": "ECOLI", "UP000002311_559292_YEAST": "YEAST",
           "UP000005640_9606_HUMAN": "HUMAN"}
FOUR = ("PLANT", "ECOLI", "YEAST", "HUMAN")
MEASURES = ((COL, 4000), ("neff", 6000), ("neff_per_len", 30), ("msa_raw_depth", 20000))
COLOR = {"PLANT": "#0EB519", "ECOLI": "#3D5BA9", "YEAST": "#F68A21", "HUMAN": "#E41E26"}


def curve(v, xmax, xmin=0, n=1000):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    grid = np.linspace(xmin, xmax, n)
    return grid, gaussian_kde(v, bw_method="silverman")(grid)


def peaks(v, xmax, prom=None, xmin=0):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    if len(v) < 30 or v.std() == 0:
        return "n < 30"
    grid, y = curve(v, xmax, xmin)
    pk, _ = find_peaks(np.r_[0, y, 0], prominence=y.max() * (prom if prom is not None else 0.03))
    pk = pk - 1
    top = y[pk].max()
    return "; ".join(f"{grid[p]:.4g} ({y[p] / top:.2f})" for p in pk)


def prom_of(measure):
    return 0.01 if measure == "neff_per_len" else 0.03


m = pd.read_parquet(D / "Dataset_S10_alignment_depth/merged_neff.parquet")
acc = pd.read_parquet(D / "Dataset_S12_disorder/merged_per_protein.parquet", columns=["db_accession", "proteome"]).dropna().drop_duplicates("db_accession")
acc["species"] = acc.proteome.map(SPECIES)
s = m[m.source.eq("proteome")].merge(acc.dropna(subset=["species"])[["db_accession", "species"]], on="db_accession")
s = s[s.neff.notna()]
print("proteins of the four species", len(s), s.species.value_counts().to_dict(), flush=True)

# ---------------------------------------------------------------- assignment
rows = []
for col, xmax in MEASURES:
    for sp in FOUR:
        for how in ("all", "is_representative", "cluster_member", "mmseqs", "without mmseqs"):
            t = s[s.species.eq(sp)]
            t = t[t.how.ne("mmseqs")] if how == "without mmseqs" else t if how == "all" else t[t.how.eq(how)]
            rows.append(dict(measure=col, species=sp, assignment=how, n=len(t), share_percent=round(len(t) / s.species.eq(sp).sum() * 100, 1),
                             median=round(t[col].median(), 2), Q1=round(t[col].quantile(.25), 2), Q3=round(t[col].quantile(.75), 2),
                             peaks=peaks(t[col], xmax, prom_of(col))))
pd.DataFrame(rows).to_csv(OUT / "by_assignment.tsv", sep="\t", index=False)

rows = []
for name, t in (("plants, four proteomes", s[s.species.eq("PLANT")]), ("human", s[s.species.eq("HUMAN")]), ("E. coli", s[s.species.eq("ECOLI")]), ("all 988,021 proteins", m[m.neff.notna()])):
    u = t.drop_duplicates("rep_id")
    for col, xmax in MEASURES[:2]:
        rows.append(dict(proteins=name, measure=col, alignments=len(u), peaks=peaks(u[col], xmax)))
pd.DataFrame(rows).to_csv(OUT / "by_alignment.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- rows
a = m[m.neff.notna()].drop_duplicates("rep_id")
h = pd.cut(a.n_sequences, np.arange(800, 1301, 25), right=False).value_counts(sort=False).rename_axis("rows").reset_index(name="alignments")
h.to_csv(OUT / "rows_per_alignment.tsv", sep="\t", index=False)
rows = []
for lo, hi in ((1, 500), (500, 1000), (1000, 1500), (1500, 2000), (2000, 3000), (3000, 10 ** 9)):
    t = a[(a.n_sequences >= lo) & (a.n_sequences < hi)]
    rows.append(dict(rows_from=lo, rows_below=hi, alignments=len(t), median_AF2_style=round(t[COL].median(), 1), median_Neff=round(t.neff.median(), 1)))
pd.DataFrame(rows).to_csv(OUT / "rows_vs_depth.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- saturation against the AFDB MSA depth
b = [0, 10, 30, 100, 300, 500, 1000, 1500, 2000, 3000, 5000, 7500, 10000, 15000, 19000, 25000]
g = s.assign(bin=pd.cut(s.msa_raw_depth, b)).groupby("bin", observed=True)
sat = pd.DataFrame({"proteins": g.size(), "median_AF2_style": g[COL].median().round(0), "Q1_AF2_style": g[COL].quantile(.25).round(0), "Q3_AF2_style": g[COL].quantile(.75).round(0),
                    "median_Neff": g.neff.median().round(0), "median_rows": g.n_sequences.median().round(0)}).reset_index()
sat.to_csv(OUT / "saturation.tsv", sep="\t", index=False)
rows = []
for lo, hi in ((0, 800), (800, 2000)):
    t = s[(s[COL] >= lo) & (s[COL] < hi) & s.msa_raw_depth.notna()]
    rows.append(dict(AF2_style_Neff_from=lo, AF2_style_Neff_below=hi, proteins=len(t), spearman_with_AFDB_MSA_depth=round(spearmanr(t[COL], t.msa_raw_depth).correlation, 3)))
pd.DataFrame(rows).to_csv(OUT / "saturation_correlation.tsv", sep="\t", index=False)
print(rows, flush=True)

# ---------------------------------------------------------------- homolog class
p = pd.read_parquet(D / "Dataset_S13_homolog_class/panel_strata.parquet", columns=["db_accession", "panel_bin"]).drop_duplicates("db_accession")
c = s.merge(p, on="db_accession")
c["band"] = pd.cut(c[COL], [-1, 100, 800, 1e9], right=False, labels=["below 100", "100 to 800", "800 or more"])
rows = []
for sp in FOUR:
    t = c[c.species.eq(sp)]
    for band, u in t.groupby("band", observed=True):
        row = dict(species=sp, AF2_style_Neff=band, proteins=len(u), share_of_species_percent=round(len(u) / len(t) * 100, 1))
        row.update({k: round((u.panel_bin == k).mean() * 100, 1) for k in ("universal", "eukaryote_wide", "kingdom_only", "family_or_genus", "species_only")})
        rows.append(row)
    for k, u in t.groupby("panel_bin"):
        rows.append(dict(species=sp, AF2_style_Neff=f"class {k}", proteins=len(u), share_of_species_percent=round(len(u) / len(t) * 100, 1), median=round(u[COL].median(), 1), peaks=peaks(u[COL], 4000)))
pd.DataFrame(rows).to_csv(OUT / "homolog_class.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- length
rows = []
for sp in ("PLANT", "HUMAN"):
    t = s[s.species.eq(sp)]
    for lo, hi in ((0, 200), (200, 400), (400, 800), (800, 40000)):
        u = t[(t.query_length >= lo) & (t.query_length < hi)]
        rows.append(dict(species=sp, query_length_from=lo, query_length_below=hi, proteins=len(u), peaks_AF2_style=peaks(u[COL], 4000), peaks_Neff_per_residue=peaks(u.neff_per_len, 30, 0.01)))
pd.DataFrame(rows).to_csv(OUT / "length_bands.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- the original data
o = pd.read_csv(D / "Dataset_S16_source_data_main_text_figures/Fig1b_Fig1h-k/MSA_depth_with_plddt_34921.csv")
o["species"] = "other"
for key, name in (("Homo sapiens", "HUMAN"), ("Escherichia coli", "ECOLI"), ("Saccharomyces", "YEAST")):
    o.loc[o.organism_scientific.str.contains(key, na=False), "species"] = name
o.loc[o.organism_scientific.str.contains("Arabidopsis|Oryza|Glycine|Zea", na=False), "species"] = "PLANT"
o["alignments"] = np.where(o.Source.astype(str).str.contains("_"), "PDB chain, three pooled", "UniClust30, single")
rows = []
for sp in FOUR:
    t = o[o.species.eq(sp)]
    for name, u in (("all rows", t), ("PDB chain, three pooled", t[t.alignments.eq("PDB chain, three pooled")]), ("UniClust30, single", t[t.alignments.eq("UniClust30, single")])):
        rows.append(dict(species=sp, rows_of=name, n=len(u), median_Neff=round(u.Neff.median(), 1) if len(u) else None, median_Neff_per_len=round(u["Neff/len"].median(), 2) if len(u) else None,
                         median_total_sequences=u["Total sequences"].median() if len(u) else None,
                         peaks_Neff_per_len=peaks(u["Neff/len"], 160) if len(u) else "", peaks_Neff=peaks(u.Neff, 25000) if len(u) else ""))
pd.DataFrame(rows).to_csv(OUT / "original_data.tsv", sep="\t", index=False)
near = pd.cut(o["Total sequences"], [800, 900, 1000, 1100, 1200, 1300], right=False).value_counts(sort=False)
print("original data, rows near 1,000:", near.to_dict(), flush=True)

# ---------------------------------------------------------------- same proteins
r = s[s.db_accession.isin(set(o.db_accession))]
rows = []
for col, xmax in MEASURES:
    for sp in FOUR:
        t, u = r[r.species.eq(sp)], s[s.species.eq(sp)]
        rows.append(dict(measure=col, species=sp, n_original_proteins=len(t), n_all_proteins=len(u), median_original=round(t[col].median(), 2), median_all=round(u[col].median(), 2),
                         peaks_original=peaks(t[col], xmax, prom_of(col)), peaks_all=peaks(u[col], xmax, prom_of(col)),
                         percent_below_100_original=round((t[col] < 100).mean() * 100, 1) if col != "neff_per_len" else None,
                         percent_below_100_all=round((u[col] < 100).mean() * 100, 1) if col != "neff_per_len" else None))
pd.DataFrame(rows).to_csv(OUT / "same_proteins.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- same chains
d = pd.read_csv(SERVER / "pdb_af2style.tsv", sep="\t")
q = d.pivot(index="chain", columns="alignment", values="af2style_neff").dropna()
print("chains with all three alignments:", len(q), flush=True)
LAB = {"bfd_uniclust_hits": "HHblits BFD + UniClust30", "uniref90_hits": "JackHMMER UniRef90", "mgnify_hits": "JackHMMER MGnify"}
rows = []
for k in LAB:
    v = q[k]
    rows.append(dict(alignment=LAB[k], chains=len(v), median=round(v.median(), 0), Q1=round(v.quantile(.25), 0), Q3=round(v.quantile(.75), 0), P90=round(v.quantile(.9), 0), maximum=round(v.max(), 0),
                     percent_800_to_1600=round(((v >= 800) & (v < 1600)).mean() * 100, 1), percent_above_2000=round((v > 2000).mean() * 100, 1), percent_above_3000=round((v > 3000).mean() * 100, 1),
                     peaks=peaks(v, 8000)))
pd.DataFrame(rows).to_csv(OUT / "same_chains.tsv", sep="\t", index=False)
top = q[q.bfd_uniclust_hits >= 800]
print(f"chains with HHblits AF2-style Neff >= 800: {len(top)}; HHblits IQR {top.bfd_uniclust_hits.quantile(.25):.0f}-{top.bfd_uniclust_hits.quantile(.75):.0f}; "
      f"UniRef90 IQR {top.uniref90_hits.quantile(.25):.0f}-{top.uniref90_hits.quantile(.75):.0f}; Spearman {spearmanr(top.bfd_uniclust_hits, top.uniref90_hits).correlation:.3f}", flush=True)
cr = pd.read_csv(SERVER / "pdb_rows.tsv", sep="\t")
edges = [0, 100, 500, 800, 900, 1000, 1100, 1200, 1500, 2000, 3000, 4999, 5000, 5001, 9999, 10000, 10001, 10 ** 9]
tab = pd.DataFrame({LAB[k]: pd.cut(cr[k][cr[k] >= 0], edges, right=False).value_counts(sort=False) for k in LAB}).rename_axis("rows").reset_index()
tab.to_csv(OUT / "same_chains_rows.tsv", sep="\t", index=False)
print(tab.to_string(index=False), flush=True)
cov = pd.read_csv(SERVER / "uniclust_sample_coverage.tsv", sep="\t")
deep = cov[cov.rows >= 1000]
print(f"coverage sample {len(cov):,}: alignments with >= 1,000 rows {len(deep):,}; median column with >= 1,000 rows in {(deep.cov_median >= 1000).mean() * 100:.1f}% of them; "
      f"alignments with < 1,000 rows: {(cov[cov.rows < 1000].cov_median >= 1000).sum()}", flush=True)
print("median-column coverage, counts:", pd.cut(cov.cov_median, [0, 200, 400, 600, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 2000, 3000, 10 ** 6], right=False).value_counts(sort=False).to_dict(), flush=True)

# ---------------------------------------------------------------- figure
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False})
fig, ax = plt.subplots(2, 3, figsize=(10.5, 6.2))
ax = ax.ravel()
x = np.arange(800, 1300, 10)
hh, _ = np.histogram(a.n_sequences, bins=np.r_[x, 1300])
ax[0].bar(x, hh, width=9, align="edge", color="#555555")
ax[0].axvline(1000, color="black", ls="--", lw=0.8)
ax[0].set(xlabel="Rows in the alignment", ylabel="Alignments (of 341,470), per 10 rows", title="A  Step at 1,000 rows: the HHblits filter (-diff 1000)")
for sp in FOUR:
    gr, y = curve(s.loc[s.species.eq(sp), COL], 4000)
    ax[1].plot(gr, y, color=COLOR[sp], lw=1.4, label={"PLANT": "Plants", "ECOLI": "E. coli", "YEAST": "Yeast", "HUMAN": "Human"}[sp])
ax[1].legend(frameon=False, loc="upper right")
ax[1].set(xlabel="AF2-style Neff", ylabel="Density", title="B  The two peaks of Fig. 1h-k")
for how, lab, col in (("is_representative", "Own alignment", "#7570b3"), ("cluster_member", "Cluster member", "#0EB519"), ("mmseqs", "MMseqs2 homolog", "#d95f02")):
    t = s[s.species.eq("PLANT") & s.how.eq(how)]
    gr, y = curve(t[COL], 4000)
    ax[2].plot(gr, y, color=col, lw=1.4, label=f"{lab} (n = {len(t):,})")
ax[2].set_yscale("log")
ax[2].set_ylim(1e-5, 6e-3)
ax[2].legend(frameon=False, loc="upper right")
ax[2].set(xlabel="AF2-style Neff", ylabel="Density (log scale)", title="C  Plants: two peaks in each class")
for k, col in (("bfd_uniclust_hits", "#111111"), ("uniref90_hits", "#d95f02"), ("mgnify_hits", "#7570b3")):
    gr, y = curve(q[k], 8000)
    ax[3].plot(gr, y, color=col, lw=1.4, label=LAB[k])
ax[3].legend(frameon=False, loc="upper right")
ax[3].set(xlabel="AF2-style Neff", ylabel="Density", title=f"D  The same {len(q)} PDB chains, one alignment at a time")
g2 = s.assign(bin=pd.cut(s.msa_raw_depth, np.arange(0, 21000, 1000))).groupby("bin", observed=True)
xs = [iv.mid for iv in g2[COL].median().index]
ax[4].plot(xs, g2[COL].median().values, color="black", lw=1.4, label="Median AF2-style Neff")
ax[4].fill_between(xs, g2[COL].quantile(.25).values, g2[COL].quantile(.75).values, color="#bbbbbb", alpha=0.6, lw=0, label="Interquartile range")
ax[4].legend(frameon=False, loc="lower right")
ax[4].set(xlabel="AFDB MSA depth of the same proteins (sequences)", ylabel="AF2-style Neff", title="E  Levels off against the AFDB MSA depth")
t = s[s.species.eq("PLANT")]
gr, y = curve(t[COL], 4000)
ax[5].plot(gr, y, color="#0EB519", lw=1.4, label=f"All (n = {len(t):,})")
t = r[r.species.eq("PLANT")]
gr, y = curve(t[COL], 4000)
ax[5].plot(gr, y, color="black", lw=1.4, label=f"In the original analysis (n = {len(t):,})")
ax[5].legend(frameon=False, loc="upper right")
ax[5].set(xlabel="AF2-style Neff", ylabel="Density", title="F  Plants: original proteins versus all")
for t in ax:
    t.title.set_fontsize(8.5)
    t.title.set_ha("left")
    t.title.set_x(0)
fig.tight_layout(w_pad=2.0, h_pad=2.0)
fig.savefig(OUT / "two_peaks_diagnostics.png", dpi=250)
print("DONE", flush=True)
