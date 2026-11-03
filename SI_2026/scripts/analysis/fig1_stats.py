"""Statistics of main-text Fig. 1b and 1h-k, with the code of the 2024 notebook
(Refer_back_to folder_names_and_calculate_MSA_depth.ipynb, cells 39, 40, 44-50):
  original   AFDB v4 global pLDDT, 34,921 chain-level rows (MSA_depth_with_plddt_34921.csv);
  recomputed AFDB v6 global pLDDT and Neff per residue of the revision alignments (Dataset S10), one row per protein of
             the same seven proteomes.
p-values that underflow are computed in log space (correlations: t distribution via exact_p.log10_p; Kruskal-Wallis:
chi-square, 3 df; Dunn: normal tail, Bonferroni x6; KS: first term of the Kolmogorov series).
Also: what produces the vertical stripes of Fig. S5B, and the points outside the plotted ranges of Fig. S5."""
import sys

import numpy as np
import pandas as pd
import scikit_posthocs as sp
from scipy import special, stats
from scipy.spatial.distance import jensenshannon
from scipy.stats import gaussian_kde

sys.path.insert(0, "/path/to/project/SI_2026/scripts")
from exact_p import log10_p  # noqa: E402

B = "/path/to/project"
L10 = np.log(10)


def kde_and_jsd(d1, d2):  # notebook cells 47-50
    grid = np.linspace(min(d1.min(), d2.min()), 160, 1000)
    p1 = gaussian_kde(d1, bw_method="silverman")(grid)
    p2 = gaussian_kde(d2, bw_method="silverman")(grid)
    p1 /= np.trapz(p1, grid)
    p2 /= np.trapz(p2, grid)
    return jensenshannon(p1, p2, base=2)


def log10_chi2_sf3(x):
    z = np.sqrt(x)
    return np.logaddexp(np.log(2) + special.log_ndtr(-z), 0.5 * np.log(2 * x / np.pi) - x / 2) / L10


def dunn_log10(df, val, grp, groups):
    """Dunn's test with tie correction and Bonferroni (6 pairs), as scikit-posthocs, returning log10 p."""
    r = stats.rankdata(df[val])
    n_all = len(r)
    _, t = np.unique(df[val], return_counts=True)
    tie = (t ** 3 - t).sum() / (12 * (n_all - 1))
    mean_rank = {g: r[(df[grp] == g).to_numpy()].mean() for g in groups}
    n = {g: int((df[grp] == g).sum()) for g in groups}
    out = {}
    for i, a in enumerate(groups):
        for b in groups[i + 1:]:
            z = abs(mean_rank[a] - mean_rank[b]) / np.sqrt((n_all * (n_all + 1) / 12 - tie) * (1 / n[a] + 1 / n[b]))
            out[(a, b)] = min(0.0, (np.log(2) + special.log_ndtr(-z)) / L10 + np.log10(6))
    return out


def fmt_p(lp):
    if lp > -4:
        return f"{10 ** lp:.2g}"
    m, e = 10 ** (lp - np.floor(lp)), int(np.floor(lp))
    return f"{m:.1f} × 10^{e}"


def report(label, df, val, plddt, neff):
    print(f"\n==== {label}")
    groups = ["PLANT", "HUMAN", "YEAST", "ECOLI"]
    for g in groups:
        s = df[df.grp == g]
        extra = (" (" + ", ".join(f"{k} {v:,}" for k, v in s["prot"].value_counts().items()) + ")") if g == "PLANT" else ""
        print(f"  n {g}: {len(s):,}{extra} | median Neff/len {s[val].median():.2f} | >30 {100 * (s[val] > 30).mean():.2f}%"
              f" | median pLDDT {s[plddt].median():.1f}")
    n_all = len(df)
    rho = stats.spearmanr(df[val], df[plddt])[0]
    r = stats.pearsonr(df[val], df[plddt])[0]
    lr = stats.linregress(df[val], df[plddt])
    rho_n = stats.spearmanr(df[neff], df[plddt])[0]
    print(f"  n all {n_all:,} | Spearman {rho:.4f} p {fmt_p(log10_p(rho, n_all))} | Pearson {r:.4f} p {fmt_p(log10_p(r, n_all))}"
          f" | R2 {lr.rvalue ** 2:.4f} | Spearman with Neff {rho_n:.4f} p {fmt_p(log10_p(rho_n, n_all))}")
    h, p = stats.kruskal(*[df.loc[df.grp == g, val] for g in groups])
    print(f"  Kruskal-Wallis H {h:.1f} p {p:.3g} (log-space {fmt_p(log10_chi2_sf3(h))})")
    if n_all < 50000:
        d = sp.posthoc_dunn(df, val_col=val, group_col="grp", p_adjust="bonferroni")
        print("  scikit-posthocs Dunn:", {f"{a}-{b}": f"{d.loc[a, b]:.3g}" for a, b in
                                         [("PLANT", "HUMAN"), ("PLANT", "ECOLI"), ("PLANT", "YEAST")]})
    dl = dunn_log10(df, val, "grp", groups)
    print("  Dunn (log space):", {f"{a}-{b}": fmt_p(v) for (a, b), v in dl.items()})
    for g in ("ECOLI", "HUMAN", "YEAST"):
        x, y = df.loc[df.grp == "PLANT", val].to_numpy(), df.loc[df.grp == g, val].to_numpy()
        ks = stats.ks_2samp(x, y)
        neff_n = len(x) * len(y) / (len(x) + len(y))
        lp = np.log10(2) - 2 * neff_n * ks.statistic ** 2 / L10
        print(f"  plant vs {g}: KS {ks.statistic:.4f} p {ks.pvalue:.3g} (asymptotic {fmt_p(lp)}) | JSD {kde_and_jsd(x, y):.4f}")


# ---------------------------------------------------------------- original table
o = pd.read_csv(f"{B}/uniclust30_MSA/MSA_depth_with_plddt_34921.csv")
o["grp"] = o.organism_scientific.map(lambda s: "HUMAN" if s.startswith("Homo") else "ECOLI" if s.startswith("Escherichia")
                                     else "YEAST" if s.startswith("Saccharomyces") else "PLANT")
o["prot"] = o.organism_scientific.str.split().str[:2].str.join(" ")
report("ORIGINAL (AFDB v4, 34,921 rows)", o, "Neff/len", "plddt", "Neff")

# ---------------------------------------------------------------- recomputed, AFDB v6
m = pd.read_parquet(f"{B}/SI_2026/datasets/Dataset_S10_alignment_depth/merged_neff.parquet")
v6 = pd.read_parquet(f"{B}/SI_2026/datasets/Dataset_S6_AFDB_v6_47_proteomes/AF2_48_proteomes_v6.parquet", columns=["db_accession", "proteome_tar"])
v6["code"] = v6.proteome_tar.str.extract(r"^UP\d+_\d+_([A-Z0-9]+)_v6")[0]
code = {"ARATH": "PLANT", "ORYSJ": "PLANT", "MAIZE": "PLANT", "SOYBN": "PLANT", "HUMAN": "HUMAN", "YEAST": "YEAST", "ECOLI": "ECOLI"}
name = {"ARATH": "Arabidopsis", "ORYSJ": "rice", "MAIZE": "maize", "SOYBN": "soybean"}
r7 = m.merge(v6[v6.code.isin(code)].drop_duplicates("db_accession")[["db_accession", "code"]], on="db_accession")
r7 = r7[np.isfinite(r7.neff_per_len)].copy()
r7["grp"], r7["prot"] = r7.code.map(code), r7.code.map(name)
report("RECOMPUTED (AFDB v6, one row per protein)", r7, "neff_per_len", "plddt_global", "neff")

# ---------------------------------------------------------------- Fig. S5: stripes in B and points outside the axes
print("\n==== Fig. S5")
h = m[np.isfinite(m.neff)].copy()
pl = h[h.tax_group7 == "Viridiplantae"]
for lab, d in (("all", h), ("plants", pl)):
    one = np.isclose(d.neff, 1.0, atol=1e-6)
    small = d.neff < 20
    print(f"{lab}: n {len(d):,} | Neff = 1: {one.sum():,} ({100 * one.mean():.1f}%) | Neff < 20 (x < 1.32): {small.sum():,}"
          f" ({100 * small.mean():.1f}%), {d.loc[small, 'neff'].round(4).nunique():,} distinct values"
          f" | Neff >= 20: {d.loc[~small, 'neff'].round(4).nunique():,} distinct values")
top = pl.neff.round(4).value_counts().head(8)
print("plants, most frequent Neff values (x = log10(Neff+1)):", "; ".join(f"{v:g} (x {np.log10(v + 1):.2f}): {c:,}" for v, c in top.items()))
per_aln = h.groupby("rep_id").size()
print(f"proteins per alignment: median {per_aln.median():.0f}, 99th pct {per_aln.quantile(.99):.0f}, max {per_aln.max():,};"
      f" proteins whose alignment is shared by >= 10 proteins: {100 * (h.rep_id.map(per_aln) >= 10).mean():.1f}%")
pl_aln = pl.groupby("rep_id").size()
print(f"plant proteins per alignment: max {pl_aln.max():,}; plant proteins sharing an alignment with >= 9 other plant proteins:"
      f" {100 * (pl.rep_id.map(pl_aln) >= 10).mean():.1f}%; of plant proteins with Neff >= 20: "
      f"{100 * (pl[pl.neff >= 20].rep_id.map(pl_aln) >= 10).mean():.1f}%")
print(f"outside axes: A depth > 20,000: {(m.msa_raw_depth > 20000).sum():,} (max {m.msa_raw_depth.max():,.0f});"
      f" C Neff/len > 60: {(h.neff_per_len > 60).sum():,}; pLDDT < 20: {(m.plddt_global < 20).sum():,}")
