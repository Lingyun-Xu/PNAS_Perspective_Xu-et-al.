#!/usr/bin/env python3
"""Alignment depth of the selected proteins and of the plant proteins solved in 2021-2026 (concise PNAS SI, Fig. S19 and
Table S33): AF2-style Neff of the OpenProteinSet alignments, with Neff kept in the tables.

Values are those of the merged AFDB v6 set (Section S4.1; Dataset S10 merged_neff.parquet); nothing is recomputed.
Proteins without a value are left out.
  neff_sets.png        log10(AF2-style Neff + 1) of sel3800, sel690, pdb3825 (full-length) and the two random baselines
  neff_sets.tsv        per-protein values
  neff_sets_tests.tsv  median, share below 10, Mann-Whitney U / Cliff's delta versus pdb3825, for both measures
Usage: python plot_neff_sets.py <outdir> <merged_neff.parquet> <per_protein_comparison.tsv>
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

OUT, NEFF, PC = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
OUT.mkdir(parents=True, exist_ok=True)
COL = "AF2-style_Neff"
SETS = [("sel3800", "AI sel3800", "#d95f02"), ("sel690", "AI sel690", "#e7298a"),
        ("pdb3825_full", "PDB 2021–26", "#1b9e77"), ("random_matched", "Random, sel3800\nspecies mix", "#7570b3"),
        ("random_pool", "Random plant\n(no structure)", "#666666")]


def style(ax):  # as in plot_sel690_ages.py
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(1.2)
    ax.tick_params(width=1.2, length=3, labelsize=8)


def cliff(a, b):
    u = mannwhitneyu(a, b, alternative="two-sided")
    return 2 * u.statistic / (len(a) * len(b)) - 1, u.pvalue


s10 = pd.read_parquet(NEFF, columns=["db_accession", "neff", COL]).drop_duplicates("db_accession")
s10["db_accession"] = s10.db_accession.astype(str)
pc = pd.read_csv(PC, sep="\t", low_memory=False)
pc = pc[pc.set.isin([s for s, _, _ in SETS])].copy()
for c, log in (("neff", "log10_neff"), (COL, f"log10_{COL}")):
    pc[c] = pc.accession.astype(str).map(s10.set_index("db_accession")[c])
    pc[log] = np.log10(pc[c].clip(lower=0) + 1)
pc[["set", "accession", "organism", "neff", "log10_neff", COL, f"log10_{COL}"]].to_csv(OUT / "neff_sets.tsv", sep="\t", index=False)

rows = []
for s, label, _ in SETS:
    row = {"set": s, "n_total": int(pc.set.eq(s).sum())}
    for c, tag in (("neff", "neff"), (COL, "af2style")):
        v = pc.loc[pc.set.eq(s), c].dropna().to_numpy()
        ref = pc.loc[pc.set.eq("pdb3825_full"), c].dropna().to_numpy()
        row[f"n_with_{tag}"] = len(v)
        row[f"median_{tag}"] = float(np.median(v))
        row[f"share_{tag}_lt_10"] = float(np.mean(v < 10))        # on the values themselves: many are exactly 10
        if s != "pdb3825_full":
            d, p = cliff(v, ref)
            row[f"cliffs_delta_{tag}_vs_pdb3825"] = d
            row[f"p_{tag}_vs_pdb3825"] = p
    rows.append(row)
tests = pd.DataFrame(rows)
tests.to_csv(OUT / "neff_sets_tests.tsv", sep="\t", index=False)
print(tests.T.to_string(), flush=True)

fig, ax = plt.subplots(figsize=(7.2, 4.6))
data = [pc.loc[pc.set.eq(s), f"log10_{COL}"].dropna().to_numpy() for s, _, _ in SETS]
bp = ax.boxplot(data, patch_artist=True, showfliers=False, widths=0.6, medianprops={"color": "black"})
for patch, (_, _, color) in zip(bp["boxes"], SETS):
    patch.set_facecolor(color)
    patch.set_alpha(0.75)
ax.set_xticks(range(1, len(SETS) + 1), [f"{label}\nn={len(v):,}" for (_, label, _), v in zip(SETS, data)], fontsize=8)
ax.set_ylabel("log10(AF2-style Neff + 1)")
ax.set_title("Alignment depth of AI-selected and recently solved plant proteins")
style(ax)
fig.tight_layout()
fig.savefig(OUT / "neff_sets.png", bbox_inches="tight", dpi=200)
plt.close(fig)
print("DONE", OUT, flush=True)
