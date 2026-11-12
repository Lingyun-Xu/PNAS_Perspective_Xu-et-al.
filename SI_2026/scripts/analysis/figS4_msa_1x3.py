#!/usr/bin/env python3
"""Fig. S4 of the SI for the rebuttal, 1 x 3 (author request, 9 Nov 2026): global pLDDT versus (A) the AFDB MSA depth,
(B) log10(AF2-style Neff + 1) and (C) AF2-style Neff per residue (AF2-style Neff divided by the length of the
alignment's query, the quantity labelled Neff/length in main-text Fig. 1C), for the merged AFDB v6 set (47 proteomes
and Swiss-Prot). Replaces the 2 x 2 version (s5() of si_pnas_server.py: MSA depth, Neff, AF2-style Neff, Neff per
residue); same plotting code (colours, opacities, drawing order, marker size, least-squares fit with 95% band), drawn
from Dataset S10 alone.
Usage: python figS4_msa_1x3.py <SI_2026/datasets> <out.png>"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

D, OUT = Path(sys.argv[1]), Path(sys.argv[2])
plt.rcParams.update({"font.family": "DejaVu Sans", "axes.titleweight": "normal"})
GROUP7 = ["Bacteria", "Archaea", "Animalia", "Fungi", "Protist", "Parasitic Worms", "Viridiplantae"]
PAL = {"Viridiplantae": "#0eb519", "Animalia": "#E41E26", "Fungi": "#F68A21", "Bacteria": "#3D5BA9",
       "Archaea": "#7f7f7f", "Protist": "#9467bd", "Parasitic Worms": "#8c564b"}
ALPHA = {"Bacteria": 0.08, "Animalia": 0.08, "Archaea": 0.2, "Fungi": 0.12, "Protist": 0.2,
         "Parasitic Worms": 0.15, "Viridiplantae": 0.35}


def frame(ax, lw, tick_len, labelsize):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(lw)
    ax.tick_params(width=lw, length=tick_len, labelsize=labelsize)


def scatter_fit(ax, df, xcol, xlabel, xlim, xticks, title, ylabel=True):
    rng = np.random.default_rng(0)
    for g in GROUP7:                                  # plants last, each group in random order
        sub = df[df.tax_group7 == g]
        sub = sub.iloc[rng.permutation(len(sub))]
        ax.scatter(sub[xcol], sub.plddt_global, s=0.3, c=PAL[g], alpha=ALPHA[g], linewidths=0, rasterized=True)
    x, y = df[xcol].to_numpy(float), df.plddt_global.to_numpy(float)
    grid = np.linspace(xlim[0], xlim[1], 200)
    res = stats.linregress(x, y)
    resid = y - (res.intercept + res.slope * x)
    s = np.sqrt((resid ** 2).sum() / (len(x) - 2))
    se = s * np.sqrt(1 / len(x) + (grid - x.mean()) ** 2 / ((x - x.mean()) ** 2).sum())
    fit = res.intercept + res.slope * grid
    ax.fill_between(grid, fit - 1.96 * se, fit + 1.96 * se, color="#434343", alpha=0.25, lw=0)
    ax.plot(grid, fit, color="#434343", lw=1.5)
    rho, _ = stats.spearmanr(x, y)
    r, _ = stats.pearsonr(x, y)
    ax.set_xlim(*xlim)
    ax.set_ylim(20, 100)
    ax.set_xticks(xticks)
    ax.set_yticks([20, 40, 60, 80, 100])
    frame(ax, 1.0, 2.5, 7)
    ax.set_xlabel(xlabel, fontsize=8)
    if ylabel:
        ax.set_ylabel("Global pLDDT (AFDB v6)", fontsize=8)
    ax.set_title(f"{title}\nρ = {rho:.3f}, r = {r:.3f}", fontsize=7)    # Spearman rho, Pearson r (legend)
    outside = int((x > xlim[1]).sum())
    print(f"{title}: n {len(x):,}; rho {rho:.4f}; r {r:.4f}; slope {res.slope:.4g}; intercept {res.intercept:.4g}; "
          f"beyond x = {xlim[1]}: {outside}", flush=True)


m = pd.read_parquet(D / "Dataset_S10_alignment_depth/merged_neff.parquet",
                    columns=["tax_group7", "plddt_global", "msa_raw_depth", "AF2-style_Neff", "query_length"])
depth = m[m.msa_raw_depth.notna()].copy()
af2 = m[m["AF2-style_Neff"].notna()].copy()
af2["log10_af2style"] = np.log10(af2["AF2-style_Neff"] + 1)
af2["af2style_per_len"] = af2["AF2-style_Neff"] / af2.query_length
xmax_b = float(af2.log10_af2style.max()) + 0.1

fig, axes = plt.subplots(1, 3, figsize=(6.0, 2.6))
scatter_fit(axes[0], depth, "msa_raw_depth", "MSA depth (sequences)", (0, 20000), [0, 5000, 10000, 15000, 20000],
            "AFDB MSA depth")
scatter_fit(axes[1], af2, "log10_af2style", "log10(AF2-style Neff + 1)", (0, xmax_b), [0, 1, 2, 3],
            "AF2-style Neff (OpenProteinSet)", ylabel=False)
scatter_fit(axes[2], af2, "af2style_per_len", "AF2-style Neff / length", (0, 40), [0, 10, 20, 30, 40],
            "AF2-style Neff per residue", ylabel=False)
# B: AF2-style Neff = 30 sequences, the depth below which AlphaFold2 accuracy drops (author request, 9 Nov 2026), and
# the share of each group below it
x30 = np.log10(30 + 1)
axes[1].axvline(x30, color="black", ls="--", lw=0.9, zorder=4)
axes[1].text(x30 + 0.06, 22, "Neff = 30", rotation=90, fontsize=5.5, va="bottom", ha="left", zorder=5,
             bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=0.6))
below = af2.assign(b=af2["AF2-style_Neff"] < 30).groupby("tax_group7").b.agg(["mean", "size"])
print("AF2-style Neff < 30, by group:", {g: (round(below.loc[g, "mean"] * 100, 2), int(below.loc[g, "size"]))
                                         for g in GROUP7},
      "| all:", round(float((af2["AF2-style_Neff"] < 30).mean() * 100), 2), f"(x of the line {x30:.3f})", flush=True)
axes[0].set_xticks([0, 5000, 10000, 15000, 20000])
axes[0].set_xticklabels(["0", "5k", "10k", "15k", "20k"])
for a, t in zip(axes, "ABC"):
    a.text(-0.26 if t == "A" else -0.16, 1.17, t, transform=a.transAxes, fontsize=9, fontweight="bold", va="bottom",
           ha="left")
handles = [plt.Line2D([0], [0], marker="o", ls="", ms=4.2, color=PAL[g], label=g) for g in GROUP7]
fig.legend(handles=handles, loc="lower center", ncol=7, frameon=False, fontsize=7, handletextpad=0.1,
           columnspacing=1.0, bbox_to_anchor=(0.5, 0.0))
fig.tight_layout(rect=(0, 0.07, 1, 1), w_pad=0.8)
fig.savefig(OUT, dpi=500)
print("written", OUT, "| groups (B, C):", af2.tax_group7.value_counts().to_dict())
