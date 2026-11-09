#!/usr/bin/env python3
"""Figs. S15-S18 of the concise PNAS SI (HHsearch and Foldseek top-hit changes), drawn from Dataset S17 (FigS15_S18).

Local version of s15_18() in si_pnas_cluster.py (same data, bins, colours and layout); only the axis names differ:
  x  top-hit distance: the Euclidean distance between the ESM-2 (650M) mean embeddings of a plant protein and of its
     top search hit (HHsearch, the best hit with probability >= 50%; Foldseek, the best hit with E <= 1e-3)
  y  average homolog distance (d_avg): mean distance to the other candidates of the candidate pool (Dataset S17 as
     corrected on 6 Nov 2026; before, the mean distance to all plant sequences of AFDB v6, itself included)
Files: <S17>/FigS15_S18/{hhsearch,foldseek}_cohort.tsv (axis limits) and {search}_{sel3800,pdb3825}_pairs.tsv
(d_before, d_after, d_average_homolog, top_is_added). Writes figS_{search}_{comp}_density_change.png.
Usage: python figS15_S18_search_density.py <Dataset_S17 dir> <output dir>"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm

S17, OUT = Path(sys.argv[1]) / "FigS15_S18", Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)
COUNT_NORM = dict(vmin=1, vmax=100)          # proteins per bin, 10^0-10^2 (log scale)
DENSITY_LABEL = "Protein density difference"  # free-energy difference of the binned densities, -6 to 6
XLABEL = "Top-hit distance"
YLABEL = "Average homolog distance"


def free_energy_difference(x1, y1, x2, y2, ax, xlim, ylim, vmax=6):
    """as select_plants.free_energy_difference: -ln of the binned densities (150 x 150), each relative to its minimum,
    first minus second"""
    hist1 = np.histogram2d(x1, y1, bins=[150, 150], range=[list(xlim), list(ylim)], density=True)
    hist2 = np.histogram2d(x2, y2, bins=[150, 150], range=[list(xlim), list(ylim)], density=True)
    xedge, yedge = hist1[1], hist1[2]
    area = (xedge[1] - xedge[0]) * (yedge[1] - yedge[0])
    with np.errstate(divide="ignore"):
        free1 = -np.log(hist1[0] * area)
        free2 = -np.log(hist2[0] * area)
    free1[np.isinf(free1)] = np.nan
    free2[np.isinf(free2)] = np.nan
    if not np.isfinite(free1).any():
        free1 = np.zeros_like(free1)
    if not np.isfinite(free2).any():
        free2 = np.zeros_like(free2)
    delta1 = free1 - np.nanmin(free1)
    delta2 = free2 - np.nanmin(free2)
    nan1, nan2 = np.isnan(delta1), np.isnan(delta2)
    delta1 = np.where(nan1 & ~nan2, vmax, delta1)
    delta2 = np.where(nan2 & ~nan1, vmax, delta2)
    xx = (xedge[:-1] + xedge[1:]) / 2
    yy = (yedge[:-1] + yedge[1:]) / 2
    levels = np.linspace(-vmax, vmax, vmax * 10 + 1)
    return ax.contourf(xx, yy, (delta1 - delta2).T, levels, vmin=-vmax, vmax=vmax, cmap="seismic")


def small_letter(ax, text, x=-0.2, y=1.03, size=10):
    ax.text(x, y, text, transform=ax.transAxes, fontsize=size, fontweight="bold", va="bottom", ha="left")


def thin(ax, lw=0.9, length=2.5, labelsize=7):
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(lw)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(width=lw, length=length, labelsize=labelsize)


def cbar_ticks(cb, size):
    cb.ax.tick_params(labelsize=size, width=0.6, length=2)


def main():
    plt.rcParams.update({"font.size": 13, "axes.titlesize": 13, "axes.labelsize": 13, "xtick.labelsize": 12,
                         "ytick.labelsize": 12})
    coh = {s: pd.read_csv(S17 / f"{s}_cohort.tsv", sep="\t") for s in ("hhsearch", "foldseek")}
    allx = np.concatenate([df[[c for c in df.columns if c.startswith("d_") and c != "d_average_homolog"]].to_numpy(float).ravel()
                           for df in coh.values()])
    ally = np.concatenate([df.d_average_homolog.to_numpy(float) for df in coh.values()])
    xlim = (0.0, float(np.nanpercentile(allx, 99.9)) * 1.05)
    ylim = (float(np.nanpercentile(ally, 0.1)) * 0.97, float(np.nanpercentile(ally, 99.9)) * 1.03)
    labels = {"sel3800": ("September 2026 PDB template pool", "September 2026 PDB template pool + sel3800"),
              "pdb3825": ("baseline template pool", "September 2026 PDB template pool")}
    for search, name in (("hhsearch", "HHsearch"), ("foldseek", "Foldseek")):
        for comp in ("sel3800", "pdb3825"):
            df = pd.read_csv(S17 / f"{search}_{comp}_pairs.tsv", sep="\t")
            x0, x1, y = df.d_before.to_numpy(float), df.d_after.to_numpy(float), df.d_average_homolog.to_numpy(float)
            added = df.top_is_added.astype(str).str.lower().eq("true").to_numpy()
            gain = x0 - x1
            moved = np.abs(gain) > 1e-6
            fig = plt.figure(figsize=(6.0, 5.2), constrained_layout=True)
            gs = fig.add_gridspec(2, 2)
            axes = [[fig.add_subplot(gs[r, k]) for k in range(2)] for r in range(2)]
            for ax, x, ttl in ((axes[0][0], x0, f"Before: {labels[comp][0]}"), (axes[0][1], x1, f"After: {labels[comp][1]}")):
                h = ax.hist2d(x, y, bins=150, range=[xlim, ylim], cmap="plasma_r", norm=LogNorm(**COUNT_NORM), cmin=1)
                ax.set_title(ttl, fontsize=7)
            cb = fig.colorbar(h[3], ax=axes[0], fraction=0.046, pad=0.02, shrink=1.0, aspect=40)   # one colour bar for A and B
            cb.set_label("Proteins per bin", fontsize=6.5)
            cbar_ticks(cb, 6)
            ax = axes[1][0]
            cf = free_energy_difference(x1, y, x0, y, ax, xlim, ylim)
            ax.set_xlim(*xlim)
            ax.set_ylim(*ylim)
            cb = fig.colorbar(cf, ax=ax, fraction=0.046, pad=0.03, shrink=1.0, aspect=40, ticks=range(-6, 7, 2))
            cb.set_label(DENSITY_LABEL, fontsize=6.5)
            cbar_ticks(cb, 6)
            ax.set_title("After minus before", fontsize=7)
            for ax in (axes[0][0], axes[0][1], axes[1][0]):
                ax.set_xlabel(XLABEL, fontsize=7)
                ax.set_ylabel(YLABEL, fontsize=7)
                thin(ax, lw=0.7, length=2, labelsize=6.5)
            ax = axes[1][1]
            g = gain[moved & added]
            if len(g):
                lo, hi = np.percentile(g, [0.5, 99.5])
                ax.hist(g, bins=np.linspace(min(lo, -0.1), hi, 60), color="#d95f02", edgecolor="white", linewidth=0.3)
                ax.axvline(0, color="black", lw=0.8)
                ax.axvline(np.median(g), color="#555555", lw=1.0, ls="--")
                ax.set_ylim(0, ax.get_ylim()[1] * 1.15)
            ax.set_xlabel("Top-hit distance before − after (positive = closer)", fontsize=7)
            ax.set_ylabel("Proteins whose top hit is added", fontsize=7)
            ax.set_title("Change when the top hit is an added structure", fontsize=7)
            thin(ax, lw=0.7, length=2, labelsize=6.5)
            for a, l in zip(axes[0] + axes[1], "ABCD"):
                small_letter(a, l, x=-0.24, y=1.04, size=10)
            fig.suptitle(f"{name}, {len(df):,} candidates", fontsize=7.5)
            fig.savefig(OUT / f"figS_{search}_{comp}_density_change.png", dpi=500)
            plt.close(fig)
            print(f"{search} {comp}: {len(df):,} proteins; top hit added {int(added.sum()):,}; switched to an added "
                  f"structure and moved {int((moved & added).sum()):,} (closer {int((gain[moved & added] > 0).sum()):,}, "
                  f"median gain {np.median(g):.2f})", flush=True)


if __name__ == "__main__":
    main()
