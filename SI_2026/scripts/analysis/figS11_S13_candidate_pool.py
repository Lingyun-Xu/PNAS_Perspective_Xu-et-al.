#!/usr/bin/env python3
"""Figs. S11-S13 of the concise PNAS SI, drawn from Dataset S17 (folders FigS11, FigS12, FigS13).

Every scatter places proteins of the candidate pool against a template library, with the definitions of the selection
(select_plants.py):
  candidate pool    the 179,596 unique sequences of the 182,501 AFDB v6 plant accessions without a SIFTS-mapped PDB
                    chain; in a selection round, the candidates still in the pool
  template library  the 186,823 unique protein sequences of the September 2026 PDB ("complete PDB"), or the same
                    without the 3,825 plant sequences first released after 3 February 2021 (pdb3825; 182,998
                    sequences), plus added proteins where stated
  x  nearest template-pool distance (d_min): Euclidean distance between ESM-2 mean embeddings from the candidate to
     the closest sequence of the template pool (raw minimum, as in the selection)
  y  average homolog distance (d_avg): mean ESM-2 distance from the candidate to the other candidates of the candidate
     pool, itself not counted (in a selection round, to the other candidates still in the candidate pool; values of
     Dataset S17 as corrected by apply_davg_self_excluded.py, 6 Nov 2026)
Fig. S11 (plant_library_ABC.png): A complete PDB, B PDB without pdb3825, C density difference B minus A; black dots in
  A, the candidates whose closest library sequence is a pdb3825 sequence.
Fig. S12 (rounds_AF.png): sel3800 rounds 0 and 37 (A, B) and sel690 rounds 0 and 68 (D, E), the first and last round
  of each; library = complete PDB + the proteins selected in earlier rounds; black dots, the proteins selected in that
  round; C, F density difference last round minus round 0.
Fig. S13 (figS13_real_vs_ai_3800.png): the 175,537 candidates selected in no arm, against the baseline library (PDB
  without pdb3825; A), the baseline + pdb3825 (= complete PDB; B) and the baseline + 3,800 proteins selected on the
  baseline library (C); D, E density differences; F share of the candidates brought closer.
Density differences: -ln of the binned density (150 x 150 bins), each relative to its minimum, first minus second (as
in select_plants.py); bins of S11 and S12 from FigS12/limits.json (round 0 of the selection).
Replaces plant_library(), rounds() and s13() of si_pnas_cluster.py and figS11_plant_library.py.
Usage: python figS11_S13_candidate_pool.py <Dataset_S17 dir> <output dir>"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm

S17, OUT = Path(sys.argv[1]), Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)
GREEN = "#0eb519"
WIDE = {"font.size": 16, "axes.titlesize": 16, "axes.labelsize": 17, "xtick.labelsize": 15, "ytick.labelsize": 15}
LIM = json.loads((S17 / "FigS12/limits.json").read_text())
XLIM, YLIM = tuple(LIM["xlim"]), tuple(LIM["ylim"])


def free_energy_difference(x1, y1, x2, y2, ax, xlim, ylim, vmax=6):
    hist1 = np.histogram2d(x1, y1, bins=[150, 150], range=[list(xlim), list(ylim)], density=True)
    hist2 = np.histogram2d(x2, y2, bins=[150, 150], range=[list(xlim), list(ylim)], density=True)
    xedge, yedge = hist1[1], hist1[2]
    area = (xedge[1] - xedge[0]) * (yedge[1] - yedge[0])
    with np.errstate(divide="ignore"):
        free1 = -np.log(hist1[0] * area)
        free2 = -np.log(hist2[0] * area)
    free1[np.isinf(free1)] = np.nan
    free2[np.isinf(free2)] = np.nan
    delta1 = free1 - np.nanmin(free1)
    delta2 = free2 - np.nanmin(free2)
    nan1, nan2 = np.isnan(delta1), np.isnan(delta2)
    delta1 = np.where(nan1 & ~nan2, vmax, delta1)
    delta2 = np.where(nan2 & ~nan1, vmax, delta2)
    xx = (xedge[:-1] + xedge[1:]) / 2
    yy = (yedge[:-1] + yedge[1:]) / 2
    levels = np.linspace(-vmax, vmax, vmax * 10 + 1)
    return ax.contourf(xx, yy, (delta1 - delta2).T, levels, vmin=-vmax, vmax=vmax, cmap="seismic")


def fixed_limits(ax):
    ax.set_xlim(0, 6)
    ax.set_ylim(2, 8)
    ax.set_xticks(range(0, 7, 2))
    ax.set_yticks(range(2, 9, 2))


def style(ax, lw=1.2, length=4, labelsize=None):
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(lw)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(width=lw, length=length, **({"labelsize": labelsize} if labelsize else {}))


def letter(ax, text, x=-0.14, y=1.06, size=22):
    ax.text(x, y, text, transform=ax.transAxes, fontsize=size, fontweight="bold", va="bottom", ha="left")


def outside(x, y):
    return int(((x < 0) | (x > 6) | (y < 2) | (y > 8)).sum())


def draw_s11(fig, axes, legend_size=14, letter_size=22, letters="ABC"):
    """candidate pool against the September 2026 PDB template pool and the baseline template pool, and B minus A"""
    t = pd.read_csv(S17 / "FigS11/candidate_pool_round0.tsv.gz", sep="\t")
    x_all, x_ex, y = (t[c].to_numpy(float) for c in ("d_nearest_complete", "d_nearest_without_pdb3825",
                                                      "d_average_homolog"))
    black = t.nearest_is_pdb3825.to_numpy(bool)
    print(f"S11: candidates {len(t):,}; closest library sequence in pdb3825 {int(black.sum()):,} ({black.mean():.1%}); "
          f"median x complete {np.median(x_all):.3f} without {np.median(x_ex):.3f}; black {np.median(x_all[black]):.3f} "
          f"-> {np.median(x_ex[black]):.3f} (median increase {np.median((x_ex - x_all)[black]):.3f}); median y black "
          f"{np.median(y[black]):.3f} others {np.median(y[~black]):.3f}; outside A {outside(x_all, y)} B "
          f"{outside(x_ex, y)}")
    for ax, x, title in ((axes[0], x_all, "September 2026 PDB\ntemplate pool"),
                         (axes[1], x_ex, "Baseline template pool\n(without pdb3825)")):
        ax.scatter(x, y, s=3, c=GREEN, alpha=0.12, linewidths=0, rasterized=True)
        fixed_limits(ax)
        ax.set_xlabel("Nearest template-pool distance")
        ax.set_ylabel("Average homolog distance")
        ax.set_title(title)
        style(ax)
    axes[0].scatter(x_all[black], y[black], s=3, c="#111111", alpha=0.3, linewidths=0, rasterized=True,
                    label="Nearest template in pdb3825")
    leg = axes[0].legend(frameon=False, loc="upper right", fontsize=legend_size, markerscale=4, handletextpad=0.2)
    for h in leg.legend_handles:
        h.set_alpha(1)
    ax = axes[2]
    cf = free_energy_difference(x_ex, y, x_all, y, ax, XLIM, YLIM)
    fixed_limits(ax)
    ax.set_xlabel("Nearest template-pool distance")
    ax.set_ylabel("Average homolog distance")
    ax.set_title(f"{letters[1]} minus {letters[0]}")
    style(ax)
    bar = fig.colorbar(cf, ax=ax, fraction=0.046, pad=0.04, ticks=range(-6, 7, 2))
    bar.set_label("Protein density difference")
    for a, t_ in zip(axes, letters):
        letter(a, t_, size=letter_size)


def draw_s12(fig, axes, legend_size=14, letter_size=22, letters="ABCDEF"):
    """selection rounds: sel3800 rounds 0 and 37 and their difference (row 0), sel690 rounds 0 and 68 (row 1)"""
    rows = (("sel3800", "sel3800_rounds_0_and_37.tsv.gz", 37, 14), ("sel690", "sel690_rounds_0_and_68.tsv.gz", 68, 40))
    for r, (name, path, last, size) in enumerate(rows):
        t = pd.read_csv(S17 / "FigS12" / path, sep="\t")
        data = {}
        for rnd in (0, last):
            s = t[t["round"] == rnd].sort_values("pool_position")
            data[rnd] = (s.d_nearest_library.to_numpy(float), s.d_average_homolog.to_numpy(float),
                         s.selected_this_round.to_numpy(bool))
            mn, av, sel = data[rnd]
            print(f"S12 {name} round {rnd}: candidates {len(mn):,}; selected {int(sel.sum())}; median x pool "
                  f"{np.median(mn):.3f} selected {np.median(mn[sel]):.3f}; median y pool {np.median(av):.3f} "
                  f"selected {np.median(av[sel]):.3f}; outside {outside(mn, av)}")
        for c, rnd in enumerate((0, last)):
            ax = axes[r][c]
            mn, av, sel = data[rnd]
            ax.scatter(mn, av, s=3, c=GREEN, alpha=0.12, linewidths=0, rasterized=True)
            ax.scatter(mn[sel], av[sel], s=size, c="#111111", alpha=0.85, linewidths=0, label="Selected")
            ax.legend(frameon=False, loc="upper right", fontsize=legend_size, markerscale=2 if size < 20 else 1.2)
            fixed_limits(ax)
            ax.set_xlabel("Nearest template-pool distance")    # September 2026 PDB + proteins selected earlier
            ax.set_ylabel("Average homolog distance")
            ax.set_title(f"{name}, round {rnd}\n{len(mn):,} remaining")
            style(ax)
        ax = axes[r][2]
        cf = free_energy_difference(data[last][0], data[last][1], data[0][0], data[0][1], ax, XLIM, YLIM)
        fixed_limits(ax)
        ax.set_xlabel("Nearest template-pool distance")
        ax.set_ylabel("Average homolog distance")
        ax.set_title(f"{name}, round {last}\nminus round 0")
        style(ax)
        bar = fig.colorbar(cf, ax=ax, fraction=0.046, pad=0.04, ticks=range(-6, 7, 2))
        bar.set_label("Protein density difference")
    for a, t_ in zip([a for row in axes for a in row], letters):
        letter(a, t_, size=letter_size)


def fig_s11():
    with plt.rc_context(WIDE):
        fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.9), gridspec_kw={"width_ratios": [1, 1, 1.18]})
        draw_s11(fig, axes)
        fig.tight_layout()
        fig.savefig(OUT / "plant_library_ABC.png", dpi=300)
        plt.close(fig)


def fig_s12():
    with plt.rc_context(WIDE):
        fig, axes = plt.subplots(2, 3, figsize=(15.5, 10.0), gridspec_kw={"width_ratios": [1, 1, 1.18]})
        draw_s12(fig, axes)
        fig.tight_layout(h_pad=2.0)
        fig.savefig(OUT / "rounds_AF.png", dpi=300)
        plt.close(fig)


def fig_s10_combined():
    """Fig. S10 of the SI for the rebuttal: the candidate pool figure (A-C) above the selection rounds (D-I), one
    figure; fonts a little larger than in the separate figures, because the combined figure is printed smaller"""
    big = {"font.size": 18, "axes.titlesize": 18, "axes.labelsize": 19, "xtick.labelsize": 17, "ytick.labelsize": 17}
    with plt.rc_context(big):
        fig, axes = plt.subplots(3, 3, figsize=(15.5, 15.6), gridspec_kw={"width_ratios": [1, 1, 1.18]})
        draw_s11(fig, axes[0], legend_size=16, letter_size=26, letters="ABC")
        draw_s12(fig, axes[1:], legend_size=16, letter_size=26, letters="DEFGHI")
        fig.tight_layout(h_pad=2.2)
        fig.savefig(OUT / "candidate_pool_and_rounds_AI.png", dpi=300)
        plt.close(fig)


def fig_s13():
    t = pd.read_parquet(S17 / "FigS13/evaluation_distances.parquet")
    base, real, ai, hom = (t[k].to_numpy(float) for k in ("base_3800", "real_3800", "ai_3800", "d_average_homolog"))
    xlim = (0.0, float(np.percentile(base, 99.9)) * 1.05)
    ylim = (float(np.percentile(hom, 0.1)) * 0.97, float(np.percentile(hom, 99.9)) * 1.03)
    print(f"S13: candidates {len(t):,}; median y {np.median(hom):.3f}; outside A {outside(base, hom)} B "
          f"{outside(real, hom)} C {outside(ai, hom)}; any decrease real {(base - real > 1e-4).mean():.1%} AI "
          f"{(base - ai > 1e-4).mean():.1%}")
    with plt.rc_context({"font.size": 6.5}):
        fig = plt.figure(figsize=(6.5, 4.9), constrained_layout=True)
        gs = fig.add_gridspec(2, 3)
        axes = [[fig.add_subplot(gs[r, k]) for k in range(3)] for r in range(2)]
        for ax, x, name in ((axes[0][0], base, "Baseline library"),
                            (axes[0][1], real, "+ 3,825 real plant PDB sequences"),
                            (axes[0][2], ai, "+ 3,800 AI-selected plant proteins")):
            h = ax.hist2d(x, hom, bins=150, range=[xlim, ylim], cmap="plasma_r", norm=LogNorm(vmin=1, vmax=100), cmin=1)
            ax.set_title(name, fontsize=6.5)
        cb = fig.colorbar(h[3], ax=axes[0], orientation="horizontal", location="bottom", fraction=0.05, aspect=60,
                          shrink=0.5)
        cb.set_label("Proteins per bin", fontsize=6.5)
        cb.ax.tick_params(labelsize=6, width=0.6, length=2)
        for ax, x, name in ((axes[1][0], real, "Real minus baseline"), (axes[1][1], ai, "AI-selected minus baseline")):
            cf = free_energy_difference(x, hom, base, hom, ax, xlim, ylim)
            ax.set_title(name, fontsize=6.5)
        cb = fig.colorbar(cf, ax=axes[1][:2], orientation="horizontal", location="bottom", fraction=0.05, aspect=40,
                          shrink=0.5, ticks=range(-6, 7, 2))
        cb.set_label("Protein density difference", fontsize=6.5)
        cb.ax.tick_params(labelsize=6, width=0.6, length=2)
        for ax in axes[0] + axes[1][:2]:
            fixed_limits(ax)
            ax.set_xlabel("Nearest template distance", fontsize=6.5)
            ax.set_ylabel("Average homolog distance", fontsize=6.5)
            style(ax, lw=0.7, length=2, labelsize=6)
        ax = axes[1][2]
        grid = np.linspace(0, 2.0, 201)
        for x, name, color in ((real, "Real (pdb3825)", "#555555"), (ai, "AI-selected", "#d95f02")):
            gain = base - x
            share = [(gain > 1e-4).mean() * 100] + [(gain >= g).mean() * 100 for g in grid[1:]]
            ax.plot(grid, share, color=color, lw=1.4, label=name)
        ax.set_xlabel("Decrease in nearest template\ndistance (at least)", fontsize=6.5)
        ax.set_ylabel("Evaluated plant proteins (%)", fontsize=6.5)
        ax.set_title("Plant proteins brought closer", fontsize=6.5)
        ax.legend(frameon=False, fontsize=6)
        ax.set_xlim(0, 2.0)
        style(ax, lw=0.7, length=2, labelsize=6)
        for a, l in zip(axes[0] + axes[1], "ABCDEF"):
            letter(a, l, x=-0.3, y=1.05, size=9)
        fig.savefig(OUT / "figS13_real_vs_ai_3800.png", dpi=500)
        plt.close(fig)


if __name__ == "__main__":
    # optional third and later arguments choose the figures (s11 s12 s13); default all three. The former Fig. S13
    # (figS13_real_vs_ai_3800.png) is now in the main text; the SI Fig. S13 is figS13_average_homolog_distance.py
    # "s10ai": the rebuttal SI's combined Fig. S10 (candidate_pool_and_rounds_AI.png, panels A-I)
    todo = [a.lower() for a in sys.argv[3:]] or ["s11", "s12", "s13"]
    for name, fn in (("s11", fig_s11), ("s12", fig_s12), ("s13", fig_s13), ("s10ai", fig_s10_combined)):
        if name in todo:
            fn()
    print("wrote", OUT)
