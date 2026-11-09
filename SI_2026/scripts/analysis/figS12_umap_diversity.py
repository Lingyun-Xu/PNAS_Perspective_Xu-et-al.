#!/usr/bin/env python3
"""Fig. S12 of the SI for the rebuttal (7 Nov 2026), from Dataset S17 folder FigS14:
  A, B  HDBSCAN clusters on the UMAP of each set (pdb3825 K = 97; sel3800 K = 109), drawn as in
        redraw_hdbscan_maps.py (the former panels A and B)
  C, D  one UMAP of both sets together (embeddings standardized over the 7,625 proteins of both sets, PCA50, UMAP
        n_neighbors 30, min_dist 0.45; coordinates of plot_dr_diversity_umap_opt.py), C pdb3825 and D sel3800 with
        the other set in grey; colour, the distance from each protein to the nearest other protein of either set in the
        standardized 1,280-dimensional space (NN1; PowerNorm gamma 0.55 between the 2nd and 98th percentiles, magma,
        as in joint_umap_colored_by_nn1.png). No titles: set names, sizes and the colour explanation are in the legend.
The colour bar has its own column at the right edge, so it does not cover panel D.
Usage: python figS12_umap_diversity.py <Dataset_S17 dir> <output .png>"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import PowerNorm

sys.path.insert(0, str(Path(__file__).resolve().parent))
import redraw_hdbscan_maps as rh  # noqa: E402

S17, OUT = Path(sys.argv[1]), Path(sys.argv[2])
F = S17 / "FigS14"
AB = [("hdbscan_on_umap2_nn30_md0.45_K97_best_mcs14_ms3_eps0.05_eom_clusters.tsv", "PDB plant set (pdb3825), n = 3,825"),
      ("hdbscan_on_umap2_nn30_md0.45_K109_best_mcs12_ms3_eps0.05_leaf_clusters.tsv", "Selection sel3800, n = 3,800")]

joint = pd.read_csv(F / "joint_umap_nn1_pdb3825_sel3800.tsv", sep="\t")
xy = joint[["umap1", "umap2"]].to_numpy(float)
nn1 = joint.nn1_joint_z.to_numpy(float)
cohort = joint.cohort.to_numpy()
norm = PowerNorm(gamma=0.55, vmin=float(np.quantile(nn1, 0.02)), vmax=float(np.quantile(nn1, 0.98)))
lim_x = (xy[:, 0].min() - 0.4, xy[:, 0].max() + 0.4)
lim_y = (xy[:, 1].min() - 0.4, xy[:, 1].max() + 0.4)

fig = plt.figure(figsize=(10.5, 8.6))
gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 0.045], wspace=0.28, hspace=0.32)
axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])]
cax = fig.add_subplot(gs[1, 2])
for ax, (name, title) in zip(axes[:2], AB):
    rh.draw(ax, pd.read_csv(F / name, sep="\t"), title)
for ax, own in ((axes[2], "pdb3825"), (axes[3], "sel3800")):
    m = cohort == own
    ax.scatter(xy[~m, 0], xy[~m, 1], s=3, c="#DDDDDD", alpha=0.25, linewidths=0, rasterized=True, zorder=0)
    sc = ax.scatter(xy[m, 0], xy[m, 1], c=nn1[m], s=6, cmap="magma", norm=norm, alpha=0.85, linewidths=0,
                    rasterized=True, zorder=2)
    ax.set_xlim(*lim_x)
    ax.set_ylim(*lim_y)
    ax.set_xlabel("Joint UMAP 1", fontsize=8)
    ax.set_ylabel("Joint UMAP 2", fontsize=8)
    ax.tick_params(labelsize=7)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    print(f"{own}: n {int(m.sum()):,}; median NN1 {np.median(nn1[m]):.2f}")
bar = fig.colorbar(sc, cax=cax, ticks=[t for t in (5, 10, 15, 20, 25, 30) if norm.vmin <= t <= norm.vmax])
bar.ax.tick_params(labelsize=7)
bar.set_label("NN1 distance", fontsize=8)
for ax, tag in zip(axes, "ABCD"):
    ax.text(-0.12, 1.04, tag, transform=ax.transAxes, fontsize=15, fontweight="bold")
fig.savefig(OUT, dpi=200, bbox_inches="tight", pad_inches=0.05)
print(f"colour scale {norm.vmin:.2f}-{norm.vmax:.2f} (2nd-98th percentile of {len(nn1):,} proteins); wrote {OUT}")
