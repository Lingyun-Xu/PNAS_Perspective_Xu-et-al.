#!/usr/bin/env python3
"""Redraw HDBSCAN-on-UMAP maps for the SI from the saved cluster tables (no re-clustering).

Same palette, point sizes and noise colour as hdbscan_*_clean.py; the per-cluster legend is
replaced by a "K clusters + noise" label so maps with ~100 clusters stay legible.
Also writes a 2x2 panel (pdb3825 K97, sel3800 K109, sel690 K5, sel690 K10).
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DR = Path(r"/path/to/project/template_esm/campaign/dr")
OUT = Path(r"/path/to/project/SI_2026/figures")
NOISE = "#B0B0B0"
CB = [
    "#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00", "#56B4E9", "#F0E442",
    "#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677",
    "#882255", "#AA4499", "#661100", "#6699CC", "#AA3377", "#5AAE61", "#FDB863",
    "#B2ABD2", "#E66101", "#5E3C99", "#1B7837", "#762A83", "#C51B7D", "#01665E",
    "#8C510A", "#35978F", "#BF812D", "#543005", "#80CDC1", "#DFC27D", "#018571",
    "#A6DBA0", "#C2A5CF", "#7B3294", "#008837", "#2C7BB6",
]
MAPS = [
    ("pdb3825", "pdb3825/hdbscan_clean/hdbscan_on_umap2_nn30_md0.45_K97_best_mcs14_ms3_eps0.05_eom_clusters.tsv",
     "PDB plant set (pdb3825), n = 3,825"),
    ("sel3800", "sel3800/hdbscan_clean_K0_200/hdbscan_on_umap2_nn30_md0.45_K109_best_mcs12_ms3_eps0.05_leaf_clusters.tsv",
     "Selection sel3800, n = 3,800"),
    ("sel690_K5", "sel690/hdbscan_clean/hdbscan_on_umap5_nn30_md0.45_K5_best_mcs13_ms13_eps0.0_eom_clusters.tsv",
     "Selection sel690, n = 690 (optimum)"),
    ("sel690_K10", "sel690/hdbscan_clean/hdbscan_on_umap2_nn30_md0.45_K10_fine_mcs17_ms3_eps0.0_leaf_clusters.tsv",
     "Selection sel690, n = 690 (finer view)"),
]


def draw(ax, df, title):
    lab = df.cluster.to_numpy()
    x, y = df.umap1.to_numpy(), df.umap2.to_numpy()
    noise = lab < 0
    ax.scatter(x[noise], y[noise], s=4, c=NOISE, alpha=0.35, linewidths=0, rasterized=True)
    clusters = sorted({int(c) for c in lab if c >= 0})
    for i, c in enumerate(clusters):
        m = lab == c
        ax.scatter(x[m], y[m], s=6, color=CB[i % len(CB)], alpha=0.75, linewidths=0, rasterized=True)
    ax.set_title(title, fontsize=10)
    ax.text(0.02, 0.98, f"{len(clusters)} clusters + noise ({noise.mean():.0%})",
            transform=ax.transAxes, va="top", fontsize=8)
    ax.set_xlabel("UMAP 1", fontsize=8)
    ax.set_ylabel("UMAP 2", fontsize=8)
    ax.tick_params(labelsize=7)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    tables = {key: pd.read_csv(DR / rel, sep="\t") for key, rel, _ in MAPS}
    for key, _, title in MAPS:
        fig, ax = plt.subplots(figsize=(6.4, 5.0))
        draw(ax, tables[key], title)
        fig.tight_layout()
        fig.savefig(OUT / f"umap_hdbscan_{key}.png", dpi=220)
        plt.close(fig)
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 8.6))
    for ax, (key, _, title), tag in zip(axes.ravel(), MAPS, "ABCD"):
        draw(ax, tables[key], title)
        ax.text(-0.12, 1.04, tag, transform=ax.transAxes, fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUT / "umap_hdbscan_panel.png", dpi=200)
    plt.close(fig)
    for key in tables:
        n = len(tables[key]); k = tables[key].cluster.max() + 1
        print(key, "n", n, "K", k, "noise", int((tables[key].cluster < 0).sum()))


if __name__ == "__main__":
    main()
