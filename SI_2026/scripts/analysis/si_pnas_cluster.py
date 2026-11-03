#!/usr/bin/env python3
"""Combined figures for the PNAS-formatted concise SI, redrawn on the cluster from the saved data (no data recomputed).

  dtemplate      d_template,min versus d_nearest (as of model building / September 2026 PDB): 2 x 3 panels A-F
                 (full SI Figs. S18-S19 in one figure; legend fonts 1.5 x the original 7 pt)
  s6, s15_18     PNAS Figs. S6 and S15-S18
Figs. S11-S13 (plant_library_ABC.png, rounds_AF.png, figS13_real_vs_ai_3800.png) are drawn by
figS11_S13_candidate_pool.py from Dataset S17.

Output: /path/to/cluster/si_pnas/out/*.png (300 dpi). Original files are not touched.
Usage: python si_pnas_cluster.py [dtemplate|s6|s15_18 ...]
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path("/path/to/cluster/redo_PDB/ESM_embedding/template_v6")
CAMP = ROOT / "campaign"
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(CAMP))
sys.path.insert(0, "/path/to/cluster/si_pnas")
OUT = Path("/path/to/cluster/si_pnas/out")
OUT.mkdir(parents=True, exist_ok=True)
GREEN = "#0eb519"
plt.rcParams.update({"font.size": 13, "axes.titlesize": 13, "axes.labelsize": 13, "xtick.labelsize": 12,
                     "ytick.labelsize": 12})
WIDE = {"font.size": 16, "axes.titlesize": 16, "axes.labelsize": 17, "xtick.labelsize": 15, "ytick.labelsize": 15}


def ptext(p):
    return "p < 10$^{-300}$" if p == 0 else f"p = {p:.2g}"


def letter(ax, text, x=-0.14, y=1.06):
    ax.text(x, y, text, transform=ax.transAxes, fontsize=22, fontweight="bold", va="bottom", ha="left")


def fixed_limits(ax, x=(0, 6), y=(2, 8)):
    """Figs. S11-S12: x 0-6, y 2-8, ticks every 2"""
    ax.set_xlim(*x)
    ax.set_ylim(*y)
    ax.set_xticks(range(x[0], x[1] + 1, 2))
    ax.set_yticks(range(y[0], y[1] + 1, 2))


def style(ax):
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(1.2)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(width=1.2, length=4)


# ---------------------------------------------------------------- PNAS Figs. S6 and S7 (template distance)
GROUPS7 = ["Bacteria", "Archaea", "Animalia", "Fungi", "Protist", "Parasitic Worms", "Viridiplantae"]


def small_letter(ax, text, x=-0.2, y=1.03, size=10):
    ax.text(x, y, text, transform=ax.transAxes, fontsize=size, fontweight="bold", va="bottom", ha="left")


def thin(ax, lw=0.9, length=2.5, labelsize=7):
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(lw)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(width=lw, length=length, labelsize=labelsize)


def group_scatter(ax, use, xcol, ycol, s):
    import plot_template_distances as ptd
    rng = np.random.default_rng(0)
    counts = {}
    for g in ptd.DRAW:
        part = use[use.tax_group7 == g]
        counts[g] = len(part)
        if part.empty:
            continue
        part = part.iloc[rng.permutation(len(part))]
        ax.scatter(part[xcol], part[ycol], s=s, c=ptd.PALETTE[g], alpha=ptd.ALPHA[g], linewidths=0, rasterized=True)
    return counts


def bottom_legend(fig, fontsize):
    import plot_template_distances as ptd
    handles = [plt.Line2D([0], [0], marker="o", ls="", ms=5, color=ptd.PALETTE[g], label=g) for g in GROUPS7]
    fig.legend(handles=handles, loc="lower center", ncol=7, frameon=False, fontsize=fontsize, handletextpad=0.1,
               columnspacing=1.0, bbox_to_anchor=(0.5, 0.0))


def panels_of(df):
    return [("Without own PDB structure", df[~df.in_pdb.astype(bool)]),
            ("With own PDB structure", df[df.in_pdb.astype(bool)]),
            ("All models with templates", df)]


def s6():
    from exact_p import p_text
    df = pd.read_parquet(ROOT / "distances/template_distances.parquet")
    df = df[df.n_templates.fillna(0) > 0].copy()
    fig, axes = plt.subplots(1, 3, figsize=(6.3, 2.75))
    out = []
    for ax, (name, sub), t in zip(axes, panels_of(df), "ABC"):
        use = sub.dropna(subset=["d_template_min", "plddt_global"])
        counts = group_scatter(ax, use, "d_template_min", "plddt_global", 0.6)
        rho, _ = stats.spearmanr(use.d_template_min, use.plddt_global)
        ax.set_title(f"{name}\nρ = {rho:.3f}, {p_text(rho, len(use))}".replace("ρ = -", "ρ = −"), fontsize=7.5)
        ax.set_xlim(0, 10)
        ax.set_xticks([0, 5, 10])
        ax.set_xlabel("Closest-template distance", fontsize=8)
        ax.set_ylabel("Global pLDDT", fontsize=8)
        thin(ax, labelsize=7.5)
        small_letter(ax, t, x=-0.27, y=1.1, size=9)
        out.append({"panel": name, "n": len(use), "spearman": rho, "p": p_text(rho, len(use), False), **counts})
    bottom_legend(fig, 7.5)
    fig.tight_layout(rect=(0, 0.07, 1, 1), w_pad=1.0)
    fig.savefig(OUT / "figS6_plddt_vs_dtemplate.png", dpi=500)
    plt.close(fig)
    print(pd.DataFrame(out).to_string(index=False), flush=True)


def dtemplate():
    from exact_p import p_text
    df = pd.read_parquet(ROOT / "distances/template_distances.parquet")
    df = df[df.n_templates.fillna(0) > 0].copy()
    rows = [("d_nearest_asof", "Nearest PDB distance (at model building)"),
            ("d_nearest_now", "Nearest PDB distance (September 2026)")]
    fig, axes = plt.subplots(2, 3, figsize=(6.5, 5.3))
    out = []
    for r, (xcol, xlabel) in enumerate(rows):
        for c, (name, sub) in enumerate(panels_of(df)):
            ax = axes[r, c]
            use = sub.dropna(subset=[xcol, "d_template_min"])
            counts = group_scatter(ax, use, xcol, "d_template_min", 1.0)
            x, y = use[xcol].to_numpy(float), use.d_template_min.to_numpy(float)
            rho, _ = stats.spearmanr(x, y)
            lo, hi = float(min(x.min(), y.min())), float(max(x.max(), y.max()))
            ax.plot([lo, hi], [lo, hi], color="black", lw=0.8)
            ax.set_title(f"{name}\nρ = {rho:.3f}, {p_text(rho, len(use))}", fontsize=6.5)
            ax.set_xlim(0, 8)
            ax.set_ylim(0, 10)
            ax.set_xticks(range(0, 9, 2))
            ax.set_yticks(range(0, 11, 2))
            ax.set_xlabel(xlabel, fontsize=6.5)
            ax.set_ylabel("Closest-template distance", fontsize=6.5)
            thin(ax, lw=0.8, length=2.2, labelsize=6)
            small_letter(ax, "ABCDEF"[3 * r + c], x=-0.25, y=1.1, size=10)
            out.append({"x": xcol, "panel": name, "n": len(use), "spearman": rho,
                        "p": p_text(rho, len(use), False), **counts})
    bottom_legend(fig, 7)
    fig.tight_layout(rect=(0, 0.04, 1, 1), h_pad=1.6, w_pad=0.8)
    fig.savefig(OUT / "dtemplate_vs_dnearest_AF.png", dpi=500)
    plt.close(fig)
    print(pd.DataFrame(out).to_string(index=False), flush=True)


# ---------------------------------------------------------------- shared by PNAS Figs. S15-S18
COUNT_NORM = dict(vmin=1, vmax=100)          # proteins per bin, 10^0-10^2 (log scale)
DENSITY_LABEL = "Protein density difference"  # free-energy difference of the binned densities, -6 to 6


def cbar_ticks(cb, size):
    cb.ax.tick_params(labelsize=size, width=0.6, length=2)


# ---------------------------------------------------------------- PNAS Figs. S15-S18 (search density change)
def s15_18():
    from matplotlib.colors import LogNorm
    import select_plants as sp
    d = CAMP / "figures" / "search_density_change_hhprob50"
    coh = {s: pd.read_csv(d / f"{s}_cohort.tsv", sep="\t") for s in ("hhsearch", "foldseek")}
    allx = np.concatenate([df[[c for c in df.columns if c.startswith("d_") and c != "d_homolog"]].to_numpy(float).ravel()
                           for df in coh.values()])
    ally = np.concatenate([df.d_homolog.to_numpy(float) for df in coh.values()])
    xlim = (0.0, float(np.nanpercentile(allx, 99.9)) * 1.05)
    ylim = (float(np.nanpercentile(ally, 0.1)) * 0.97, float(np.nanpercentile(ally, 99.9)) * 1.03)
    labels = {"sel3800": ("September 2026 PDB", "PDB + sel3800"),
              "pdb3825": ("PDB without the 3,825\nrecent plant sequences", "September 2026 PDB")}
    for search, name in (("hhsearch", "HHsearch"), ("foldseek", "Foldseek")):
        for comp in ("sel3800", "pdb3825"):
            df = pd.read_csv(d / f"{search}_{comp}_pairs.tsv", sep="\t")
            x0, x1, y = df.d_before.to_numpy(float), df.d_after.to_numpy(float), df.d_homolog.to_numpy(float)
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
            cf = sp.free_energy_difference(x1, y, x0, y, ax, xlim, ylim)
            ax.set_xlim(*xlim)
            ax.set_ylim(*ylim)
            cb = fig.colorbar(cf, ax=ax, fraction=0.046, pad=0.03, shrink=1.0, aspect=40, ticks=range(-6, 7, 2))
            cb.set_label(DENSITY_LABEL, fontsize=6.5)
            cbar_ticks(cb, 6)
            ax.set_title("After minus before", fontsize=7)
            for ax in (axes[0][0], axes[0][1], axes[1][0]):
                ax.set_xlabel("ESM-2 distance to top search hit", fontsize=7)
                ax.set_ylabel("Average homolog distance", fontsize=7)
                thin(ax, lw=0.7, length=2, labelsize=6.5)
            ax = axes[1][1]
            g = gain[moved & added]
            if len(g):
                lo, hi = np.percentile(g, [0.5, 99.5])
                ax.hist(g, bins=np.linspace(min(lo, -0.1), hi, 60), color="#d95f02", edgecolor="white", linewidth=0.3)
                ax.axvline(0, color="black", lw=0.8)
                ax.axvline(np.median(g), color="#555555", lw=1.0, ls="--")
                ax.set_ylim(0, ax.get_ylim()[1] * 1.15)
            ax.set_xlabel("Distance before − after (positive = closer)", fontsize=7)
            ax.set_ylabel("Proteins whose top hit is added", fontsize=7)
            ax.set_title("Change when the top hit is an added structure", fontsize=7)
            thin(ax, lw=0.7, length=2, labelsize=6.5)
            for a, l in zip(axes[0] + axes[1], "ABCD"):
                small_letter(a, l, x=-0.24, y=1.04, size=10)
            fig.suptitle(f"{name}, {len(df):,} plant proteins without a structure", fontsize=7.5)
            fig.savefig(OUT / f"figS_{search}_{comp}_density_change.png", dpi=500)
            plt.close(fig)

if __name__ == "__main__":
    for job in sys.argv[1:] or ["dtemplate"]:
        {"dtemplate": dtemplate, "s6": s6, "s15_18": s15_18}[job]()
        print("DONE", job, flush=True)
