#!/usr/bin/env python3
"""Analysis server: figures of the PNAS-formatted concise SI redrawn from the original data with the original plotting logic.
Figures are drawn at the size they are printed (width in inches as placed in the SI), so point sizes are final.

  s1   Fig. S1  pLDDT density, 47 proteomes + Swiss-Prot, 7 groups (afdb_v6/plot_merged_7groups.py data)
  s3   Fig. S3  mapped pLDDT versus lDDT-Ca (sifts_accuracy/09_analyze.py data)
  s4   Fig. S4  AFDB v6 global pLDDT versus ESMFold pLDDT (esm_atlas/plot_precise_esm.py data)
  s5   Fig. S5  global pLDDT versus AFDB MSA depth (A) and log10(Neff + 1) (B) (plot_scatter_paper_style / plot_neff)
  s9   Fig. S9  homolog class within the 47 proteomes (conservation/plot_panel_bins.py data)
  s10  Fig. S10 published GenEra ages joined to AFDB v6 (conservation/plot_published_suite.plot_bins data)

Output: /path/to/server_data/si_pnas/out/*.png (600 dpi for line art, 400 dpi for scatters). Original figures are not touched.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

OUT = Path("/path/to/server_data/si_pnas/out")
OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "axes.titleweight": "normal"})
GROUP7 = ["Bacteria", "Archaea", "Animalia", "Fungi", "Protist", "Parasitic Worms", "Viridiplantae"]


def letter(ax, text, x=-0.12, y=1.02, size=12):
    ax.text(x, y, text, transform=ax.transAxes, fontsize=size, fontweight="bold", va="bottom", ha="left")


def frame(ax, lw, tick_len, labelsize):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(lw)
    ax.tick_params(width=lw, length=tick_len, labelsize=labelsize)


def bottom_legend(fig, colors, names, fontsize, ncol, y=0.0, ms=None):
    ms = fontsize * 0.6 if ms is None else ms
    handles = [plt.Line2D([0], [0], marker="o", ls="", ms=ms, color=colors[g], label=g) for g in names]
    rows = -(-len(handles) // ncol)
    if rows > 1:   # matplotlib fills legends column by column; reorder so the entries read left to right, row by row
        grid = [handles[r * ncol:(r + 1) * ncol] for r in range(rows)]
        handles = [grid[r][c] for c in range(ncol) for r in range(rows) if c < len(grid[r])]
    fig.legend(handles=handles, loc="lower center", ncol=ncol, frameon=False, fontsize=fontsize, handletextpad=0.1,
               columnspacing=1.0, bbox_to_anchor=(0.5, y))


# ---------------------------------------------------------------- Fig. S1
def s1():
    sys.path.insert(0, "/path/to/server_data/afdb_v6/scripts")
    import plot_merged_7groups as pm
    import plot_plddt_density as ppd
    merged = pm.build_merged(verbose=False)
    rng = np.random.default_rng(0)
    fig, ax = plt.subplots(figsize=(6.0, 4.8))
    peak = 0
    for g in [g for g in ppd.ORDER if g in set(merged.tax_group7.dropna())]:
        vals = merged.loc[merged.tax_group7 == g, "plddt_global"].to_numpy(float)
        med = float(np.median(vals))
        dens = ppd.density(vals, rng)
        peak = max(peak, dens.max())
        ax.plot(ppd.GRID, dens, color=ppd.COLORS[g], lw=2, label=g)
        ax.axvline(med, color=ppd.COLORS[g], lw=1, ls="--", alpha=0.8)
        print(g, len(vals), round(med, 2), flush=True)
    ax.set_xlim(10, 110)
    ax.set_xticks(range(10, 101, 10))
    ax.set_ylim(0, max(0.12, peak * 1.03))     # the Archaea peak (0.131) is above 0.12; labels stop at 0.12
    ax.set_yticks(np.round(np.arange(0, 0.121, 0.02), 2))
    ax.set_xlabel("Global pLDDT (AFDB v6, per-residue mean)", fontsize=12)
    ax.set_ylabel("Density", fontsize=12)
    ax.set_title("47 model-organism proteomes + Swiss-Prot (7 groups)", fontsize=12)
    frame(ax, 1, 4, 12)
    names = [g for g in ppd.ORDER if g in set(merged.tax_group7.dropna())]
    handles = [plt.Line2D([0], [0], color=ppd.COLORS[g], lw=2, label=g) for g in names]
    ncol, rows = 4, 2
    grid = [handles[r * ncol:(r + 1) * ncol] for r in range(rows)]
    handles = [grid[r][c] for c in range(ncol) for r in range(rows) if c < len(grid[r])]   # read left to right
    fig.legend(handles=handles, loc="lower center", ncol=ncol, frameon=False, fontsize=12, handlelength=1.6,
               columnspacing=1.2, bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(rect=(0, 0.12, 1, 1))
    fig.savefig(OUT / "figS1_density_merged_7groups.png", dpi=600)
    plt.close(fig)


# ---------------------------------------------------------------- Fig. S3
LDDT_COLORS = {"Archaea": "#7f7f7f", "Bacteria": "#1f77b4", "Protist": "#9467bd", "Fungi": "#ff7f0e",
               "Animalia": "#d62728", "Parasitic Worms": "#8c564b", "Viridiplantae": "#2ca02c"}


def s3():
    df = pd.read_csv("/path/to/server_data/sifts_accuracy/lddt/lddt_primary.tsv", sep="\t")
    ok = df[df.status == "ok"].copy()
    fs = 12   # every text the size of the title
    fig, ax = plt.subplots(figsize=(6.0, 5.0))
    for g in ["Archaea", "Bacteria", "Protist", "Fungi", "Animalia", "Parasitic Worms", "Viridiplantae"]:
        s = ok[ok.tax_group7 == g]
        if s.empty:
            continue
        ax.scatter(s.plddt_mapped, s.lddt_ca, s=8.5, alpha=0.35, c=LDDT_COLORS[g], linewidths=0)
    ax.plot([0, 100], [0, 1], ls="--", c="0.5", lw=1)
    rho, _ = stats.spearmanr(ok.plddt_mapped, ok.lddt_ca, nan_policy="omit")
    ax.set_xlabel("Mean pLDDT on scored residues", fontsize=fs)
    ax.set_ylabel("lDDT-Cα", fontsize=fs)
    ax.set_xlim(10, 100)
    ax.set_xticks(range(10, 101, 10))
    ax.set_ylim(0, 1)
    ax.set_yticks(np.round(np.arange(0, 1.01, 0.1), 1))
    ax.set_title(f"AFDB v6 versus experimental structures (n = {len(ok):,})", fontsize=fs)
    frame(ax, 1, 4, fs)
    # one row of 7 names at 12 pt is wider than the page: two rows, read left to right
    bottom_legend(fig, LDDT_COLORS, ["Archaea", "Bacteria", "Protist", "Fungi", "Animalia", "Parasitic Worms",
                                     "Viridiplantae"], 10, 4, ms=3.6)
    print("S3 n", len(ok), "spearman", round(rho, 4), ok.tax_group7.value_counts().to_dict(), flush=True)
    fig.tight_layout(rect=(0, 0.1, 1, 1))
    fig.savefig(OUT / "figS3_plddt_vs_lddt.png", dpi=500)
    plt.close(fig)


# ---------------------------------------------------------------- Fig. S4
def s4():
    sys.path.insert(0, "/path/to/server_data/esm_atlas_plddt/scripts")
    import plot_precise_esm as pe
    hashes = pd.read_csv(pe.ROOT / "hashes.tsv", sep="\t",
                         dtype={"sequence_hash": str, "db_accession": str, "tax_group7": str})
    scores = pe.load_means()
    finite = scores.mean_plddt.dropna()
    scale = 100.0 if len(finite) and float(finite.max()) <= 1.5 else 1.0
    scores["esm_plddt"] = scores.mean_plddt * scale
    df = hashes.merge(scores[["sequence_hash", "esm_plddt"]], on="sequence_hash", how="left")
    plotted = df[df.esm_plddt.notna() & df.plddt_global.notna()].copy()
    rng = np.random.default_rng(1)
    drawn = plotted.iloc[rng.permutation(len(plotted))]
    fig, ax = plt.subplots(figsize=(6.0, 5.0))
    ax.scatter(drawn.esm_plddt, drawn.plddt_global, s=3, c=drawn.tax_group7.map(pe.PALETTE).fillna("#666666"),
               alpha=0.4, linewidths=0, rasterized=True)
    x, y = plotted.esm_plddt.to_numpy(float), plotted.plddt_global.to_numpy(float)
    ax.plot([0, 100], [0, 100], color="#b0b0b0", lw=1.2, zorder=3)
    rho, _ = stats.spearmanr(x, y)
    r, _ = stats.pearsonr(x, y)
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.set_xticks(range(0, 101, 20))
    ax.set_yticks(range(0, 101, 20))
    frame(ax, 1, 4, 12)
    ax.set_xlabel("ESMFold mean pLDDT", fontsize=12)
    ax.set_ylabel("Global pLDDT (AFDB v6)", fontsize=12)
    ax.set_title(f"n = {len(plotted):,}; Spearman ρ = {rho:.3f}, Pearson r = {r:.3f}", fontsize=12)
    bottom_legend(fig, pe.PALETTE, [g for g in pe.DRAW_ORDER if (plotted.tax_group7 == g).any()], 10, 4, ms=3.6)
    print("S4 n", len(plotted), "rho", round(rho, 4), "r", round(r, 4),
          plotted.tax_group7.value_counts().to_dict(), flush=True)
    fig.tight_layout(rect=(0, 0.1, 1, 1))
    fig.savefig(OUT / "figS4_af2_vs_esmfold.png", dpi=500)
    plt.close(fig)


# ---------------------------------------------------------------- Fig. S5
PAL = {"Viridiplantae": "#0eb519", "Animalia": "#E41E26", "Fungi": "#F68A21", "Bacteria": "#3D5BA9",
       "Archaea": "#7f7f7f", "Protist": "#9467bd", "Parasitic Worms": "#8c564b"}
ALPHA = {"Bacteria": 0.08, "Animalia": 0.08, "Archaea": 0.2, "Fungi": 0.12, "Protist": 0.2,
         "Parasitic Worms": 0.15, "Viridiplantae": 0.35}


def scatter_fit(ax, df, xcol, xlabel, xlim, xticks, title):
    rng = np.random.default_rng(0)
    for g in GROUP7:
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
    ax.set_ylabel("Global pLDDT (AFDB v6)", fontsize=8)
    ax.set_title(f"{title}\nSpearman ρ = {rho:.3f}\nPearson r = {r:.3f}", fontsize=7.5)
    print(title, "n", len(x), "rho", round(rho, 4), "r", round(r, 4), "slope", res.slope,
          df.tax_group7.value_counts().to_dict(), flush=True)


def s5():
    sys.path.insert(0, "/path/to/server_data/afdb_v6/scripts")
    import plot_merged_7groups as pm
    depth = pm.build_merged(verbose=False)
    depth = depth[depth.msa_raw_depth.notna()].copy()
    neff = pd.read_parquet("/path/to/server_data/neff_v6/merged_neff.parquet")
    neff = neff[neff.neff_per_len.notna()].copy()
    neff["log10_neff"] = np.log10(neff.neff + 1)
    fig, axes = plt.subplots(1, 3, figsize=(6.0, 2.45))
    scatter_fit(axes[0], depth, "msa_raw_depth", "MSA depth (sequences)", (0, 20000),
                [0, 10000, 20000], "AFDB MSA depth")
    scatter_fit(axes[1], neff, "log10_neff", "log10(Neff + 1)", (0, float(neff.log10_neff.max()) + 0.1),
                [0, 1, 2, 3, 4], "Neff (OpenProteinSet)")
    scatter_fit(axes[2], neff, "neff_per_len", "Neff / length", (0, 60), [0, 20, 40, 60],
                "Neff per residue")
    for a, t in zip(axes, "ABC"):
        letter(a, t, x=-0.2, y=1.04, size=9)
    bottom_legend(fig, PAL, GROUP7, 7, 7)
    fig.tight_layout(rect=(0, 0.08, 1, 1), w_pad=1.0)
    fig.savefig(OUT / "figS5_msa_AB.png", dpi=500)
    plt.close(fig)


# ---------------------------------------------------------------- Fig. S9
P_GROUPS = ["Viridiplantae", "Animalia", "Fungi"]
P_LABEL = {"Viridiplantae": "Plants", "Animalia": "Animals", "Fungi": "Fungi"}
P_PAL = {"Viridiplantae": "#0eb519", "Animalia": "#E41E26", "Fungi": "#F68A21"}
P_BINS = ["universal", "eukaryote_wide", "kingdom_only", "family_or_genus", "species_only"]
P_BIN_LABEL = {"universal": "All life", "eukaryote_wide": "Eukaryotes", "kingdom_only": "Kingdom",
               "family_or_genus": "Family", "species_only": "Species"}
P_BIN_COLOR = {"universal": "#1b4f72", "eukaryote_wide": "#2e86ab", "kingdom_only": "#7f8c8d",
               "family_or_genus": "#f5b041", "species_only": "#922b21"}


def s9():
    df = pd.read_parquet("/path/to/server_data/conservation/panel_strata.parquet")
    df = df[df.tax_group7.isin(P_GROUPS)].copy()
    df["log10_neff"] = np.log10(df.neff.clip(lower=0) + 1)
    fig, axes = plt.subplots(1, 3, figsize=(6.5, 3.0), gridspec_kw={"width_ratios": [1.05, 1.3, 1.3]})
    ax = axes[0]
    bottom = np.zeros(len(P_GROUPS))
    for b in P_BINS:
        h = [float((df[df.tax_group7 == g].panel_bin == b).mean()) for g in P_GROUPS]
        ax.bar(range(len(P_GROUPS)), h, bottom=bottom, color=P_BIN_COLOR[b], width=0.72, label=P_BIN_LABEL[b])
        bottom += np.array(h)
    ax.set_xticks(range(len(P_GROUPS)))
    ax.set_xticklabels([P_LABEL[g] for g in P_GROUPS], rotation=30, ha="right")
    ax.set_ylim(0, 1)
    ax.set_ylabel("Fraction of proteins", fontsize=7.5)
    frame(ax, 1.0, 2.5, 7)
    for ax, col, ylabel, ylim in ((axes[1], "plddt_global", "Global pLDDT", (20, 100)),
                                  (axes[2], "log10_neff", "log10(Neff + 1)", (0, 4.2))):
        width = 0.24
        for gi, g in enumerate(P_GROUPS):
            data, pos = [], []
            for bi, b in enumerate(P_BINS):
                v = df.loc[(df.tax_group7 == g) & (df.panel_bin == b), col].dropna().to_numpy()
                if len(v) < 20:
                    continue
                data.append(v)
                pos.append(bi + (gi - 1) * width)
            bp = ax.boxplot(data, positions=pos, widths=width * 0.9, showfliers=False, patch_artist=True,
                            medianprops={"color": "black", "linewidth": 0.8}, whiskerprops={"linewidth": 0.6},
                            capprops={"linewidth": 0.6}, boxprops={"linewidth": 0.4})
            for box in bp["boxes"]:
                box.set_facecolor(P_PAL[g])
                box.set_alpha(0.9)
        ax.set_xticks(range(len(P_BINS)))
        ax.set_xticklabels([P_BIN_LABEL[b] for b in P_BINS], rotation=30, ha="right")
        ax.set_ylabel(ylabel, fontsize=7.5)
        ax.set_ylim(*ylim)
        frame(ax, 1.0, 2.5, 7)
    fig.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=P_BIN_COLOR[b], label=P_BIN_LABEL[b]) for b in P_BINS],
               frameon=False, fontsize=7, loc="lower left", ncol=5, bbox_to_anchor=(0.01, 0.0),
               title="A: homolog class", title_fontsize=7, handlelength=1, columnspacing=0.8)
    fig.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=P_PAL[g], label=P_LABEL[g]) for g in P_GROUPS],
               frameon=False, fontsize=7, loc="lower right", ncol=3, bbox_to_anchor=(0.99, 0.0),
               title="B and C: kingdom", title_fontsize=7, handlelength=1, columnspacing=0.8)
    for a, t in zip(axes, "ABC"):
        letter(a, t, x=-0.3 if t == "A" else -0.22, y=1.03, size=11)
    fig.tight_layout(rect=(0, 0.13, 1, 1), w_pad=1.2)
    fig.savefig(OUT / "figS9_homolog_class.png", dpi=600)
    plt.close(fig)
    for g in P_GROUPS:
        s = df[df.tax_group7 == g]
        print(g, len(s), s.panel_bin.value_counts(normalize=True).round(3).to_dict(),
              s.groupby("panel_bin").plddt_global.median().round(1).to_dict(),
              s.groupby("panel_bin").neff.median().round(0).to_dict(), flush=True)


# ---------------------------------------------------------------- Fig. S10
def s10():
    sys.path.insert(0, "/path/to/server_data/conservation/scripts")
    import plot_published_suite as ps
    parent, rank, names = ps.mga.load_taxdump(ps.TAXDUMP)
    meta = ps.load_meta()
    sgd = ps.load_sgd()
    ages = ps.assign_names(ps.barrera_genes(sgd), parent, rank, names)
    parts = []
    for proteome, *_rest in ps.BARRERA:
        hit = ps.join_meta(ages[ages.proteome.eq(proteome)], meta, proteome)
        parts.append(hit.sort_values("rank").drop_duplicates("db_accession"))
    parts.append(ps.join_meta(ps.load_human(), meta, "HUMAN"))
    df = pd.concat(parts, ignore_index=True)
    plot = df[df.age_bin.isin(ps.BINS)].copy()
    order = [row[0] for row in ps.BARRERA if row[0] in set(plot.proteome)]
    fig, axes = plt.subplots(2, 1, figsize=(9.0, 9.6))
    ax = axes[0]
    bottom = np.zeros(len(order))
    for b in ps.BINS:
        heights = [(plot[plot.proteome.eq(p)].age_bin.eq(b)).mean() for p in order]
        ax.bar(np.arange(len(order)), heights, bottom=bottom, color=ps.BIN_COLOR[b], width=0.78, label=ps.BIN_LABEL[b])
        bottom += np.array(heights)
    labels = [f"{label}\nn={int(df.proteome.eq(p).sum()):,}" for p, _k, _s, _t, label, _g in ps.BARRERA if p in order]
    ax.set_xticks(np.arange(len(order)))
    ax.set_xticklabels(labels, fontsize=10, rotation=45, ha="right")
    ax.set_ylim(0, 1)
    ax.set_ylabel("Fraction of matched proteins", fontsize=11)
    frame(ax, 1.2, 3, 10)
    ax.legend(frameon=False, fontsize=10, loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=5,
              title="Age class", title_fontsize=10, columnspacing=1.2, handlelength=1.4)
    ax = axes[1]
    step, width = 1.7, 0.15
    for gi, proteome in enumerate(order):
        data, pos = [], []
        for bi, b in enumerate(ps.BINS):
            v = plot.loc[(plot.proteome.eq(proteome)) & (plot.age_bin.eq(b)), "plddt_global"].dropna()
            if len(v) < 20:
                continue
            data.append(v.to_numpy())
            pos.append(bi * step + (gi - (len(order) - 1) / 2) * width * 1.15)
        if not data:
            continue
        bp = ax.boxplot(data, positions=pos, widths=width, showfliers=False, patch_artist=True,
                        medianprops={"color": "black", "linewidth": 0.8}, boxprops={"linewidth": 0.4},
                        whiskerprops={"linewidth": 0.6}, capprops={"linewidth": 0.6})
        for box in bp["boxes"]:
            box.set_facecolor(ps.SPECIES_COLOR[proteome])
            box.set_alpha(0.85)
    ax.set_xticks([bi * step for bi in range(len(ps.BINS))])
    ax.set_xticklabels([ps.BIN_LABEL[b] for b in ps.BINS], fontsize=10)
    ax.set_xlim(-0.8, (len(ps.BINS) - 1) * step + 0.8)
    ax.set_ylabel("Global pLDDT", fontsize=11)
    ax.set_xlabel("Age class", fontsize=11)
    ax.set_ylim(20, 100)
    frame(ax, 1.2, 3, 10)
    handles = [plt.Rectangle((0, 0), 1, 1, color=ps.SPECIES_COLOR[p], label=lab)
               for p, _k, _s, _t, lab, _g in ps.BARRERA if p in order]
    ax.legend(handles=handles, frameon=False, fontsize=9.5, loc="lower center", bbox_to_anchor=(0.5, 1.01),
              ncol=len(handles), title="Species", title_fontsize=10, columnspacing=0.9, handlelength=1.2,
              handletextpad=0.4)
    letter(axes[0], "A", x=-0.08, y=1.16, size=16)
    letter(axes[1], "B", x=-0.08, y=1.16, size=16)
    fig.tight_layout(h_pad=3.0)
    fig.savefig(OUT / "figS10_published_ages_2x1.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    for p in order:
        s = plot[plot.proteome.eq(p)]
        print(p, int(df.proteome.eq(p).sum()), s.age_bin.value_counts(normalize=True).round(3).to_dict(),
              s.groupby("age_bin").plddt_global.median().round(1).to_dict(), flush=True)


if __name__ == "__main__":
    for job in sys.argv[1:] or ["s1", "s3", "s4", "s5", "s9", "s10"]:
        {"s1": s1, "s3": s3, "s4": s4, "s5": s5, "s9": s9, "s10": s10}[job]()
        print("DONE", job, flush=True)
