#!/usr/bin/env python3
"""Local figures for the PNAS-formatted concise SI (data already on B:).

  ages      2 x 2: A sel3800 and B sel690 global pLDDT versus log10(Neff + 1) by GenEra age bin; C age composition
            and D global pLDDT by age bin (median, interquartile range) for sel690, sel3800 and the plant proteins
            solved in 2021-2026 (pdb3825, full-length). One age-bin legend for all panels, across the top.
  crops     5 x 2: continuous measures for Arabidopsis, rice, maize and soybean (compare_sets.plot_panels data),
            panels A-J with x and y axis labels
  msa       1 x 2: the two existing MSA-depth scatters (AFDB MSA depth; Neff) side by side, panels A-B
  disorder  2 x 1 stack of the two existing disorder figures (D1 over D2), panels A-D

Usage: python si_pnas_local.py [ages crops msa disorder]   -> SI_2026/pnas_figures/parts/*.png
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont

B = Path(r"/path/to/project")
OUT = B / "SI_2026" / "pnas_figures" / "parts"
OUT.mkdir(parents=True, exist_ok=True)
FIG = B / "template_esm" / "campaign" / "figures"
PC = B / "SI_2026" / "pdb3825_vs_selections" / "per_protein_comparison.tsv"
IMG = Path("/path/to/scratch/docx_build/img_cache")

BINS = ["all_life", "eukaryotes", "kingdom_only", "family_or_genus"]                 # as plot_sel690_ages.py
BIN_LABEL = {"all_life": "All life", "eukaryotes": "Eukaryotes", "kingdom_only": "Kingdom",
             "family_or_genus": "Family or genus"}
BIN_COLOR = {"all_life": "#1b4f72", "eukaryotes": "#2e86ab", "kingdom_only": "#7f8c8d", "family_or_genus": "#f5b041"}


def style(ax, lw=1.2):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(lw)
    ax.tick_params(width=lw, length=3)


def letter(ax, text, x=-0.12, y=1.03, size=18):
    ax.text(x, y, text, transform=ax.transAxes, fontsize=size, fontweight="bold", va="bottom", ha="left")


# ---------------------------------------------------------------- 2 x 2 gene age
def ages():
    s3800 = pd.read_csv(FIG / "genera_ages_3800" / "sel3800_genera_ages.csv")
    s690 = pd.read_csv(FIG / "genera_ages_690" / "sel690_genera_ages.csv")
    pc = pd.read_csv(PC, sep="\t", low_memory=False)
    pdb = pc[pc.set.eq("pdb3825_full")].rename(columns={"plddt": "plddt_global"})
    plt.rcParams.update({"font.size": 12, "axes.labelsize": 13, "xtick.labelsize": 11.5, "ytick.labelsize": 11.5})
    fig, axes = plt.subplots(2, 2, figsize=(14.4, 11.6))
    counts = {}
    for ax, (tag, df) in zip(axes[0], (("sel3800", s3800), ("sel690", s690))):
        use = df[df.age_bin.isin(BINS) & df.neff.notna() & df.plddt_global.notna()].copy()
        use["log10_neff"] = np.log10(use.neff.clip(lower=0) + 1)
        counts[tag] = {b: int(use.age_bin.eq(b).sum()) for b in BINS}
        for b in BINS:
            piece = use[use.age_bin.eq(b)]
            if len(piece) < 8:
                continue
            ax.scatter(piece.log10_neff, piece.plddt_global, s=22, alpha=0.7, c=BIN_COLOR[b], linewidths=0)
        ax.set_xlabel("log10(Neff + 1)")
        ax.set_ylabel("Global pLDDT")
        ax.set_ylim(20, 100)
        ax.set_title(f"AI {tag}: GenEra age and Neff (n = {len(use):,})", fontsize=13)
        style(ax)
    print("A/B counts", counts, flush=True)

    sets = [("AI sel690", s690[s690.age_bin.isin(BINS)], "o"), ("AI sel3800", s3800[s3800.age_bin.isin(BINS)], "s"),
            ("PDB 2021–26", pdb[pdb.age_bin.isin(BINS)], "^")]
    ax = axes[1, 0]
    bottom = np.zeros(len(sets))
    for b in BINS:
        heights = [float(sub.age_bin.eq(b).mean()) if len(sub) else 0.0 for _, sub, _ in sets]
        ax.bar(range(len(sets)), heights, bottom=bottom, color=BIN_COLOR[b], width=0.62)
        for i, h in enumerate(heights):
            if h >= 0.06:
                ax.text(i, bottom[i] + h / 2, f"{h:.0%}", ha="center", va="center", fontsize=11, color="white")
        bottom += np.asarray(heights)
    ax.set_xticks(range(len(sets)))
    ax.set_xticklabels([f"{name}\nn = {len(sub):,}" for name, sub, _ in sets])
    ax.set_ylim(0, 1)
    ax.set_ylabel("Fraction of proteins")
    ax.set_title("GenEra age composition", fontsize=13)
    style(ax)

    ax = axes[1, 1]
    for j, (name, sub, marker) in enumerate(sets):
        sub = sub[sub.plddt_global.notna()]
        for i, b in enumerate(BINS):
            vals = sub.loc[sub.age_bin.eq(b), "plddt_global"]
            if len(vals) < 8:
                continue
            xpos = i + (j - 1) * 0.25
            ax.vlines(xpos, vals.quantile(0.25), vals.quantile(0.75), color=BIN_COLOR[b], lw=2.6, alpha=0.75)
            ax.scatter([xpos], [vals.median()], s=70, color=BIN_COLOR[b], marker=marker, edgecolors="black",
                       linewidths=0.6, zorder=3)
        ax.scatter([], [], s=70, color="white", edgecolors="black", marker=marker, label=name)
    ax.set_xticks(range(len(BINS)))
    ax.set_xticklabels([BIN_LABEL[b] for b in BINS])
    ax.set_ylim(20, 100)
    ax.set_ylabel("Global pLDDT (median and IQR)")
    ax.set_title("AlphaFold2 confidence by age", fontsize=13)
    ax.legend(frameon=False, fontsize=12, loc="lower left")
    style(ax)
    for a, t in zip(axes.flat, "ABCD"):
        letter(a, t)
    handles = [plt.Line2D([0], [0], marker="o", ls="", ms=12, color=BIN_COLOR[b], label=BIN_LABEL[b]) for b in BINS]
    fig.legend(handles=handles, loc="upper center", ncol=4, frameon=False, fontsize=16,   # 2 x the original 8 pt
               bbox_to_anchor=(0.5, 1.0), columnspacing=2.0, handletextpad=0.4)
    fig.tight_layout(rect=(0, 0, 1, 0.955), h_pad=2.5, w_pad=3)
    fig.savefig(OUT / "ages_ABCD.png", dpi=300)
    plt.close(fig)
    summary = {name: {b: int(sub.age_bin.eq(b).sum()) for b in BINS} for name, sub, _ in sets}
    plddt_n = {name: int(sub.plddt_global.notna().sum()) for name, sub, _ in sets}
    print("C composition", summary, "D with pLDDT", plddt_n, flush=True)
    pd.DataFrame(summary).to_csv(OUT / "ages_ABCD_counts.tsv", sep="\t")


# ---------------------------------------------------------------- 5 x 2 continuous measures
SETS = ["sel3800", "sel690", "pdb3825_full", "pdb3825_construct", "random_matched", "random_pool"]
SHORT = {"sel3800": "AI\nsel3800", "sel690": "AI\nsel690", "pdb3825_full": "PDB full-\nlength",
         "pdb3825_construct": "PDB\nconstruct", "random_matched": "Random,\nmatched", "random_pool": "Random,\nall"}
COLOR = {"sel3800": "#d95f02", "sel690": "#e6ab02", "pdb3825_full": "#1b9e77", "pdb3825_construct": "#66a61e",
         "random_matched": "#7570b3", "random_pool": "#999999"}
CONT = [("frac_disordered", "Disordered fraction\n(AIUPred > 0.5)"), ("mean_disorder", "Mean AIUPred score"),
        ("longest_idr", "Longest disordered\nregion (residues)"), ("length", "Length (residues)"),
        ("n_plant_homologs", "log10(plant homologs + 1)\nboth ≥ 50% covered"),
        ("n_plant_homologs_qcov", "log10(plant homologs + 1)\nquery ≥ 50% covered"),
        ("n_plant_homologs_id40", "log10(plant homologs + 1)\n≥ 40% identity"),
        ("n_crop4_homologs", "log10(homologs + 1) in\nthe four crop proteomes"),
        ("n_within_1.5", "log10(plant proteins + 1)\nwithin ESM-2 distance 1.5"),
        ("d_homolog_mean", "Average homolog distance")]


def crops():
    df = pd.read_csv(PC, sep="\t", low_memory=False)
    df = df[df.crop.notna()]
    plt.rcParams.update({"font.size": 11.8, "axes.labelsize": 12.3, "xtick.labelsize": 11.2, "ytick.labelsize": 11.2})
    fig, axes = plt.subplots(5, 2, figsize=(10.8, 13.2))
    rows = []
    for k, (ax, (c, ylabel)) in enumerate(zip(axes.flat, CONT)):
        sets = [s for s in SETS if s in set(df.set)]
        data = [df.loc[df.set.eq(s), c].dropna().to_numpy(float) for s in sets]
        log = c.startswith("n_")
        plot_data = [np.log10(d + 1) if log else d for d in data]
        keep = [i for i, d in enumerate(plot_data) if len(d)]
        bp = ax.boxplot([plot_data[i] for i in keep], showfliers=False, patch_artist=True, widths=0.6,
                        medianprops={"color": "black", "linewidth": 1.2})
        for patch, i in zip(bp["boxes"], keep):
            patch.set_facecolor(COLOR[sets[i]])
            patch.set_alpha(0.8)
        ax.set_xticks(range(1, len(keep) + 1))
        ax.set_xticklabels([SHORT[sets[i]] for i in keep], fontsize=10.2)
        ax.set_xlabel("Protein set")
        ax.set_ylabel(ylabel)
        top = ax.get_ylim()[1]
        for j, i in enumerate(keep):
            ax.text(j + 1, top, f"{len(data[i]):,}", ha="center", va="bottom", fontsize=9.1)
        style(ax)
        letter(ax, "ABCDEFGHIJ"[k], x=-0.24, y=1.07, size=17)
        rows.append({"panel": "ABCDEFGHIJ"[k], "measure": c, **{SHORT[sets[i]].replace("\n", " "): len(data[i]) for i in keep}})
    fig.tight_layout(h_pad=1.6, w_pad=2.5)
    fig.savefig(OUT / "crops_AJ.png", dpi=300)
    plt.close(fig)
    print(pd.DataFrame(rows).to_string(index=False), flush=True)


# ---------------------------------------------------------------- Fig. S2 (AFDB v4 versus v6), from the saved outputs
V4V6 = B / "SI_2026" / "plddt_v4_vs_v6"
V4V6_ORDER = ["Archaea", "Bacteria", "Protist", "Fungi", "Animalia", "Parasitic worms", "Viridiplantae"]
V4V6_COLORS = {"Archaea": "#8c564b", "Bacteria": "#1f77b4", "Protist": "#9467bd", "Fungi": "#ff7f0e",
               "Animalia": "#d62728", "Parasitic worms": "#7f7f7f", "Viridiplantae": "#2ca02c"}


def v4v6():
    """same panels as template_esm/campaign/v4v6/plddt_v4_vs_v6.py, uppercase bold panel letters"""
    from matplotlib.colors import LogNorm
    df = pd.read_csv(V4V6 / "plddt_v4_vs_v6_entries.tsv.gz", sep="\t")
    summary = pd.read_csv(V4V6 / "plddt_v4_vs_v6_summary.tsv", sep="\t")
    a = summary.iloc[0]
    plt.rcParams.update({"font.size": 7, "axes.labelsize": 7.5, "xtick.labelsize": 7, "ytick.labelsize": 7})
    fig, axes = plt.subplots(2, 2, figsize=(6.5, 5.6))
    ax = axes[0][0]
    h2 = ax.hist2d(df.plddt_v4, df.plddt_v6, bins=200, range=[[20, 100], [20, 100]], cmap="viridis", norm=LogNorm(), cmin=1)
    ax.plot([20, 100], [20, 100], color="black", lw=0.8, ls="--")
    fig.colorbar(h2[3], ax=ax, fraction=0.046, pad=0.03, label="Entries per bin")
    ax.set_xlabel("Global pLDDT, AFDB v4")
    ax.set_ylabel("Global pLDDT, AFDB v6")
    ax.set_title("Same entries in AFDB v4 and v6", fontsize=7.5)
    # agreement statistics are given in the figure legend (caption), not in the panel
    ax = axes[0][1]
    mean = (df.plddt_v4 + df.plddt_v6) / 2
    h2 = ax.hist2d(mean, df["diff"], bins=200, range=[[20, 100], [-30, 30]], cmap="viridis", norm=LogNorm(), cmin=1)
    for yv, st in ((a.mean_diff_v6_minus_v4, "-"), (a.loa_low, "--"), (a.loa_high, "--")):
        ax.axhline(yv, color="black", lw=0.8, ls=st)
    fig.colorbar(h2[3], ax=ax, fraction=0.046, pad=0.03, label="Entries per bin")
    ax.set_xlabel("Mean of v4 and v6 pLDDT")
    ax.set_ylabel("v6 − v4 pLDDT")
    ax.set_title("Bland–Altman: difference versus mean", fontsize=7.5)
    ax = axes[1][0]
    groups = [g for g in V4V6_ORDER if (df.tax_group7.eq(g) & df.same_sequence).sum() > 20]
    bp = ax.boxplot([df.loc[df.tax_group7.eq(g), "diff"] for g in groups], showfliers=False, patch_artist=True,
                    medianprops={"color": "black", "linewidth": 0.8})
    for patch, g in zip(bp["boxes"], groups):
        patch.set_facecolor(V4V6_COLORS[g])
        patch.set_alpha(0.7)
    ax.axhline(0, color="black", lw=0.6)
    ax.set_xticks(range(1, len(groups) + 1))
    ax.set_xticklabels(groups, rotation=25, ha="right")
    ax.set_ylabel("v6 − v4 pLDDT")
    ax.set_title("Change by taxonomic group", fontsize=7.5)
    ax = axes[1][1]
    sub = summary[summary.subset.str.startswith(("group:", "plant:"))].reset_index(drop=True)
    labels = sub.subset.str.replace("group: ", "").str.replace("plant: ", "plant: ")
    ypos = np.arange(len(sub))[::-1]
    cols = ["#2ca02c" if s.startswith("plant") or "Viridiplantae" in s else "#555555" for s in sub.subset]
    ax.scatter(sub.spearman_rho, ypos, c=cols, s=18, zorder=3, label="Spearman ρ")
    ax.scatter(sub.lin_ccc, ypos, facecolors="none", edgecolors=cols, s=18, zorder=3, label="Lin's CCC")
    for yy, nn in zip(ypos, sub.n):
        ax.text(1.0005, yy, f"n={int(nn):,}", va="center", fontsize=5.5)
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels, fontsize=6)
    ax.set_xlim(max(0.97, min(sub.spearman_rho.min(), sub.lin_ccc.min()) - 0.002), 1.0)
    ax.legend(frameon=False, fontsize=6, loc="lower left")
    ax.set_xlabel("Agreement, v4 versus v6")
    ax.set_title("Agreement by group and plant proteome", fontsize=7.5)
    for axx, t in zip(axes.flat, "ABCD"):
        axx.spines[["top", "right"]].set_visible(False)
        letter(axx, t, x=-0.2 if t != "D" else -0.62, y=1.04, size=10)
    fig.tight_layout(h_pad=1.5, w_pad=1.2)
    fig.savefig(OUT / "figS2_v4_vs_v6.png", dpi=500)
    plt.close(fig)
    print("S2", {k: round(float(a[k]), 4) for k in ("n", "pearson_r", "spearman_rho", "lin_ccc", "mean_diff_v6_minus_v4",
                                                   "loa_low", "loa_high", "frac_abs_diff_gt5", "frac_cross_70")})
    print(sub[["subset", "n", "spearman_rho", "lin_ccc", "mean_diff_v6_minus_v4"]].round(4).to_string(index=False))


# ---------------------------------------------------------------- Fig. S14 (UMAP / HDBSCAN maps), letters 2 pt larger
def umap():
    sys.path.insert(0, str(B / "SI_2026" / "scripts"))
    import redraw_hdbscan_maps as rh
    tables = {key: pd.read_csv(rh.DR / rel, sep="\t") for key, rel, _ in rh.MAPS}
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 8.6))
    for ax, (key, _, title), tag in zip(axes.ravel(), rh.MAPS, "ABCD"):
        rh.draw(ax, tables[key], title)
        ax.text(-0.12, 1.04, tag, transform=ax.transAxes, fontsize=15, fontweight="bold")   # was 13 pt
    fig.tight_layout()
    fig.savefig(OUT / "umap_ABCD.png", dpi=200)
    plt.close(fig)


# ---------------------------------------------------------------- image composites
FONT = ImageFont.truetype("arialbd.ttf", 40)     # half the earlier 80 px (Fig. S8 letters)


def labelled(im, letters, strip=55):
    """add a white strip above the image and bold panel letters at the given fractional x positions"""
    out = Image.new("RGB", (im.width, im.height + strip), "white")
    out.paste(im, (0, strip))
    d = ImageDraw.Draw(out)
    for text, fx in letters:
        d.text((int(fx * im.width) + 8, 8), text, fill="black", font=FONT)
    return out


def msa():
    a = Image.open(IMG / "figures__plddt_vs_msa_depth_scatter_linear.jpg").convert("RGB")
    b = Image.open(IMG / "report_assets__plddt_vs_msa_depth_scatter_log10_neff.jpg").convert("RGB")
    b = b.resize((int(b.width * a.height / b.height), a.height), Image.LANCZOS)
    a, b = labelled(a, [("A", 0)]), labelled(b, [("B", 0)])
    gap = 60
    out = Image.new("RGB", (a.width + gap + b.width, a.height), "white")
    out.paste(a, (0, 0))
    out.paste(b, (a.width + gap, 0))
    out.save(OUT / "msa_AB.png", dpi=(300, 300))
    print("msa", out.size, flush=True)


def disorder():
    a = Image.open(IMG / "report_assets__D1_plddt_vs_disorder.png").convert("RGB")
    b = Image.open(IMG / "report_assets__D2_all_vs_ordered_plddt.png").convert("RGB")
    b = b.resize((a.width, int(b.height * a.width / b.width)), Image.LANCZOS)
    a = labelled(a, [("A", 0.0), ("B", 0.575)])
    b = labelled(b, [("C", 0.0), ("D", 0.575)])
    out = Image.new("RGB", (a.width, a.height + 40 + b.height), "white")
    out.paste(a, (0, 0))
    out.paste(b, (0, a.height + 40))
    out.save(OUT / "disorder_ABCD.png", dpi=(300, 300))
    print("disorder", out.size, flush=True)


if __name__ == "__main__":
    for job in sys.argv[1:] or ["ages", "crops", "msa", "disorder"]:
        {"ages": ages, "crops": crops, "msa": msa, "disorder": disorder, "v4v6": v4v6, "umap": umap}[job]()
