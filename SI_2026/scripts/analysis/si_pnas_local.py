#!/usr/bin/env python3
"""Local figures for the PNAS-formatted concise SI (data already on B:).

  ages      2 x 2: A sel3800 and B sel690 global pLDDT versus log10(AF2-style Neff + 1) by GenEra age bin; C age composition
            and D global pLDDT by age bin (median, interquartile range) for sel690, sel3800 and the plant proteins
            solved in 2021-2026 (pdb3825, full-length). One age-bin legend for all panels, across the top.
  crops     5 x 2: continuous measures for Arabidopsis, rice, maize and soybean (compare_sets.plot_panels data),
            panels A-J with x and y axis labels
  msa       1 x 2: the two existing MSA-depth scatters (AFDB MSA depth; Neff) side by side, panels A-B
  disorder  2 x 1 stack of the two existing disorder figures (D1 over D2), panels A-D
  annotation  A UniProt annotation score, B annotation indicators, C Pfam families per protein and D the most
            frequent best Pfam families of sel3800, sel690 and the full-length pdb3825 proteins (needs
            parts/annotation_pfam_names.tsv, the InterPro names of the Pfam families shown, written by
            annotation_pfam_names.py)

Usage: python si_pnas_local.py [ages crops msa disorder annotation]   -> SI_2026/pnas_figures/parts/*.png
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
FIGS20 = B / "SI_2026" / "datasets" / "Dataset_S17_source_data_SI_figures" / "FigS20"
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
    s3800 = pd.read_csv(FIGS20 / "sel3800_genera_ages.csv")
    s690 = pd.read_csv(FIGS20 / "sel690_genera_ages.csv")
    pc = pd.read_csv(PC, sep="\t", low_memory=False)
    pdb = pc[pc.set.eq("pdb3825_full")].rename(columns={"plddt": "plddt_global"})
    plt.rcParams.update({"font.size": 12, "axes.labelsize": 13, "xtick.labelsize": 11.5, "ytick.labelsize": 11.5})
    fig, axes = plt.subplots(2, 2, figsize=(14.4, 11.6))
    counts = {}
    for ax, (tag, df) in zip(axes[0], (("sel3800", s3800), ("sel690", s690))):
        use = df[df.age_bin.isin(BINS) & df["AF2-style_Neff"].notna() & df.plddt_global.notna()].copy()
        use["log10_neff"] = np.log10(use["AF2-style_Neff"].clip(lower=0) + 1)
        counts[tag] = {b: int(use.age_bin.eq(b).sum()) for b in BINS}
        for b in BINS:
            piece = use[use.age_bin.eq(b)]
            if len(piece) < 8:
                continue
            ax.scatter(piece.log10_neff, piece.plddt_global, s=22, alpha=0.7, c=BIN_COLOR[b], linewidths=0)
        ax.set_xlabel("log10(AF2-style Neff + 1)")
        ax.set_ylabel("Global pLDDT")
        ax.set_ylim(20, 100)
        ax.set_title(f"AI {tag}: GenEra age and AF2-style Neff (n = {len(use):,})", fontsize=13)
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


# ---------------------------------------------------------------- annotation and characterization of the three sets
ANN_SETS = ["sel3800", "sel690", "pdb3825_full"]
ANN_SHORT = {"sel3800": "AI\nsel3800", "sel690": "AI\nsel690", "pdb3825_full": "PDB\n2021–26"}
ANN_LABEL = {"sel3800": "AI sel3800", "sel690": "AI sel690", "pdb3825_full": "PDB 2021–26"}
SCORE_COLOR = ["#c6dbef", "#9ecae1", "#6baed6", "#3182bd", "#08519c"]            # UniProt annotation score 1-5
NPF_CATS, NPF_COLOR = ["0", "1", "2", "≥ 3"], ["#d9d9d9", "#c7e9c0", "#74c476", "#238b45"]
INDICATORS = [("is_reviewed", "Swiss-Prot reviewed"), ("has_function_cc", "Function comment"),
              ("has_interpro", "InterPro entry"), ("has_pfam", "≥ 1 Pfam domain"), ("has_duf", "DUF domain"),
              ("uncharacterized", "Uncharacterized"), ("poorly_characterized", "Poorly characterized")]
N_TOP = 15                       # panel D: up to 15 most frequent best Pfam families per set, never splitting a tie


def top_families(counts, max_n=N_TOP):
    """the most frequent families (value_counts, descending), up to max_n, cut where no tie is split"""
    n = max_n
    while 0 < n < len(counts) and counts.iloc[n - 1] == counts.iloc[n]:
        n -= 1
    return counts.iloc[:n]


def dark(color, alpha=1.0):
    """True when white text reads better than black on this colour (alpha-blended on white)"""
    r, g, b = (alpha * int(color[i:i + 2], 16) / 255 + 1 - alpha for i in (1, 3, 5))
    return 0.2126 * r + 0.7152 * g + 0.0722 * b < 0.5


def wrap(text, limit=40):
    """break a long label once, at the space or slash nearest its middle"""
    if len(text) <= limit:
        return text
    i = min((i for i, ch in enumerate(text) if ch in " /"), key=lambda i: abs(i - len(text) / 2))
    return text[:i + 1] + "\n" + text[i + 1:] if text[i] == "/" else text[:i] + "\n" + text[i + 1:]


def annotation():
    """row 1: A UniProt annotation score, B Pfam families per protein, C annotation features (with the set legend);
    row 2: D the most frequent best Pfam families of each set; for sel3800, sel690 and the full-length pdb3825
    proteins (drawn at print size: the 6.0 in text width of the PNAS SI, and short enough to share the page with its
    legend)"""
    pc = pd.read_csv(PC, sep="\t", low_memory=False)
    pf_name = pd.read_csv(OUT / "annotation_pfam_names.tsv", sep="\t").set_index("pfam_name").interpro_name
    d = {s: pc[pc.set.eq(s)] for s in ANN_SETS}
    n = {s: len(g) for s, g in d.items()}
    W, H = 6.0, 5.9
    plt.rcParams.update({"font.size": 6.5, "axes.labelsize": 7, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5})
    fig = plt.figure(figsize=(W, H))

    def box(x, y, w, h):                                           # axes position in inches
        return fig.add_axes([x / W, y / H, w / W, h / H])

    def heading(letter, x, y, title, dx=0.18):                    # bold panel letter, title right after it
        fig.text(x / W, y / H, letter, fontsize=10, fontweight="bold", ha="left", va="baseline")
        fig.text((x + dx) / W, y / H, title, fontsize=7, ha="left", va="baseline")

    values = []                                                    # everything plotted -> annotation_ABCD_values.tsv
    ax_y, ax_h, head_y = 3.65, 1.8, 5.75                           # row-1 axes bottom and height; row-1 headings

    def stacked(ax, shares, cats, colors, panel, title):
        x = np.arange(len(ANN_SETS))
        bottom = np.zeros(len(ANN_SETS))
        for k, (cat, col) in enumerate(zip(cats, colors)):
            v = np.array([shares[s][k] for s in ANN_SETS])
            ax.bar(x, v, 0.64, bottom=bottom, color=col, edgecolor="white", linewidth=0.5, label=cat)
            for i in range(len(x)):
                if v[i] >= 5:
                    ax.text(x[i], bottom[i] + v[i] / 2, f"{v[i]:.0f}%", ha="center", va="center", fontsize=5.8,
                            color="white" if dark(col) else "black")
            values.extend({"panel": panel, "measure": f"{title}: {cat}", "set": s, "percent": round(vi, 3)}
                          for s, vi in zip(ANN_SETS, v))
            bottom += v
        ax.set_xticks(x, [ANN_SHORT[s] for s in ANN_SETS])
        ax.tick_params(axis="x", labelsize=6.2)
        ax.set_xlim(-0.55, len(x) - 0.45)
        ax.set_ylim(0, 100)
        ax.set_yticks(range(0, 101, 25))
        ax.set_ylabel("Proteins (%)")
        ax.legend(loc="lower left", bbox_to_anchor=(0, 1.02), ncol=len(cats), frameon=False, fontsize=6.0,
                  handlelength=0.8, handleheight=0.8, handletextpad=0.3, columnspacing=0.6, borderaxespad=0)
        style(ax, lw=0.8)

    # A, UniProt annotation score (1-5)
    score = {s: [(g.annotation_score_num == k).mean() * 100 for k in range(1, 6)] for s, g in d.items()}
    stacked(box(0.47, ax_y, 1.15, ax_h), score, ["1", "2", "3", "4", "5"], SCORE_COLOR, "A",
            "UniProt annotation score")
    heading("A", 0.02, head_y, "UniProt annotation score")

    # B, distinct Pfam families per protein
    v = {s: g.pfam_unique_families for s, g in d.items()}
    npf = {s: [(x == 0).mean() * 100, (x == 1).mean() * 100, (x == 2).mean() * 100, (x >= 3).mean() * 100]
           for s, x in v.items()}
    stacked(box(2.07, ax_y, 1.15, ax_h), npf, NPF_CATS, NPF_COLOR, "B", "Pfam families per protein")
    heading("B", 1.66, head_y, "Pfam families per protein")

    # C, annotation features (% of proteins), with the set legend under its heading
    c = box(4.35, ax_y, 1.57, ax_h)
    off = {"sel3800": 0.28, "sel690": 0.0, "pdb3825_full": -0.28}
    for gi, (col, lab) in enumerate(INDICATORS):
        for s in ANN_SETS:
            v = d[s][col].astype(float).mean() * 100
            c.barh(-gi + off[s], v, 0.26, color=COLOR[s], edgecolor="none")
            c.text(v + 1.5, -gi + off[s], f"{v:.1f}", va="center", ha="left", fontsize=5.5)
            values.append({"panel": "C", "measure": lab, "set": s, "percent": round(v, 3)})
    c.set_yticks([-i for i in range(len(INDICATORS))], [lab for _, lab in INDICATORS])
    c.set_ylim(-len(INDICATORS) + 0.5, 0.5)
    c.axhline(-3.5, color="#bbbbbb", lw=0.6, ls=(0, (3, 2)))          # presence of annotation | lack of it
    c.set_xlim(0, 118)
    c.set_xticks(range(0, 101, 25))
    c.set_xlabel("Proteins (%)")
    style(c, lw=0.8)
    c.tick_params(axis="y", length=0)
    heading("C", 3.3, head_y, "Annotation features")
    fig.legend([plt.Rectangle((0, 0), 1, 1, color=COLOR[s]) for s in ANN_SETS],
               [f"{ANN_LABEL[s]} ({n[s]:,} proteins)" for s in ANN_SETS], loc="upper left",
               bbox_to_anchor=(3.3 / W, (head_y - 0.06) / H), ncol=2, frameon=False, fontsize=5.5, handlelength=0.9,
               handleheight=0.7, columnspacing=0.8, borderaxespad=0, handletextpad=0.35, labelspacing=0.25)

    # D, the most frequent best Pfam families of each set: a ranked list per set, names on light bars scaled to the
    # most frequent family of the set, protein counts on the right
    top = {s: top_families(g.pfam_best_family.dropna().value_counts()) for s, g in d.items()}
    dx = box(0, 0, W, 3.25)                                        # data coordinates = inches
    dx.set_xlim(0, W)
    dx.set_ylim(0, 3.25)
    dx.set_axis_off()
    heading("D", 0.02, 3.08, "Most frequent Pfam families of each set (family of the best-scoring domain of each "
            "protein)")
    rh, y0, cw, bw = 0.18, 2.78, 1.9, 1.6                           # row height, top of the rows, column, bar width
    for j, s in enumerate(ANN_SETS):
        x0 = 0.05 + j * 2.0
        dx.add_patch(plt.Rectangle((x0, 2.835), 0.09, 0.09, color=COLOR[s], lw=0))
        dx.text(x0 + 0.14, 2.88, ANN_LABEL[s], fontsize=6.5, fontweight="bold", ha="left", va="center")
        dx.text(x0 + cw, 2.88, "Proteins", fontsize=5.8, color="#444444", ha="right", va="center")
        kmax = top[s].max()
        for i, (fam, k) in enumerate(top[s].items()):
            yc = y0 - (i + 0.5) * rh
            dx.add_patch(plt.Rectangle((x0, yc - 0.39 * rh), bw * k / kmax, 0.78 * rh, color=COLOR[s], alpha=0.32,
                                       lw=0))
            dx.text(x0 + 0.03, yc, wrap(pf_name[fam], 40), fontsize=5.6, ha="left", va="center", linespacing=0.9)
            dx.text(x0 + cw, yc, f"{k}", fontsize=5.6, ha="right", va="center")
            values.append({"panel": "D", "measure": f"best Pfam family {fam} ({pf_name[fam]})", "set": s,
                           "percent": round(k / n[s] * 100, 3), "proteins": int(k), "rank": i + 1})

    fig.savefig(OUT / "annotation_ABCD.png", dpi=300)
    fig.savefig(OUT / "annotation_ABCD.pdf")
    plt.close(fig)
    pd.DataFrame(values).to_csv(OUT / "annotation_ABCD_values.tsv", sep="\t", index=False)

    # numbers for the caption
    for s, g in d.items():
        withpf = g.pfam_best_family.notna()
        print(f"{s}: n={n[s]:,}  score 1 {score[s][0]:.1f}%  median score {g.annotation_score_num.median():.0f}  "
              f">= 2 families {npf[s][2] + npf[s][3]:.1f}%  with Pfam {withpf.sum():,}  distinct best families "
              f"{g.pfam_best_family.nunique():,}  shown {len(top[s])} (>= {top[s].min()} proteins), "
              f"{top[s].sum() / withpf.sum() * 100:.1f}% of those with Pfam", flush=True)
        print("   ", list(top[s].index), flush=True)


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
        {"ages": ages, "crops": crops, "msa": msa, "disorder": disorder, "v4v6": v4v6, "umap": umap,
         "annotation": annotation}[job]()
