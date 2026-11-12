#!/usr/bin/env python3
"""Main-text Fig. 4 redrawn on the AFDB v6 plant proteins of main-text Fig. 3 (computed on the cluster with
fig4_v6_selection.py; SI Section S3). Panels in the style of fig3_s13_panels.py (free-energy maps, 200 x 200 bins over
x 0-6 and y 2-7, levels 0-6; differences 150 x 150 bins, levels -6 to 6, seismic; dataset 1 minus dataset 2, so red
marks a density loss and blue a gain).
Rows (template pool): A-E thermostability, F-J MaveDB, K-O mdCATH. Columns: (i) template pool alone (round 0), (ii)
after the first round of AI selections (round 1), (iii) after all rounds (round 20; 40 for MaveDB), (ii) minus (i),
(iii) minus (i). x, nearest template distance (to the nearest sequence of that template pool, the proteins selected
so far included); y, average homolog distance (mean distance to the other candidates still in the pool).
Outputs (Fig4/): Fig4A.png ... Fig4O.png without labels (frames and ticks only, as the Fig. 3 panels) and
Fig4A_labeled.png ... with axis labels, tick labels, title and colour bar; colour bars; Fig4_overview_labeled.png;
statistics/fig4_v6_summary.tsv; data/Fig4_v6_source_data.tsv.gz.
Usage: python fig4_v6_panels.py <results folder from the cluster> <output root (updated_maintext_figure)>"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fig3_s13_panels import (VMAX, colorbar, custom_cmap, density, difference, free_energy_plot,  # noqa: E402
                             free_energy_plot_difference)

ROWS = [("thermostability", "Thermostability", "Nearest Thermostability Template Distance", 20, 500, "ABCDE"),
        ("mavedb", "MaveDB", "Nearest Function Template Distance", 40, 50, "FGHIJ"),
        ("mdcath", "mdCATH", "Nearest Dynamic Template Distance", 20, 500, "KLMNO")]
YLABEL = "Homologs Distance"
UNIQUE = {"thermostability": 1778010, "thermostability_2025": 781273, "mavedb": 1170, "mdcath": 5398}   # sequences
# (the thermostability pools hold 47 repeated sequences, shared between FireProtDB, Meltome Atlas and mega-scale)


def labeled(kind, a, b, xlabel, title, path):
    """one panel with tick labels, axis labels, title and its colour bar"""
    fig, ax = plt.subplots(figsize=(7.2, 5.4), dpi=300)
    if kind == "density":
        cf = free_energy_plot(a[0], a[1], ax=ax, cmap=custom_cmap, vmax=VMAX)
        ticks, clabel = list(range(0, VMAX + 1)), "-log(probability)"
    else:
        cf = free_energy_plot_difference(a[0], a[1], b[0], b[1], ax=ax, cmap="seismic", vmax=VMAX)
        ticks, clabel = list(range(-VMAX, VMAX + 1, 2)), "Difference of -log(probability)"
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("bottom", "left"):
        ax.spines[side].set_linewidth(1.6)
    ax.set_xlim(0, 6)
    ax.set_ylim(2, 7)
    ax.set_xticks(range(0, 7))
    ax.set_yticks(range(2, 8))
    ax.tick_params(width=1.6, length=5, labelsize=12)
    ax.set_xlabel(xlabel, fontsize=13)
    ax.set_ylabel(YLABEL, fontsize=13)
    ax.set_title(title, fontsize=12)
    cb = fig.colorbar(cf, ax=ax, ticks=ticks, fraction=0.05, pad=0.03)
    cb.set_label(clabel, fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def stats(row, panel, what, x, y, base=None):
    r = {"row": row, "panel": panel, "what": what, "candidates": int(len(x)),
         "median nearest template distance": float(np.median(x)), "mean nearest template distance": float(np.mean(x)),
         "within 1.9 (%)": float((x <= 1.9).mean() * 100), "within 1.0 (%)": float((x <= 1.0).mean() * 100),
         "median homolog distance": float(np.median(y)),
         "outside x 0-6 or y 2-7": int(((x < 0) | (x > 6) | (y < 2) | (y > 7)).sum())}
    if base is not None:
        r["median nearest distance, same candidates in round 0"] = float(np.median(base))
        r["closer than in round 0 (%)"] = float(((base - x) > 1e-4).mean() * 100)
    return r


def main():
    RES, ROOT = Path(sys.argv[1]), Path(sys.argv[2])
    OUT, STATS, DATA = ROOT / "Fig4", ROOT / "statistics", ROOT / "data"
    for d in (OUT, STATS, DATA):
        d.mkdir(parents=True, exist_ok=True)
    cand = pd.read_csv(RES / "candidates.tsv.gz", sep="\t")
    rows, src = [], []
    fig_o, axes = plt.subplots(3, 5, figsize=(26, 14), dpi=150)
    todo = list(enumerate(ROWS))
    if (RES / "thermostability_2025_states.npz").exists():    # comparison: the thermostability pool of the first version
        todo.append((0, ("thermostability_2025", "Thermostability (Datasets 2 and 3 only)", ROWS[0][2], 20, 500,
                         "ABCDE")))
    for i, (tag, name, xlabel, R, per, letters) in todo:
        if not (RES / f"{tag}_states.npz").exists():
            print("no results yet for", tag)
            continue
        out_dir = OUT / "thermostability_pool_first_version" if tag == "thermostability_2025" else OUT
        out_dir.mkdir(parents=True, exist_ok=True)
        z = np.load(RES / f"{tag}_states.npz")
        summ = json.loads((RES / f"summary_{tag}.json").read_text())
        st = {r: (z[f"slots_r{r}"], z[f"d_min_r{r}"].astype(float), z[f"d_avg_r{r}"].astype(float)) for r in (0, 1, R)}
        lib = UNIQUE[tag]
        assert summ["template_pool"] - lib in (0, 47), (tag, summ["template_pool"])
        base_of = dict(zip(st[0][0], st[0][1]))
        titles = [f"{name} template pool ({lib:,} sequences)",
                  f"+ first round of AI selections ({per:,})",
                  f"+ all {R} rounds of AI selections ({R * per:,})",
                  f"{letters[1]} minus {letters[0]}", f"{letters[2]} minus {letters[0]}"]
        for k, r in enumerate((0, 1, R)):
            slots, x, y = st[r]
            density(x, y, out_dir / f"Fig4{letters[k]}.png")
            labeled("density", (x, y), None, xlabel, titles[k], out_dir / f"Fig4{letters[k]}_labeled.png")
            if out_dir == OUT:
                free_energy_plot(x, y, ax=axes[i][k], cmap=custom_cmap, vmax=VMAX)
            rows.append(stats(name, letters[k], f"round {r}", x, y,
                              base=None if r == 0 else np.array([base_of[s] for s in slots])))
            src.append(pd.DataFrame({"row": name, "round": r, "candidate_slot": slots,
                                     "seq_id": cand.seq_id.to_numpy()[slots],
                                     "nearest_template_distance": np.round(x, 5), "average_homolog_distance": np.round(y, 5)}))
        for k, r in ((3, 1), (4, R)):
            difference(st[r][1], st[r][2], st[0][1], st[0][2], out_dir / f"Fig4{letters[k]}.png")
            labeled("diff", (st[r][1], st[r][2]), (st[0][1], st[0][2]), xlabel, titles[k],
                    out_dir / f"Fig4{letters[k]}_labeled.png")
            if out_dir == OUT:
                free_energy_plot_difference(st[r][1], st[r][2], st[0][1], st[0][2], ax=axes[i][k], cmap="seismic",
                                            vmax=VMAX)
        for k in range(5 if out_dir == OUT else 0):
            ax = axes[i][k]
            ax.set_xlim(0, 6)
            ax.set_ylim(2, 7)
            ax.set_title(f"{letters[k]}  {titles[k]}", loc="left", fontsize=11)
            ax.set_xlabel(xlabel, fontsize=10)
            ax.set_ylabel(YLABEL, fontsize=10)
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)
    colorbar("free", OUT / "Fig4_colorbar_free_energy.png")
    colorbar("diff", OUT / "Fig4_colorbar_difference.png")
    fig_o.suptitle("Main-text Fig. 4 on the AFDB v6 candidate pool (179,596 plant sequences without a structure)",
                   fontsize=14)
    fig_o.tight_layout()
    fig_o.savefig(OUT / "Fig4_overview_labeled.png", dpi=150, bbox_inches="tight")
    plt.close(fig_o)
    S = pd.DataFrame(rows)
    S.to_csv(STATS / "fig4_v6_summary.tsv", sep="\t", index=False)
    pd.concat(src).to_csv(DATA / "Fig4_v6_source_data.tsv.gz", sep="\t", index=False, compression="gzip")
    pd.set_option("display.width", 250)
    print(S.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
