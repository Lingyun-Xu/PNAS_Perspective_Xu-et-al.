#!/usr/bin/env python3
"""SI Fig. S13 (author request, 6 Nov 2026; one panel since 6 Nov 2026, afternoon): change of the average homolog
distance of the candidates when the proteins of sel3800 or sel690 leave the candidate pool.

Average homolog distance (d_avg) of a candidate: the mean ESM-2 embedding distance to the other candidates of the
candidate pool (179,596 unique AFDB v6 plant sequences without a PDB chain); the candidate itself is not counted. It
depends only on which proteins are candidates: adding sel3800 or sel690 to a template pool removes these candidates
from the candidate pool, whereas adding pdb3825 (PDB sequences, not candidates) changes nothing.

Input (Dataset S17, FigS13): davg_change.tsv.gz, the 175,537 candidates in neither selection nor in the AI arms of
Table S30: d_avg, d_avg_without_sel3800, d_avg_without_sel690, change_sel3800, change_sel690
(apply_davg_self_excluded.py, from davg_self_excluded.py and davg_change.py on the cluster).
Writes <out>/figS13_average_homolog_distance.png and davg_change_stats.tsv.
Usage: python figS13_average_homolog_distance.py <Dataset_S17/FigS13 folder> <out folder>"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SRC, OUT = Path(sys.argv[1]), Path(sys.argv[2])
COLOR = {"pdb3825": "#66a61e", "sel3800": "#d95f02", "sel690": "#e6ab02"}
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 7.5,
                     "ytick.labelsize": 7.5, "legend.fontsize": 7.5, "axes.linewidth": 0.8})

t = pd.read_csv(SRC / "davg_change.tsv.gz", sep="\t")
ch = {tag: t[f"change_{tag}"].to_numpy(float) for tag in ("sel3800", "sel690")}
lo = min(np.percentile(d, 0.05) for d in ch.values())
hi = max(np.percentile(d, 99.95) for d in ch.values())
bins = np.linspace(min(lo, 0), max(hi, 0), 120)

fig, ax = plt.subplots(figsize=(4.6, 3.0))
for tag, ls in (("sel3800", "-"), ("sel690", "--")):
    ax.hist(ch[tag], bins=bins, histtype="step", color=COLOR[tag], lw=1.2, ls=ls,
            label=f"{tag} leaves the candidate pool")
ax.axvline(0, color=COLOR["pdb3825"], lw=1.5, label="pdb3825 added: no change")
ax.set_xlabel("Change in average homolog distance")
ax.set_ylabel("Candidates")
handles, labels = ax.get_legend_handles_labels()
fig.legend(handles, labels, frameon=False, loc="upper center", ncol=2, bbox_to_anchor=(0.55, 1.0),
           columnspacing=1.2, handlelength=1.8)
for side in ("top", "right"):
    ax.spines[side].set_visible(False)
ax.tick_params(width=0.8, length=2.5)
fig.tight_layout(rect=(0, 0, 1, 0.86))
OUT.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT / "figS13_average_homolog_distance.png", dpi=500)
plt.close(fig)

rows = [{"set": "none (full candidate pool)", "candidates": len(t), "median d_avg": t.d_avg.median()}]
for tag in ch:
    rows.append({"set": tag, "candidates": len(t), "median d_avg": t[f"d_avg_without_{tag}"].median(),
                 "median change": np.median(ch[tag]), "min change": ch[tag].min(), "max change": ch[tag].max(),
                 "share increased": float((ch[tag] > 0).mean())})
st = pd.DataFrame(rows)
st.to_csv(OUT / "davg_change_stats.tsv", sep="\t", index=False)
print(st.round(5).to_string(index=False))
