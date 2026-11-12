#!/usr/bin/env python3
"""Fig. S8 of the SI for the rebuttal with human added (author request, 9 Nov 2026): published gene ages of the eight
GenEra species (Barrera-Redondo et al. 2023), as in s10() of si_pnas_server.py, plus the human map of Lotharukpong (2023;
phylomapr, Homo_sapiens.PhyloMap). The human map gives numbered phylostrata only (1-28). Rank k is given the name of
the k-th GenEra phylostratum of the human NCBI lineage (GenEra omits Craniata, Teleostomi and the genus Homo, leaving
28 strata; the lineage is the same in every NCBI Taxonomy release from 2019 to 2026), and the names are binned with
the same code as the other species (plot_published_suite.assign_names: map_gene_ages.classify, then refine).
Writes (new files only): /path/to/server_data/si_pnas/out/figS8_published_ages_9species.png and figS8_human_binned.tsv."""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, "/path/to/server_data/conservation/scripts")
import plot_published_suite as ps

OUT = Path("/path/to/server_data/si_pnas/out")
plt.rcParams.update({"font.family": "DejaVu Sans", "axes.titleweight": "normal"})
HUMAN_STRATA = ["cellular organisms", "Eukaryota", "Opisthokonta", "Metazoa", "Eumetazoa", "Bilateria", "Deuterostomia",
                "Chordata", "Vertebrata", "Gnathostomata", "Euteleostomi", "Sarcopterygii", "Dipnotetrapodomorpha",
                "Tetrapoda", "Amniota", "Mammalia", "Theria", "Eutheria", "Boreoeutheria", "Euarchontoglires",
                "Primates", "Haplorrhini", "Simiiformes", "Catarrhini", "Hominoidea", "Hominidae", "Homininae",
                "Homo sapiens"]
ORDER = ["ARATH", "ORYSJ", "SOYBN", "MOUSE", "DROME", "CAEEL", "HUMAN", "YEAST", "SCHPO"]
LABEL = {row[0]: row[4] for row in ps.BARRERA}
LABEL["HUMAN"] = "Human*"


def frame(ax, lw, tick_len, labelsize):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(lw)
    ax.tick_params(width=lw, length=tick_len, labelsize=labelsize)


def letter(ax, text, x, y, size):
    ax.text(x, y, text, transform=ax.transAxes, fontsize=size, fontweight="bold", va="bottom", ha="left")


parent, rank, names = ps.mga.load_taxdump(ps.TAXDUMP)
meta = ps.load_meta()
sgd = ps.load_sgd()
ages = ps.assign_names(ps.barrera_genes(sgd), parent, rank, names)
parts = []
for proteome, *_rest in ps.BARRERA:
    hit = ps.join_meta(ages[ages.proteome.eq(proteome)], meta, proteome)
    parts.append(hit.sort_values("rank").drop_duplicates("db_accession"))

human = ps.load_human()
assert human["rank"].between(1, 28).all(), human["rank"].describe()
human["phylostratum"] = human["rank"].map(lambda k: HUMAN_STRATA[k - 1])
human = ps.assign_names(human, parent, rank, names)
print("human rank -> name -> bin:")
print(human.drop_duplicates("rank").sort_values("rank")[["rank", "phylostratum", "age_bin_lumped", "age_bin"]]
      .to_string(index=False), flush=True)
hh = ps.join_meta(human, meta, "HUMAN").sort_values("rank").drop_duplicates("db_accession")
parts.append(hh)
df = pd.concat(parts, ignore_index=True)
plot = df[df.age_bin.isin(ps.BINS)].copy()
order = [p for p in ORDER if p in set(plot.proteome)]

fig, axes = plt.subplots(2, 1, figsize=(9.0, 9.6))
ax = axes[0]
bottom = np.zeros(len(order))
for b in ps.BINS:
    heights = [(plot[plot.proteome.eq(p)].age_bin.eq(b)).mean() for p in order]
    ax.bar(np.arange(len(order)), heights, bottom=bottom, color=ps.BIN_COLOR[b], width=0.78, label=ps.BIN_LABEL[b])
    bottom += np.array(heights)
ax.set_xticks(np.arange(len(order)))
ax.set_xticklabels([f"{LABEL[p]}\nn={int(df.proteome.eq(p).sum()):,}" for p in order], fontsize=10, rotation=45,
                   ha="right")
ax.set_ylim(0, 1)
ax.set_ylabel("Fraction of matched proteins", fontsize=11)
frame(ax, 1.2, 3, 10)
ax.legend(frameon=False, fontsize=10, loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=5, title="Age class",
          title_fontsize=10, columnspacing=1.2, handlelength=1.4)
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
handles = [plt.Rectangle((0, 0), 1, 1, color=ps.SPECIES_COLOR[p], label=LABEL[p]) for p in order]
ax.legend(handles=handles, frameon=False, fontsize=9.5, loc="lower center", bbox_to_anchor=(0.5, 1.01),
          ncol=len(handles), title="Species", title_fontsize=10, columnspacing=0.9, handlelength=1.2,
          handletextpad=0.4)
letter(axes[0], "A", -0.08, 1.16, 16)
letter(axes[1], "B", -0.08, 1.16, 16)
fig.tight_layout(h_pad=3.0)
fig.savefig(OUT / "figS8_published_ages_9species.png", dpi=300, bbox_inches="tight")
plt.close(fig)

for p in order:
    s = plot[plot.proteome.eq(p)]
    shares = (s.age_bin.value_counts(normalize=True) * 100).round(1).reindex(ps.BINS).to_dict()
    med = s.groupby("age_bin").plddt_global.median().round(1).reindex(ps.BINS).to_dict()
    n_bin = s.age_bin.value_counts().reindex(ps.BINS).fillna(0).astype(int).to_dict()
    print(f"{p}: matched {int(df.proteome.eq(p).sum()):,}; binned {len(s):,}; shares % {shares}; n {n_bin}; "
          f"median pLDDT {med}", flush=True)
cols = ["proteome", "species_label", "source", "gene", "acc", "db_accession", "query_taxid", "phylostratum", "rank",
        "taxonomic_representativeness", "age_bin_lumped", "age_bin", "age_taxid", "plddt_global", "neff",
        "frac_disordered"]
hh[[c for c in cols if c in hh.columns]].to_csv(OUT / "figS8_human_binned.tsv", sep="\t", index=False)
print("human rows written:", len(hh), "| human proteins in the v6 proteome (meta):", int(meta.proteome.eq("HUMAN").sum()))
