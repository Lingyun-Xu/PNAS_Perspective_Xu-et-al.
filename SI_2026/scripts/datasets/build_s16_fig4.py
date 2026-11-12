#!/usr/bin/env python3
"""Add the AFDB v6 version of main-text Fig. 4 to Dataset S16 (folder Fig4_AFDB_v6; the first-version folder Fig3 is
kept) and describe it in the Dataset S16 README.
python build_s16_fig4.py <updated_maintext_figure> <SI_2026/datasets/Dataset_S16_source_data_main_text_figures>"""
import shutil
import sys
from pathlib import Path

import pandas as pd

R, S16 = Path(sys.argv[1]), Path(sys.argv[2])
RES = R / "data/fig4_v6_results"
OUT = S16 / "Fig4_AFDB_v6"
OUT.mkdir(exist_ok=True)
shutil.copy2(R / "data/Fig4_v6_source_data.tsv.gz", OUT / "Fig4_v6_source_data.tsv.gz")
shutil.copy2(R / "statistics/fig4_v6_summary.tsv", OUT / "Fig4_v6_summary.tsv")
names = {"thermostability": "Thermostability", "thermostability_2025": "Thermostability (Datasets 2 and 3 only)",
         "mavedb": "MaveDB", "mdcath": "mdCATH"}
sel = pd.concat([pd.read_csv(RES / f"{k}_selected.tsv", sep="\t").assign(row=v) for k, v in names.items()])
sel = sel[["row"] + [c for c in sel.columns if c != "row"]]
sel.to_csv(OUT / "Fig4_v6_selected_candidates.tsv.gz", sep="\t", index=False, compression="gzip", float_format="%.6g")
print({k: int((sel.row == v).sum()) for k, v in names.items()})

readme = S16 / "README.txt"
text = readme.read_text(encoding="utf-8")
old_idx = ("  Fig. 4            Fig3/    Fig3f-j_thermostability.tsv.gz (a-e), Fig3k-o_MaveDB.tsv.gz (f-j), Fig3p-t_mdCATH.tsv.gz\n"
           "                    (k-o): the thermostability, MaveDB and mdCATH rows of the first version of Fig. 3\n")
new_idx = ("  Fig. 4 (AFDB v6)  Fig4_AFDB_v6/  Fig4_v6_source_data.tsv.gz (a-o), Fig4_v6_selected_candidates.tsv.gz,\n"
           "                    Fig4_v6_summary.tsv (November 2026; candidate pool of main-text Fig. 3)\n"
           "  Fig. 4 (first version, AFDB v4 plant sequences)  Fig3/  Fig3f-j_thermostability.tsv.gz (a-e),\n"
           "                    Fig3k-o_MaveDB.tsv.gz (f-j), Fig3p-t_mdCATH.tsv.gz (k-o): the thermostability, MaveDB and\n"
           "                    mdCATH rows of the first version of Fig. 3\n")
assert text.count(old_idx) == 1, "index line not found"
text = text.replace(old_idx, new_idx)
section = """
Fig4_AFDB_v6/   main-text Fig. 4 recomputed on the AFDB v6 candidate pool of main-text Fig. 3 (179,596 unique sequences of
                the AFDB v6 Viridiplantae proteins without a SIFTS-mapped PDB chain; SI Section S3). Selection as in the
                first version: score = nearest template distance / average homolog distance (no pLDDT term); per round
                the top 500 (thermostability, mdCATH) or 50 (MaveDB) candidates move into the template pool; panels show
                rounds 0, 1 and 20 (40 for MaveDB). Distances: Euclidean between ESM-2 (650M) mean embeddings (layer
                33, sequences truncated at 1,022 residues). Template pools:
                  thermostability  1,778,010 unique sequences: FireProtDB wild types (242), every protein sequence of
                                   the 1,841,285 rows of mega-scale Dataset 1 (DNA translated, SAGG linker padding
                                   removed; 1,734,297) with the 3,138 sequences found only in Datasets 2 and 3, and the
                                   Meltome Atlas proteins (40,381)
                  thermostability (Datasets 2 and 3 only)  the pool of the first version: FireProtDB, mega-scale
                                   Datasets 2 and 3 (740,697) and Meltome Atlas, 781,273 unique sequences (comparison)
                  MaveDB           1,170 unique target sequences (export of 14 November 2024)
                  mdCATH           5,398 domains (ATOM-record sequences of the CATH v4.2.0 S20 non-redundant set)
  Fig4_v6_source_data.tsv.gz  row, round (0, 1, last), candidate_slot (position in candidates of the run), seq_id,
                              nearest_template_distance (x), average_homolog_distance (y; mean distance to the other
                              candidates still in the pool, own zero distance excluded), for every candidate still in
                              the pool in that round
  Fig4_v6_selected_candidates.tsv.gz  the selected candidates: row, round, rank_in_round, seq_id,
                              representative_accession, organism_scientific, plddt (AFDB v6 global), score, d_min,
                              d_avg (at selection)
  Fig4_v6_summary.tsv         per panel: candidates, median and mean nearest template distance, shares within 1.9 and
                              1.0, median homolog distance, points outside the plotted range, share brought closer
                              than in round 0
"""
if "Fig4_AFDB_v6/   main-text Fig. 4" not in text:
    text = text.rstrip("\n") + "\n" + section
readme.write_text(text, encoding="utf-8", newline="")
print("Dataset S16 updated:", sorted(p.name for p in OUT.iterdir()))
