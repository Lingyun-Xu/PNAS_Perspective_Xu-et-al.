"""Source data for Figs. S11-S13 on the selection's definitions (author request, 1 Nov 2026: state who is measured
against whom; same style as Fig. S12 for sel690).
Candidate pool: the 179,596 unique sequences of the 182,501 AFDB v6 plant accessions without a SIFTS-mapped PDB chain,
in the pool order of select_plants.py (pool_position). Distances: raw Euclidean minima / means between ESM-2 mean
embeddings, as in the selection.
  FigS11/candidate_pool_round0.tsv.gz   d_nearest_complete (to the 186,823 sequences of the September 2026 PDB),
      d_nearest_without_pdb3825 (to the 182,998 without pdb3825), d_average_homolog (mean distance to all 179,596
      candidates, including itself), nearest_is_pdb3825 (the closest library sequence is a pdb3825 sequence)
      [counterfactual/mins_round0.npz and pool_sums_round0.npy of counterfactual_selection.py]
  FigS12/sel690_rounds_0_and_68.tsv.gz  as sel3800_rounds_0_and_37.tsv.gz, for sel690 (selection_690/round_000, 068)
  FigS13/evaluation_distances.parquet   + pool_position, d_average_homolog (pool-based, replaces d_homolog, which was
      the mean distance to all 182,286 plant sequences)
Usage: python s11_s13_source_data.py <pull dir> <stage dir>"""
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

P, OUT = Path(sys.argv[1]), Path(sys.argv[2])
S17 = Path("/path/to/project/SI_2026/datasets/Dataset_S17_source_data_SI_figures")
for d in ("FigS11", "FigS12", "FigS13"):
    (OUT / d).mkdir(parents=True, exist_ok=True)

z = np.load(P / "cf_mins_round0.npz")
d_l0, d_690, d_rest = z["d_l0"], z["d_690"], z["d_rest"]
d_new = np.minimum(d_690, d_rest)
d_full = np.minimum(d_l0, d_new)
s0 = np.load(P / "cf_pool_sums_round0.npy")
avg = (s0 / len(s0)).astype(np.float32)
s11 = pd.DataFrame({"pool_position": np.arange(len(s0)), "d_nearest_complete": d_full,
                    "d_nearest_without_pdb3825": d_l0, "d_average_homolog": avg, "nearest_is_pdb3825": d_new < d_l0})
s11.to_csv(OUT / "FigS11/candidate_pool_round0.tsv.gz", sep="\t", index=False, float_format="%.6g")
print("FigS11", len(s11), "nearest_is_pdb3825", int(s11.nearest_is_pdb3825.sum()))

rows = []
for r in ("000", "068"):
    mn, av, pk = (np.load(P / f"sel690_round_{r}_{k}.npy") for k in ("min_dist", "avg_dist", "picked"))
    sel = np.zeros(len(mn), dtype=bool)
    sel[pk] = True
    rows.append(pd.DataFrame({"round": int(r), "pool_position": np.arange(len(mn)), "d_nearest_library": mn,
                              "d_average_homolog": av, "selected_this_round": sel}))
s12 = pd.concat(rows, ignore_index=True)
s12.to_csv(OUT / "FigS12/sel690_rounds_0_and_68.tsv.gz", sep="\t", index=False, float_format="%.6g")
shutil.copy2(P / "sel690_selected_accessions.tsv", OUT / "FigS12/selected_accessions_sel690.tsv")
print("FigS12 sel690", s12.groupby("round").agg(n=("pool_position", "size"), picked=("selected_this_round", "sum")).to_dict())

ev = pd.read_parquet(S17 / "FigS13/evaluation_distances.parquet")
pos = np.load(P / "s13_pool_positions.npy")
assert len(pos) == len(ev) and np.array_equal(ev.base_3800.to_numpy(np.float32), d_l0[pos])
assert np.array_equal(ev.real_3800.to_numpy(np.float32), d_full[pos])
ev.insert(0, "pool_position", pos)
ev["d_average_homolog"] = avg[pos]
ev = ev.drop(columns="d_homolog")
ev.to_parquet(OUT / "FigS13/evaluation_distances.parquet", index=False)
print("FigS13", len(ev), list(ev.columns))
