"""Checks before redrawing Figs. S11-S13 on the selection's definitions (candidate pool = 179,596 unique plant
sequences without a PDB structure; library = 186,823 unique sequences of the September 2026 PDB, or without pdb3825):
1. counterfactual round-0 minima and pool sums agree with sel3800 / sel690 round 0 (Dataset S17 FigS12/raw)
2. S11 numbers: nearest library sequence in pdb3825 (black dots), medians with and without pdb3825
3. S13: map the 175,537 evaluated rows of evaluation_distances.parquet to pool positions (rows are in pool order;
   matched on the exact float32 values of base_3800 = d_l0 and real_3800 = d_full), then the pool-based average
   homolog distance versus the d_homolog column (all 182,286 plant sequences)
Usage: python s11_s13_harmonize_check.py <pull dir>"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

P = Path(sys.argv[1])
S17 = Path("/path/to/project/SI_2026/datasets/Dataset_S17_source_data_SI_figures")
z = np.load(P / "cf_mins_round0.npz")
d_l0, d_690, d_rest = z["d_l0"], z["d_690"], z["d_rest"]
d_full = np.minimum(np.minimum(d_l0, d_690), d_rest)
s0 = np.load(P / "cf_pool_sums_round0.npy")
n = len(s0)
avg = (s0 / n).astype(np.float32)
print("pool", n, "arrays", len(d_l0), len(d_690), len(d_rest))

rt = pd.read_csv(S17 / "FigS12/selection_rounds_0_and_37.tsv.gz", sep="\t")
print("FigS12 table columns", list(rt.columns), rt["round"].value_counts().to_dict())
r0 = rt[rt["round"] == 0].sort_values("pool_position")
assert (r0.pool_position.to_numpy() == np.arange(len(r0))).all()
r0_min = r0.d_nearest_library.to_numpy(np.float32)
r0_avg = r0.d_average_homolog.to_numpy(np.float32)
p3800 = np.flatnonzero(r0.selected_this_round.astype(bool).to_numpy())
print("sel3800 round 0 vs cf: min max|diff|", float(np.abs(r0_min - d_full).max()),
      "avg max rel diff", float((np.abs(r0_avg - avg) / r0_avg).max()))
m690, a690 = np.load(P / "sel690_round_000_min_dist.npy"), np.load(P / "sel690_round_000_avg_dist.npy")
print("sel690 round 0 vs sel3800 round 0: min identical", bool(np.array_equal(m690, r0_min)),
      "avg identical", bool(np.array_equal(a690, r0_avg)))
p690 = np.load(P / "sel690_round_000_picked.npy")
print("sel690 round-0 picks within sel3800 round-0 picks:", len(set(p690) & set(p3800)), "of", len(p690))
for r in ("000", "068"):
    mn, av, pk = (np.load(P / f"sel690_round_{r}_{k}.npy") for k in ("min_dist", "avg_dist", "picked"))
    print(f"sel690 round {r}: pool {len(mn):,} picked {len(pk)} median min {np.median(mn):.3f} avg {np.median(av):.3f}; "
          f"picked median min {np.median(mn[pk]):.3f} avg {np.median(av[pk]):.3f}")
print("limits", json.loads((P / "sel690_limits.json").read_text()), json.loads((S17 / "FigS12/limits.json").read_text()))

# S11 on the selection's definitions
black = np.minimum(d_690, d_rest) < d_l0
ties = np.minimum(d_690, d_rest) == d_l0
print(f"S11: nearest library sequence in pdb3825 {int(black.sum()):,} ({black.mean():.1%}); ties {int(ties.sum())}; "
      f"median nearest complete {np.median(d_full):.3f} without pdb3825 {np.median(d_l0):.3f}; black "
      f"{np.median(d_full[black]):.3f} -> {np.median(d_l0[black]):.3f} (median increase "
      f"{np.median((d_l0 - d_full)[black]):.3f}); avg homolog median all {np.median(avg):.3f} black "
      f"{np.median(avg[black]):.3f} others {np.median(avg[~black]):.3f}")
outside = lambda x, y: int(((x > 6) | (y < 2) | (y > 8)).sum())
print("S11 outside 0-6 x 2-8:", outside(d_full, avg), outside(d_l0, avg), "x < 0.05:", int((d_full < 0.05).sum()),
      int((d_l0 < 0.05).sum()))

# S13 row -> pool position
ev = pd.read_parquet(S17 / "FigS13/evaluation_distances.parquet")
b, r = ev.base_3800.to_numpy(np.float32), ev.real_3800.to_numpy(np.float32)
pos = np.empty(len(ev), dtype=np.int64)
j = 0
for i in range(len(ev)):
    while not (d_l0[j] == b[i] and d_full[j] == r[i]):
        j += 1
    pos[i] = j
    j += 1
print("S13: evaluated", len(ev), "mapped, positions strictly increasing", bool((np.diff(pos) > 0).all()),
      "skipped", n - len(ev))
assert np.array_equal(ev.base_690.to_numpy(np.float32), np.minimum(d_l0, d_rest)[pos]), "base_690 mismatch"
y_pool = avg[pos]
y_all = ev.d_homolog.to_numpy(np.float32)
print(f"S13 y: pool-based median {np.median(y_pool):.3f} vs all-plant {np.median(y_all):.3f}; max |diff| "
      f"{np.abs(y_pool - y_all).max():.3f}; median |diff| {np.median(np.abs(y_pool - y_all)):.4f}; spearman "
      f"{pd.Series(y_pool).corr(pd.Series(y_all), method='spearman'):.5f}")
np.save(P / "s13_pool_positions.npy", pos)
