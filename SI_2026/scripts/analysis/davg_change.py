#!/usr/bin/env python3
"""Average homolog distance (d_avg) of the candidate pool with the baseline library and after adding pdb3825, sel3800 or
sel690 to it (new SI Fig. S13; author request, 6 Nov 2026). Runs on the cluster, next to counterfactual_selection.py.

d_avg of a candidate = mean ESM-2 distance to all candidates still in the pool, itself included (select_plants.py).
Adding sequences to the template library changes d_avg only when they leave the pool:
  pdb3825 (PDB sequences)  no candidate leaves the pool, so d_avg is unchanged;
  sel3800, sel690           these pool proteins leave it, and d_avg of the remaining candidates is taken over the
                            smaller pool: d_after(i) = (S0(i) - sum_{j in X} d(i, j)) / (N - |X|), S0 the round-0
                            sums of distances (counterfactual/pool_sums_round0.npy).
The picks selected on the baseline library (cf3800, cf690; the AI arm of the counterfactual comparison) are done too.
Outputs (campaign/davg_change/):
  davg_eval.tsv   the 175,537 evaluated candidates (selected in no arm): d_avg with the baseline library and after
                  each addition
  davg_added.tsv  each added protein: its d_avg against the pool (mean ESM-2 distance to all pool candidates) and the
                  number of pool candidates within ESM-2 distance 1.0 / 1.5 / 2.0 (a pool protein itself not counted)
  summary.json    sizes, medians and checks
"""
import json
import sys
import time

sys.path.insert(0, "/path/to/cluster/redo_PDB/ESM_embedding/template_v6")
sys.path.insert(0, "/path/to/cluster/redo_PDB/ESM_embedding/template_v6/campaign")

import numpy as np
import pandas as pd

import counterfactual_selection as cs
import plant_library_distances as pl
import select_plants as sp

OUT = cs.CAMP / "davg_change"
RADII = (1.0, 1.5, 2.0)


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def sums_and_counts(query, q_norm, lib, lib_norm):
    """per query row: sum of the distances to all library rows, and the number of library rows within each radius"""
    chunk = max(64, int(1e8 // max(1, len(lib))))
    total = np.zeros(len(query), dtype=np.float64)
    counts = {r: np.zeros(len(query), dtype=np.int64) for r in RADII}
    for s in range(0, len(query), chunk):
        e = min(len(query), s + chunk)
        d = cs.distances(query[s:e], lib, q_norm[s:e], lib_norm)
        total[s:e] = d.sum(axis=1, dtype=np.float64)
        for r in RADII:
            counts[r][s:e] = (d < r).sum(axis=1)
    return total, counts


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rep, accs = sp.load_pool()
    emb, index = pl.load_subset(sp.AF_DIRS, set(rep.seq_id))
    meta = sp.align_meta(rep, accs, index)
    slot = {sid: i for i, sid in enumerate(meta["seq_id"])}
    n = len(emb)
    norm = cs.sq_norm(emb)
    s0 = np.load(cs.OUT / "pool_sums_round0.npy")
    assert len(s0) == n, (len(s0), n)
    d0 = s0 / n
    ev_tab = pd.read_parquet(cs.OUT / "evaluation_distances.parquet", columns=["seq_id", "d_homolog"])
    ev = np.array([slot[s] for s in ev_tab.seq_id.astype(str)], dtype=np.int64)
    check = float(np.nanmax(np.abs(d0[ev] - ev_tab.d_homolog.to_numpy(np.float64))))
    log(f"pool {n:,}; evaluated {len(ev):,}; largest |S0/N - d_homolog| {check:.2e}")

    sets = {}
    for tag, folder in (("sel3800", "selection_3800"), ("sel690", "selection_690")):
        sets[tag] = np.array(sorted({s for v in cs.saved_selection(folder, slot).values() for s in v}), dtype=np.int64)
    for tag, name in (("cf3800", "cf_sel3800"), ("cf690", "cf_sel690")):
        t = pd.read_csv(cs.OUT / f"{name}_selected.tsv", sep="\t", dtype={"seq_id": str})
        sets[tag] = np.array(sorted(slot[s] for s in t.seq_id), dtype=np.int64)
    log({k: len(v) for k, v in sets.items()})
    assert not np.isin(ev, np.concatenate(list(sets.values()))).any()

    table = pd.DataFrame({"seq_id": ev_tab.seq_id.astype(str), "d_base": d0[ev]})
    table["d_plus_pdb3825"] = table.d_base          # pdb3825 are PDB sequences: no candidate leaves the pool
    added = []
    q, qn = emb[ev], norm[ev]
    for tag, rows in sets.items():
        removed, _ = sums_and_counts(q, qn, emb[rows], norm[rows])
        table[f"d_plus_{tag}"] = (s0[ev] - removed) / (n - len(rows))
        total, counts = sums_and_counts(emb[rows], norm[rows], emb, norm)
        rec = pd.DataFrame({"set": tag, "id": [meta["seq_id"][r] for r in rows], "d_avg_to_pool": total / n})
        for r in RADII:
            rec[f"n_pool_within_{r}"] = counts[r] - 1          # the protein itself (distance 0) is in the pool
        added.append(rec)
        log(f"{tag}: |X| {len(rows):,}; median change {np.median(table[f'd_plus_{tag}'] - table.d_base):+.5f}")

    rec = pd.read_csv(cs.OUT / "recent_plant_pdb_sequences.tsv", sep="\t", dtype={"embed_id": str})
    pemb, pindex = pl.load_subset(sp.PDB_DIRS, set(rec.embed_id))
    ids = [e for e in rec.embed_id if e in pindex]
    x = np.asarray(pemb[[pindex[e] for e in ids]], dtype=np.float32)
    total, counts = sums_and_counts(x, cs.sq_norm(x), emb, norm)
    r = pd.DataFrame({"set": "pdb3825", "id": ids, "d_avg_to_pool": total / n})
    for rad in RADII:
        r[f"n_pool_within_{rad}"] = counts[rad]
    added.append(r)
    log(f"pdb3825: {len(ids):,} of {len(rec):,} with an embedding")

    added = pd.concat(added, ignore_index=True)
    added.to_csv(OUT / "davg_added.tsv", sep="\t", index=False)
    table.to_csv(OUT / "davg_eval.tsv", sep="\t", index=False)
    summary = {"pool": int(n), "evaluated": int(len(ev)), "check_max_abs_S0N_vs_d_homolog": check,
               "sizes": {k: int(len(v)) for k, v in sets.items()} | {"pdb3825": len(ids)},
               "median_d_avg_pool": float(np.median(d0)), "median_d_avg_evaluated": float(np.median(d0[ev])),
               "median_d_avg_added": added.groupby("set").d_avg_to_pool.median().to_dict(),
               "median_n_pool_within_1.5": added.groupby("set")["n_pool_within_1.5"].median().to_dict()}
    for tag in sets:
        delta = table[f"d_plus_{tag}"] - table.d_base
        summary[f"change_{tag}"] = {"median": float(delta.median()), "min": float(delta.min()), "max": float(delta.max()),
                                    "share_up": float((delta > 0).mean())}
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2), flush=True)
    (OUT / "DAVG_DONE").write_text("done\n")
    log("DAVG_DONE")


if __name__ == "__main__":
    main()
