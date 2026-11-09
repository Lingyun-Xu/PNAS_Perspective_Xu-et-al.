#!/usr/bin/env python3
"""Average homolog distance without the candidate's zero distance to itself (author request, 6 Nov 2026).

The selection code (select_plants.py, counterfactual_selection.py) and the original selection
(legacy_850/analysis_avg_vs_min_L2_scaled_plddt.py) divide the sum of a candidate's distances to the candidate pool by the
number of candidates N, which counts the candidate itself (distance 0). Within a round this multiplies every candidate's
value by the same factor (N - 1)/N, so rankings and selections do not change; the values reported in the SI are now taken
over the other candidates only: d_avg(i) = sum over j != i of d(i, j) / (N - 1).

  A) candidate pool (179,596 unique sequences): d_avg of every candidate in round 0, from counterfactual/pool_sums_round0.npy
     minus the self term, with its accessions (to map the search cohorts of Figs. S15-S18, which are keyed by accession).
  B) the compared sets of Section S7.5 (queries of pdb3825/family_metrics.py) against the 182,286 unique plant sequences of
     the AFDB v6 set: the mean distance to all plant sequences other than identical ones (family_metrics.py included the
     identical sequence, distance 0, in d_homolog_mean).
Distances: Euclidean between ESM-2 (650M) mean embeddings, computed as in counterfactual_selection.distances.
Outputs (campaign/davg_self_excluded/): candidate_pool_davg.tsv.gz, sets_mean_plant.tsv, summary.json, DAVGX_DONE
"""
import json
import sys
import time

sys.path.insert(0, "/path/to/cluster/redo_PDB/ESM_embedding/template_v6")
sys.path.insert(0, "/path/to/cluster/redo_PDB/ESM_embedding/template_v6/campaign")
sys.path.insert(0, "/path/to/cluster/redo_PDB/ESM_embedding/template_v6/campaign/pdb3825")

import numpy as np
import pandas as pd

import counterfactual_selection as cs
import family_metrics as fm
import plant_library_distances as pl
import select_plants as sp

OUT = cs.CAMP / "davg_self_excluded"


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def self_terms(emb, norm, chunk=2048):
    """the diagonal of counterfactual_selection.distances (float32 rounding leaves it slightly above 0)"""
    out = np.empty(len(emb), dtype=np.float64)
    for s in range(0, len(emb), chunk):
        e = min(len(emb), s + chunk)
        out[s:e] = np.diag(cs.distances(emb[s:e], emb[s:e], norm[s:e], norm[s:e]))
    return out


def part_a():
    rep, accs = sp.load_pool()
    emb, index = pl.load_subset(sp.AF_DIRS, set(rep.seq_id))
    meta = sp.align_meta(rep, accs, index)
    n = len(emb)
    norm = cs.sq_norm(emb)
    s0 = np.load(cs.OUT / "pool_sums_round0.npy")
    assert len(s0) == n, (len(s0), n)
    self_d = self_terms(emb, norm)
    other = s0 - self_d
    t = pd.DataFrame({"pool_position": np.arange(n), "seq_id": meta["seq_id"], "rep_accession": meta["rep_acc"],
                      "accessions": [";".join(a) if isinstance(a, (list, tuple)) else str(a) for a in meta["all_acc"]],
                      "d_avg": other / (n - 1), "d_avg_with_self_old": s0 / n})
    t.to_csv(OUT / "candidate_pool_davg.tsv.gz", sep="\t", index=False)
    ratio = (t.d_avg / t.d_avg_with_self_old).to_numpy()
    log(f"part A: candidates {n:,}; self term max {self_d.max():.2e}; ratio new/old {ratio.min():.9f}-{ratio.max():.9f}")
    return t, self_d


def part_b():
    af = pd.read_parquet(fm.ROOT / "af_queries.parquet",
                         columns=["db_accession", "seq_id", "sequence", "tax_group7", "organism_scientific", "in_pdb",
                                  "plddt_global"])
    af["seq_id"] = af.seq_id.astype(str)
    af["db_accession"] = af.db_accession.astype(str)
    plants = af[af.tax_group7.eq("Viridiplantae")].drop_duplicates("seq_id").copy()
    plants["md5"] = plants.sequence.map(fm.md5)
    q = fm.queries(af)
    lib, lib_index = pl.load_subset(sp.AF_DIRS, set(plants.seq_id))
    lib_md5 = np.array([None] * len(lib), dtype=object)
    for sid, h in zip(plants.seq_id, plants.md5):
        if sid in lib_index:
            lib_md5[lib_index[sid]] = h
    lib_norm = cs.sq_norm(lib)
    m = len(lib)
    af_emb, af_index = pl.load_subset(sp.AF_DIRS, set(q.loc[q.emb_kind.eq("af"), "emb_key"]))
    pdb_emb, pdb_index = pl.load_subset(sp.PDB_DIRS, set(q.loc[q.emb_kind.eq("pdb"), "emb_key"]))
    rows = []
    todo = list(zip(q.set, q.query_id, q.emb_kind, q.emb_key, q.sequence.map(fm.md5)))
    for start in range(0, len(todo), 512):
        block = []
        for st, qid, kind, key, h in todo[start:start + 512]:
            if kind == "af" and key in af_index:
                block.append((st, qid, af_emb[af_index[key]], h))
            elif kind == "pdb" and key in pdb_index:
                block.append((st, qid, pdb_emb[pdb_index[key]], h))
        if not block:
            continue
        x = np.stack([b[2] for b in block]).astype(np.float32)
        d = cs.distances(x, lib, cs.sq_norm(x), lib_norm)
        for i, (st, qid, _, h) in enumerate(block):
            same = lib_md5 == h
            other = d[i][~same]
            rows.append((st, qid, float(other.mean()), int(same.sum()), float(d[i].sum() / m)))
    out = pd.DataFrame(rows, columns=["set", "query_id", "d_mean_plant", "n_identical", "d_mean_plant_with_identical"])
    old = pd.read_csv(fm.OUT / "family_metrics.tsv", sep="\t", dtype={"query_id": str},
                      usecols=["set", "query_id", "d_homolog_mean"])
    out = out.merge(old, on=["set", "query_id"], how="left")
    out.to_csv(OUT / "sets_mean_plant.tsv", sep="\t", index=False)
    log(f"part B: {len(out):,} proteins; plant sequences {m:,}")
    return out, m


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    a, self_d = part_a()
    b, m = part_b()
    ok = b.dropna(subset=["d_homolog_mean"])
    summary = {"candidates": int(len(a)), "self_term_max": float(self_d.max()), "self_term_median": float(np.median(self_d)),
               "median_d_avg": float(a.d_avg.median()), "median_d_avg_with_self_old": float(a.d_avg_with_self_old.median()),
               "max_abs_change": float((a.d_avg - a.d_avg_with_self_old).abs().max()),
               "plant_sequences": int(m), "sets_rows": b.groupby("set").size().to_dict(),
               "sets_n_identical": b.groupby("set").n_identical.sum().to_dict(),
               "sets_median_d_mean_plant": b.groupby("set").d_mean_plant.median().to_dict(),
               "check_with_identical_vs_family_metrics_max_abs":
                   float((ok.d_mean_plant_with_identical - ok.d_homolog_mean).abs().max())}
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2), flush=True)
    (OUT / "DAVGX_DONE").write_text("done\n")
    log("DAVGX_DONE")


if __name__ == "__main__":
    main()
