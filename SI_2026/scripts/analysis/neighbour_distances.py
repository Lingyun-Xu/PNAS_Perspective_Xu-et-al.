#!/usr/bin/env python3
"""Cutoff-free neighbour measures (author request, 6 Nov 2026). Runs on the cluster, next to davg_change.py.

Pairwise distance: the Euclidean distance between the ESM-2 (650M) mean embeddings of two proteins, the distance that
is averaged in the average homolog distance.
  1) The proteins added to the template library in Fig. S13 (pdb3825, sel3800, sel690, and cf3800, cf690) against the
     candidate pool (179,596 unique sequences): the number of pool candidates within pairwise distance r, for every r
     from 0.1 to 5.0 in steps of 0.1, and the distance to the nearest other pool candidate (a pool protein itself is not
     counted).
  2) Every protein of the compared sets of Section S7.5 (the queries of pdb3825/family_metrics.py) against the 182,286
     unique plant sequences of the AFDB v6 set: the distance to the nearest other plant sequence (identical sequences
     not counted, as in family_metrics.py), and the number within 1.5 as a check against family_metrics.tsv.
Outputs (campaign/neighbour_distances/): added_radii.tsv, sets_nearest_plant.tsv, summary.json
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

OUT = cs.CAMP / "neighbour_distances"
RADII = np.round(np.arange(0.1, 5.0001, 0.1), 1)


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def radii_counts(x, own_slot, pool, pool_norm, chunk=256):
    """per row of x: the number of pool rows with distance < r for every r in RADII, and the nearest distance; the row's
    own slot in the pool (or -1) is left out"""
    k = len(RADII) + 1
    counts = np.zeros((len(x), len(RADII)), dtype=np.int64)
    nearest = np.zeros(len(x), dtype=np.float64)
    xn = cs.sq_norm(x)
    for s in range(0, len(x), chunk):
        e = min(len(x), s + chunk)
        d = cs.distances(x[s:e], pool, xn[s:e], pool_norm)
        for i, sl in enumerate(own_slot[s:e]):
            if sl >= 0:
                d[i, sl] = np.inf
        nearest[s:e] = d.min(axis=1)
        idx = np.searchsorted(RADII, d, side="right")          # number of radii <= d; d < RADII[j] if idx <= j
        idx[~np.isfinite(d)] = k - 1
        flat = (idx + (np.arange(e - s)[:, None] * k)).ravel()
        hist = np.bincount(flat, minlength=(e - s) * k).reshape(e - s, k)
        counts[s:e] = np.cumsum(hist, axis=1)[:, :len(RADII)]
    return counts, nearest


def part1():
    rep, accs = sp.load_pool()
    emb, index = pl.load_subset(sp.AF_DIRS, set(rep.seq_id))
    meta = sp.align_meta(rep, accs, index)
    slot = {sid: i for i, sid in enumerate(meta["seq_id"])}
    norm = cs.sq_norm(emb)
    sets = {}
    for tag, folder in (("sel3800", "selection_3800"), ("sel690", "selection_690")):
        sets[tag] = np.array(sorted({s for v in cs.saved_selection(folder, slot).values() for s in v}), dtype=np.int64)
    for tag, name in (("cf3800", "cf_sel3800"), ("cf690", "cf_sel690")):
        t = pd.read_csv(cs.OUT / f"{name}_selected.tsv", sep="\t", dtype={"seq_id": str})
        sets[tag] = np.array(sorted(slot[s] for s in t.seq_id), dtype=np.int64)
    frames = []
    for tag, rows in sets.items():
        counts, nearest = radii_counts(emb[rows], rows, emb, norm)
        f = pd.DataFrame(counts, columns=[f"n_within_{r:.1f}" for r in RADII])
        f.insert(0, "d_nearest_pool", nearest)
        f.insert(0, "id", [meta["seq_id"][r] for r in rows])
        f.insert(0, "set", tag)
        frames.append(f)
        log(f"part 1 {tag}: {len(rows):,}")
    rec = pd.read_csv(cs.OUT / "recent_plant_pdb_sequences.tsv", sep="\t", dtype={"embed_id": str})
    pemb, pindex = pl.load_subset(sp.PDB_DIRS, set(rec.embed_id))
    ids = [e for e in rec.embed_id if e in pindex]
    x = np.asarray(pemb[[pindex[e] for e in ids]], dtype=np.float32)
    counts, nearest = radii_counts(x, np.full(len(ids), -1), emb, norm)
    f = pd.DataFrame(counts, columns=[f"n_within_{r:.1f}" for r in RADII])
    f.insert(0, "d_nearest_pool", nearest)
    f.insert(0, "id", ids)
    f.insert(0, "set", "pdb3825")
    frames.append(f)
    log(f"part 1 pdb3825: {len(ids):,}")
    out = pd.concat(frames, ignore_index=True)
    out.to_csv(OUT / "added_radii.tsv", sep="\t", index=False)
    return out


def part2():
    af = pd.read_parquet(fm.ROOT / "af_queries.parquet",
                         columns=["db_accession", "seq_id", "sequence", "tax_group7", "organism_scientific", "in_pdb",
                                  "plddt_global"])
    af["seq_id"] = af.seq_id.astype(str)
    af["db_accession"] = af.db_accession.astype(str)
    plants = af[af.tax_group7.eq("Viridiplantae")].drop_duplicates("seq_id").copy()
    plants["md5"] = plants.sequence.map(fm.md5)
    log(f"part 2: plant sequences {len(plants):,}")
    q = fm.queries(af)
    lib, lib_index = pl.load_subset(sp.AF_DIRS, set(plants.seq_id))
    lib_md5 = np.array([None] * len(lib), dtype=object)
    for sid, h in zip(plants.seq_id, plants.md5):
        if sid in lib_index:
            lib_md5[lib_index[sid]] = h
    lib_norm = cs.sq_norm(lib)
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
            other = d[i][lib_md5 != h]
            rows.append((st, qid, float(other.min()), int((other < 1.5).sum())))
    out = pd.DataFrame(rows, columns=["set", "query_id", "d_nearest_plant", "n_within_1.5_check"])
    old = pd.read_csv(fm.OUT / "family_metrics.tsv", sep="\t", dtype={"query_id": str}, usecols=["set", "query_id", "n_within_1.5"])
    out = out.merge(old, on=["set", "query_id"], how="left")
    out.to_csv(OUT / "sets_nearest_plant.tsv", sep="\t", index=False)
    log(f"part 2: {len(out):,} proteins with an embedding")
    return out, q


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    a = part1()
    b, q = part2()
    agree = b.dropna(subset=["n_within_1.5"])
    summary = {"part1_rows": a.groupby("set").size().to_dict(),
               "part1_median_d_nearest_pool": a.groupby("set").d_nearest_pool.median().to_dict(),
               "part1_median_n_within_1.5": a.groupby("set")["n_within_1.5"].median().to_dict(),
               "part2_rows": b.groupby("set").size().to_dict(), "part2_queries": q.groupby("set").size().to_dict(),
               "part2_median_d_nearest_plant": b.groupby("set").d_nearest_plant.median().to_dict(),
               "check_n_within_1.5_equal_share": float((agree["n_within_1.5"] == agree["n_within_1.5_check"]).mean()),
               "check_n_within_1.5_max_abs_diff": float((agree["n_within_1.5"] - agree["n_within_1.5_check"]).abs().max())}
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2), flush=True)
    (OUT / "NBR_DONE").write_text("done\n")
    log("NBR_DONE")


if __name__ == "__main__":
    main()
