#!/usr/bin/env python3
"""Main-text Fig. 4 recomputed on the AFDB v6 plant proteins of main-text Fig. 3.

Candidate pool: the 179,596 unique sequences of the 182,501 AFDB v6 Viridiplantae accessions without a SIFTS-mapped
PDB chain (select_plants.load_pool, as main-text Fig. 3 and the selections); ESM-2 650M mean embeddings (layer 33,
residues only, truncation at 1,022 tokens) of template_v6.
Template pools, one run each:
  thermostability  FireProtDB wild types, mega-scale Dataset 2/3 sequences and Meltome Atlas proteins
                   (inputs/thermostability_library.pt: the February 2025 embeddings, packed with their ids)
  mavedb           the 1,170 unique MaveDB target sequences (November 2024 export), v6 extractor (emb_mavedb)
  mdcath           the 5,398 mdCATH domains (CATH v4.2.0 domain sequences), v6 extractor (emb_mdcath)
Selection as in the first version of Fig. 4 (analysis_avg_vs_min_L2_scaled.py, no pLDDT term): score =
(d_min / ||d_min||) / (d_avg / ||d_avg||), whose ranking equals that of d_min / d_avg; the top N candidates of each
round move from the candidate pool into the template pool and all distances are updated (incrementally, as
counterfactual_selection.greedy: the pool sums and template minima are computed once and updated with the distances
to each round's picks). Saved for rounds 0 (template pool alone), 1 (after the first round) and R (after all R
rounds: 20 x 500 thermostability and mdCATH, 40 x 50 MaveDB): the candidates still in the pool, d_min (distance to the
nearest sequence of the template pool) and d_avg (mean distance to the other candidates still in the pool, own zero
distance excluded, as in main-text Fig. 3). Usage: python fig4_v6_selection.py [run ...]"""
import json
import sys
import time
from pathlib import Path

ROOT = Path("/path/to/cluster/redo_PDB/ESM_embedding/template_v6")
CAMP = ROOT / "campaign"
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(CAMP))

import numpy as np
import pandas as pd
import torch

import compute_distances as cd
import plant_library_distances as pl
import select_plants as sp
from counterfactual_selection import distances, min_to, sq_norm

W = CAMP / "fig4_v6"
OUT = W / "results"
RUNS = {"thermostability": (20, 500), "thermostability_2025": (20, 500), "mavedb": (40, 50), "mdcath": (20, 500)}


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def min_blocked(query, lib, tag, qchunk=2048, lchunk=131072):
    """minimum distance from every query row to the template pool, in blocks of both (memory ~1 GB per block)"""
    out = np.full(len(query), np.inf, dtype=np.float32)
    qn, ln = sq_norm(query), sq_norm(lib)
    for k, s in enumerate(range(0, len(query), qchunk)):
        e = min(len(query), s + qchunk)
        for t in range(0, len(lib), lchunk):
            u = min(len(lib), t + lchunk)
            out[s:e] = np.minimum(out[s:e], distances(query[s:e], lib[t:u], qn[s:e], ln[t:u]).min(axis=1))
        if k % 10 == 0:
            log(f"  min {tag}: {e:,}/{len(query):,}")
    return out


def load_library(tag):
    if tag.startswith("thermostability"):
        obj = torch.load(W / "inputs/thermostability_library_2025.pt", map_location="cpu", weights_only=False)
        emb, ids = obj["emb"].float().numpy(), list(obj["ids"])
        if tag == "thermostability":                 # + every other protein sequence of mega-scale Dataset 1
            # embedded in 10 parts; slow parts were re-split into _s0/_s1 folders, so take every folder with shards
            dirs = sorted(str(d) for d in W.glob("emb_megaD1_p*") if d.is_dir() and any(d.glob("shard_*.pt")))
            log(f"Dataset 1 embedding folders: {len(dirs)}")
            d1, idx = cd.load_emb(dirs)
            if len(idx) != 996737:
                raise SystemExit(f"expected 996,737 Dataset 1 embeddings, found {len(idx)}")
            emb = np.concatenate([emb, d1])
            ids += sorted(idx, key=idx.get)
        return emb, ids
    folder = "emb_mdcath_atom" if tag == "mdcath" else f"emb_{tag}"      # mdCATH: CATH S20 ATOM-record sequences
    emb, index = cd.load_emb([str(W / folder)])
    return emb, sorted(index, key=index.get)


def greedy_states(emb, norm, s0, min0, rounds, per_round, meta):
    n = len(emb)
    alive = np.ones(n, dtype=bool)
    sums = s0.copy()
    mins = min0.astype(np.float32).copy()
    n_alive = n
    states, picks = {}, []

    def record(r):
        slots = np.flatnonzero(alive)
        states[r] = {"slots": slots.astype(np.int32), "d_min": mins[slots].copy(),
                     "d_avg": (sums[slots] / (n_alive - 1)).astype(np.float32)}

    for rnd in range(rounds):
        if rnd in (0, 1):
            record(rnd)
        slots = np.flatnonzero(alive)
        avg = sums[slots] / n_alive                              # the selection code's mean (own zero included)
        score = mins[slots].astype(np.float64) / avg             # same ranking as the L2-normalised ratio
        picked = np.argsort(-score, kind="stable")[:per_round]
        chosen = slots[picked]
        for rank, j in enumerate(picked, start=1):
            s = int(slots[j])
            picks.append((rnd, rank, meta["seq_id"][s], meta["rep_acc"][s], meta["organism"][s],
                          float(meta["plddt"][s]), float(score[j]), float(mins[s]),
                          float(sums[s] / (n_alive - 1))))
        alive[chosen] = False
        n_alive -= len(chosen)
        rest = np.flatnonzero(alive)
        p, pn = emb[chosen], norm[chosen]
        for s in range(0, len(rest), 8192):
            idx = rest[s:s + 8192]
            d = distances(emb[idx], p, norm[idx], pn)
            sums[idx] -= d.sum(axis=1, dtype=np.float64)
            mins[idx] = np.minimum(mins[idx], d.min(axis=1))
        log(f"  round {rnd}: picked {len(chosen)}, pool {n_alive:,}, median d_min of picks "
            f"{float(np.median(mins[chosen])):.3f}")
    record(rounds)
    cols = ["round", "rank_in_round", "seq_id", "representative_accession", "organism_scientific", "plddt", "score",
            "d_min", "d_avg"]
    return states, pd.DataFrame(picks, columns=cols)


def main():
    todo = sys.argv[1:] or list(RUNS)
    OUT.mkdir(parents=True, exist_ok=True)
    rep, accs = sp.load_pool()
    plant_emb, plant_index = pl.load_subset(sp.AF_DIRS, set(rep.seq_id))
    meta = sp.align_meta(rep, accs, plant_index)
    n = len(plant_emb)
    log(f"candidate pool embeddings {plant_emb.shape}")
    if n != 179596:
        raise SystemExit(f"expected 179,596 candidates, found {n}")
    norm = sq_norm(plant_emb)
    s0 = np.load(CAMP / "counterfactual/pool_sums_round0.npy")
    rows = np.random.default_rng(0).choice(n, 32, replace=False)
    chk = distances(plant_emb[rows], plant_emb, norm[rows], norm).sum(axis=1, dtype=np.float64)
    rel = float(np.max(np.abs(chk - s0[rows]) / s0[rows]))
    log(f"pool sums of the selections reused; check on 32 rows: max relative difference {rel:.2e}")
    if len(s0) != n or rel > 1e-5:
        raise SystemExit("pool sums do not match this candidate order")
    pd.DataFrame({"slot": np.arange(n), "seq_id": meta["seq_id"], "representative_accession": meta["rep_acc"],
                  "organism_scientific": meta["organism"], "plddt": meta["plddt"]}).to_csv(
        OUT / "candidates.tsv.gz", sep="\t", index=False)
    summary = {}
    for tag in todo:
        rounds, per_round = RUNS[tag]
        lib, ids = load_library(tag)
        log(f"{tag}: template pool {lib.shape}, {len(set(ids)):,} ids")
        min0 = min_blocked(plant_emb, lib, f"candidates -> {tag}")
        if not tag.startswith("thermostability"):
            chk = min_to(plant_emb[:4096], lib, f"check {tag}")
            log(f"{tag}: blocked vs unblocked minimum, max difference {float(np.max(np.abs(chk - min0[:4096]))):.2e}")
        states, picks = greedy_states(plant_emb, norm, s0, min0, rounds, per_round, meta)
        np.savez_compressed(OUT / f"{tag}_states.npz", **{f"{k}_r{r}": v for r, st in states.items() for k, v in st.items()})
        picks.to_csv(OUT / f"{tag}_selected.tsv", sep="\t", index=False)
        summary[tag] = {"template_pool": int(len(lib)), "rounds": rounds, "per_round": per_round,
                        **{f"round_{r}": {"candidates": int(len(st["slots"])),
                                          "median_d_min": float(np.median(st["d_min"])),
                                          "median_d_avg": float(np.median(st["d_avg"])),
                                          "d_min_le_1.9_pct": float((st["d_min"] <= 1.9).mean() * 100)}
                           for r, st in states.items()}}
        log(f"{tag}: {json.dumps(summary[tag])}")
        (OUT / f"summary_{tag}.json").write_text(json.dumps(summary[tag], indent=1))
        del lib
    (OUT / f"DONE_{'_'.join(todo)}").write_text("done\n")
    log("FIG4_DONE")


if __name__ == "__main__":
    main()
