#!/usr/bin/env python3
"""Counterfactual comparison: real recent plant PDB structures versus AI-selected plant proteins.

Comparison 3800
    baseline library L0 = September 2026 PDB embedding library without the 3,825 plant sequences
    first released after 2021-02-03 (pdb3825, distances/plant_library/excluded_embed_ids.txt).
    real arm = L0 + pdb3825.  AI arm = L0 + 3,800 proteins chosen on L0 (38 rounds x 100).
Comparison 690
    baseline library L1 = full library without the 690 most recently first-released plant
    sequences (recent690, about the last year).  real arm = L1 + recent690.
    AI arm = L1 + 690 proteins chosen on L1 (69 rounds x 10).

Both arms are scored on the same plant proteins: the no-structure pool minus every protein picked
by any selection (original or counterfactual), so no evaluated protein is itself a new template in
either arm. Distances are raw Euclidean minima between ESM-2 mean embeddings, as in the selection.

The selection is the greedy rule of select_plants.py, computed incrementally: the pool sum of
distances and the library minimum are computed once, then updated with the distances to each
round's picks. Runs 1-2 repeat the original sel3800 and sel690 on the full library and are
compared with the saved selections as a check of the incremental computation.
"""
import json
import sys
import time
from datetime import date
from pathlib import Path

sys.path.insert(0, "/path/to/cluster/redo_PDB/ESM_embedding/template_v6")
sys.path.insert(0, "/path/to/cluster/redo_PDB/ESM_embedding/template_v6/campaign")

import numpy as np
import pandas as pd

import compute_distances as cd
import plant_library_distances as pl
import select_plants as sp

ROOT = sp.ROOT
CAMP = ROOT / "campaign"
OUT = CAMP / "counterfactual"
CHUNK = 2048
WALL = date(2021, 2, 3)
N_RECENT = 690
THRESHOLDS = (1.7, 1.9)


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def distances(q, lib, q_norm, lib_norm):
    d2 = q @ lib.T
    d2 *= np.float32(-2)
    d2 += q_norm[:, None]
    d2 += lib_norm
    np.maximum(d2, 0, out=d2)
    return np.sqrt(d2, out=d2)


def sq_norm(x):
    return np.einsum("ij,ij->i", x, x).astype(np.float32)


def min_to(query, lib, tag):
    out = np.full(len(query), np.inf, dtype=np.float32)
    if len(lib) == 0:
        return out
    qn, ln = sq_norm(query), sq_norm(lib)
    for s in range(0, len(query), CHUNK):
        e = min(len(query), s + CHUNK)
        out[s:e] = distances(query[s:e], lib, qn[s:e], ln).min(axis=1)
    log(f"min {tag}: {len(query):,} x {len(lib):,}")
    return out


def pool_sums(emb):
    n = len(emb)
    norm = sq_norm(emb)
    out = np.empty(n, dtype=np.float64)
    for k, s in enumerate(range(0, n, CHUNK)):
        e = min(n, s + CHUNK)
        out[s:e] = distances(emb[s:e], emb, norm[s:e], norm).sum(axis=1, dtype=np.float64)
        if k % 10 == 0:
            log(f"pool sums {e:,}/{n:,}")
    return out


def recent_plant_sequences(pdb_index):
    """pdb3825 with first plant-chain release dates (same rule as plant_library_distances.py)."""
    kinds = pd.read_csv(ROOT / "plant_org_kind.tsv", sep="\t", dtype=str).fillna("")
    plant_org = set(kinds.loc[kinds.kind.isin(["land", "algae"]), "org"])
    chains = pd.read_parquet(ROOT / "pdb_fasta/pdb_chain_map.parquet",
                             columns=["embed_id", "source_organism", "initial_release_date"])
    org = chains.source_organism.fillna("").astype(str).str.strip()
    plants = chains.loc[org.isin(plant_org)].copy()
    plants["day"] = plants.initial_release_date.map(pl.parse_day)
    plants = plants[plants.day.notna() & (plants.day <= date(2026, 9, 30))]
    first = plants.groupby("embed_id").day.min()
    recent = first[first > WALL]
    saved = set((ROOT / "distances/plant_library/excluded_embed_ids.txt").read_text().split())
    if set(recent.index.astype(str)) != saved:
        raise SystemExit(f"recomputed recent set ({len(recent)}) differs from excluded_embed_ids.txt ({len(saved)})")
    table = pd.DataFrame({"embed_id": recent.index.astype(str), "first_release": recent.values})
    missing = [e for e in table.embed_id if e not in pdb_index]
    if missing:
        raise SystemExit(f"{len(missing)} recent plant sequences have no PDB embedding, first {missing[:5]}")
    table = table.sort_values(["first_release", "embed_id"], ascending=[False, True]).reset_index(drop=True)
    table["recent690"] = table.index < N_RECENT
    table["row"] = table.embed_id.map(pdb_index)
    return table


def greedy(name, emb, norm, plddt, s0, min0, rounds, per_round):
    n = len(emb)
    alive = np.ones(n, dtype=bool)
    sums = s0.copy()
    mins = min0.astype(np.float32).copy()
    n_alive = n
    target = rounds * per_round
    rows = []
    for rnd in range(rounds):
        take = min(per_round, target - (n - n_alive), n_alive)
        if take <= 0:
            break
        slots = np.flatnonzero(alive)
        avg = (sums[slots] / n_alive).astype(np.float32)
        score = sp.score_rows(mins[slots], avg, plddt[slots])
        picked = np.argsort(-score)[:take]
        chosen = slots[picked]
        for rank, j in enumerate(picked, start=1):
            rows.append((rnd, rank, int(slots[j]), float(score[j]), float(mins[slots[j]]), float(avg[j])))
        alive[chosen] = False
        n_alive -= len(chosen)
        rest = np.flatnonzero(alive)
        picks, pick_norm = emb[chosen], norm[chosen]
        for s in range(0, len(rest), 8192):
            idx = rest[s:s + 8192]
            d = distances(emb[idx], picks, norm[idx], pick_norm)
            sums[idx] -= d.sum(axis=1, dtype=np.float64)
            mins[idx] = np.minimum(mins[idx], d.min(axis=1))
    log(f"greedy {name}: {len(rows)} picks in {rnd + 1} rounds")
    return pd.DataFrame(rows, columns=["round", "rank_in_round", "slot", "score", "min_dist", "avg_dist"])


def saved_selection(folder, seq_slot):
    rounds = {}
    for path in sorted((CAMP / folder).glob("round_*/selected_seq_ids.txt")):
        rnd = int(path.parent.name.split("_")[1])
        rounds[rnd] = [seq_slot[s] for s in path.read_text().split()]
    return rounds


def check_reproduction(name, run, saved):
    got = run.groupby("round").slot.apply(list).to_dict()
    same_rounds = sum(set(got.get(r, [])) == set(v) for r, v in saved.items())
    all_got = set(run.slot)
    all_saved = {s for v in saved.values() for s in v}
    overlap = len(all_got & all_saved)
    log(f"check {name}: identical rounds {same_rounds}/{len(saved)}  overlap {overlap}/{len(all_saved)}")
    return {"rounds_identical": same_rounds, "rounds": len(saved), "overlap": overlap, "n_saved": len(all_saved)}


def arm_metrics(label, base, arm):
    gain = base - arm
    row = {"arm": label, "n": int(len(base)), "mean_dist": float(arm.mean()), "median_dist": float(np.median(arm)),
           "frac_closer": float((gain > 1e-4).mean()), "frac_closer_by_0.1": float((gain >= 0.1).mean()),
           "frac_closer_by_0.5": float((gain >= 0.5).mean()),
           "median_gain_if_closer": float(np.median(gain[gain > 1e-4])) if (gain > 1e-4).any() else 0.0,
           "mean_gain": float(gain.mean())}
    for t in THRESHOLDS:
        row[f"frac_within_{t}"] = float((arm < t).mean())
    return row


def genus_group(org):
    g = str(org).split()[0] if str(org).split() else ""
    return g if g in ("Oryza", "Zea", "Glycine", "Arabidopsis") else "other"


def plot_comparison(tag, n_new_real, n_new_ai, base, real, ai, homolog, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    xlim = (0.0, float(np.percentile(base, 99.9)) * 1.05)
    ylim = (float(np.percentile(homolog, 0.1)) * 0.97, float(np.percentile(homolog, 99.9)) * 1.03)
    fig, axes = plt.subplots(2, 3, figsize=(16, 9.6))
    panels = [(base, "a  Baseline library"),
              (real, f"b  + {n_new_real:,} real plant PDB sequences"),
              (ai, f"c  + {n_new_ai:,} AI-selected plant proteins")]
    for ax, (x, name) in zip(axes[0], panels):
        h = ax.hist2d(x, homolog, bins=150, range=[xlim, ylim], cmap="plasma_r", norm=LogNorm(), cmin=1)
        ax.set_title(name, loc="left", fontsize=11)
        ax.set_xlabel("Nearest template distance")
        ax.set_ylabel("Average homolog distance")
        fig.colorbar(h[3], ax=ax, fraction=0.046, pad=0.03, label="proteins per bin")
    for ax, x, name in ((axes[1][0], real, "d  Real minus baseline"), (axes[1][1], ai, "e  AI-selected minus baseline")):
        cf = sp.free_energy_difference(x, homolog, base, homolog, ax, xlim, ylim)
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_title(name, loc="left", fontsize=11)
        ax.set_xlabel("Nearest template distance")
        ax.set_ylabel("Average homolog distance")
        fig.colorbar(cf, ax=ax, fraction=0.046, pad=0.03, label="Free-energy difference (red loss, blue gain)")
    ax = axes[1][2]
    grid = np.linspace(0, 2.0, 201)
    for x, name, color in ((real, "real", "#555555"), (ai, "AI-selected", "#d95f02")):
        gain = base - x
        share = [(gain > 1e-4).mean() * 100] + [(gain >= g).mean() * 100 for g in grid[1:]]
        ax.plot(grid, share, color=color, lw=2, label=name)
    ax.set_xlabel("Decrease in nearest template distance (at least)")
    ax.set_ylabel("Evaluated plant proteins (%)")
    ax.set_title("f  Plant proteins brought closer to a template", loc="left", fontsize=11)
    ax.legend(frameon=False)
    ax.set_xlim(0, 2.0)
    fig.suptitle(title, fontsize=12)
    fig.tight_layout()
    fig.savefig(OUT / f"real_vs_ai_{tag}.png", dpi=200, bbox_inches="tight")
    fig.savefig(OUT / f"real_vs_ai_{tag}.pdf", bbox_inches="tight")
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rep, accs = sp.load_pool()
    plant_emb, plant_index = pl.load_subset(sp.AF_DIRS, set(rep.seq_id))
    meta = sp.align_meta(rep, accs, plant_index)
    seq_slot = {sid: i for i, sid in enumerate(meta["seq_id"])}
    plddt = meta["plddt"]
    norm = sq_norm(plant_emb)
    log(f"plant pool embeddings {plant_emb.shape}")
    pdb_emb, pdb_index = cd.load_emb(sp.PDB_DIRS)
    log(f"PDB embeddings {pdb_emb.shape}")

    recent = recent_plant_sequences(pdb_index)
    recent.to_csv(OUT / "recent_plant_pdb_sequences.tsv", sep="\t", index=False)
    r690 = recent[recent.recent690]
    log(f"pdb3825 {len(recent):,}  first release {recent.first_release.min()} .. {recent.first_release.max()}")
    log(f"recent690 first release {r690.first_release.min()} .. {r690.first_release.max()}")
    in_recent = np.zeros(len(pdb_emb), dtype=bool)
    in_recent[recent.row.to_numpy()] = True
    rows_690 = r690.row.to_numpy()
    rows_rest = recent.loc[~recent.recent690, "row"].to_numpy()

    s0_path = OUT / "pool_sums_round0.npy"
    if s0_path.exists():
        s0 = np.load(s0_path)
        log("pool sums loaded")
    else:
        s0 = pool_sums(plant_emb)
        np.save(s0_path, s0)
    mins_path = OUT / "mins_round0.npz"
    if mins_path.exists():
        z = np.load(mins_path)
        d_l0, d_690, d_rest = z["d_l0"], z["d_690"], z["d_rest"]
        log("round-0 minima loaded")
    else:
        d_l0 = min_to(plant_emb, pdb_emb[~in_recent], "pool -> L0")
        d_690 = min_to(plant_emb, pdb_emb[rows_690], "pool -> recent690")
        d_rest = min_to(plant_emb, pdb_emb[rows_rest], "pool -> pdb3825 minus recent690")
        np.savez(mins_path, d_l0=d_l0, d_690=d_690, d_rest=d_rest)
    del pdb_emb
    d_full = np.minimum(np.minimum(d_l0, d_690), d_rest)
    d_l1 = np.minimum(d_l0, d_rest)

    checks = {}
    saved_min = np.load(CAMP / "selection_3800/round_000/min_dist.npy")
    saved_avg = np.load(CAMP / "selection_3800/round_000/avg_dist.npy")
    if len(saved_min) == len(d_full):
        checks["round0_min_max_abs_diff"] = float(np.max(np.abs(saved_min - d_full)))
        checks["round0_avg_max_rel_diff"] = float(np.max(np.abs(saved_avg - s0 / len(s0)) / saved_avg))
    log(f"round-0 check {checks}")

    runs = {
        "repro_sel3800": (d_full, 38, 100),
        "repro_sel690": (d_full, 69, 10),
        "cf_sel3800": (d_l0, 38, 100),
        "cf_sel690": (d_l1, 69, 10),
    }
    picks = {}
    for name, (min0, rounds, per_round) in runs.items():
        run = greedy(name, plant_emb, norm, plddt, s0, min0, rounds, per_round)
        run["seq_id"] = [meta["seq_id"][s] for s in run.slot]
        run["representative_accession"] = [meta["rep_acc"][s] for s in run.slot]
        run["organism_scientific"] = [meta["organism"][s] for s in run.slot]
        run["plddt"] = plddt[run.slot.to_numpy()]
        run.drop(columns="slot").to_csv(OUT / f"{name}_selected.tsv", sep="\t", index=False)
        picks[name] = run
    for name, folder in (("repro_sel3800", "selection_3800"), ("repro_sel690", "selection_690")):
        checks[name] = check_reproduction(name, picks[name], saved_selection(folder, seq_slot))
    orig3800 = {s for v in saved_selection("selection_3800", seq_slot).values() for s in v}
    orig690 = {s for v in saved_selection("selection_690", seq_slot).values() for s in v}
    cf3800, cf690 = set(picks["cf_sel3800"].slot), set(picks["cf_sel690"].slot)
    overlaps = {"cf3800_vs_orig3800": len(cf3800 & orig3800), "cf690_vs_orig690": len(cf690 & orig690),
                "cf690_in_cf3800": len(cf690 & cf3800)}
    log(f"overlaps {overlaps}")

    excluded = orig3800 | orig690 | cf3800 | cf690 | set(picks["repro_sel3800"].slot) | set(picks["repro_sel690"].slot)
    ev = np.array(sorted(set(range(len(plant_emb))) - excluded), dtype=np.int64)
    q = plant_emb[ev]
    to_cf3800 = min_to(q, plant_emb[sorted(cf3800)], "eval -> cf3800")
    to_cf690 = min_to(q, plant_emb[sorted(cf690)], "eval -> cf690")
    to_o3800 = min_to(q, plant_emb[sorted(orig3800)], "eval -> orig3800")
    to_o690 = min_to(q, plant_emb[sorted(orig690)], "eval -> orig690")
    table = pd.DataFrame({
        "seq_id": [meta["seq_id"][s] for s in ev],
        "representative_accession": [meta["rep_acc"][s] for s in ev],
        "organism_scientific": [meta["organism"][s] for s in ev],
        "plddt": plddt[ev],
        "base_3800": d_l0[ev],
        "real_3800": d_full[ev],
        "ai_3800": np.minimum(d_l0[ev], to_cf3800),
        "base_690": d_l1[ev],
        "real_690": d_full[ev],
        "ai_690": np.minimum(d_l1[ev], to_cf690),
        "full_2026": d_full[ev],
        "full_plus_sel3800": np.minimum(d_full[ev], to_o3800),
        "full_plus_sel690": np.minimum(d_full[ev], to_o690),
    })
    hom = pd.read_parquet(ROOT / "distances/plant_library/plant_library_distances.parquet", columns=["seq_id", "d_homolog"])
    hom = hom.assign(seq_id=hom.seq_id.astype(str)).drop_duplicates("seq_id").set_index("seq_id").d_homolog
    table["d_homolog"] = table.seq_id.map(hom).astype(np.float32)
    table["genus_group"] = table.organism_scientific.map(genus_group)
    table.to_parquet(OUT / "evaluation_distances.parquet", index=False)
    log(f"evaluation set {len(table):,} proteins (homolog distance missing for {int(table.d_homolog.isna().sum())})")

    rows = []
    for comp, base, arms in (
            ("3800", "base_3800", [("real: +3,825 plant PDB sequences released after 2021-02-03", "real_3800"),
                                   ("AI: +3,800 selected on the baseline library", "ai_3800")]),
            ("690", "base_690", [("real: +690 most recent plant PDB sequences", "real_690"),
                                 ("AI: +690 selected on the baseline library", "ai_690")]),
            ("on_top_2026", "full_2026", [("+ original sel3800", "full_plus_sel3800"),
                                          ("+ original sel690", "full_plus_sel690")])):
        b = table[base].to_numpy()
        rows.append({"comparison": comp, **arm_metrics("baseline", b, b)})
        for label, col in arms:
            rows.append({"comparison": comp, **arm_metrics(label, b, table[col].to_numpy())})
            for grp, sub in table.groupby("genus_group"):
                rows.append({"comparison": comp, "group": grp,
                             **arm_metrics(label, sub[base].to_numpy(), sub[col].to_numpy())})
    metrics = pd.DataFrame(rows)
    metrics.to_csv(OUT / "real_vs_ai_metrics.tsv", sep="\t", index=False)
    near_identical = {c: int((table[c] < 0.05).sum()) for c in ("base_3800", "real_3800", "ai_3800", "base_690", "ai_690")}

    ok = table.d_homolog.notna().to_numpy()
    plot_comparison("3800", len(recent), len(cf3800), table.base_3800.to_numpy()[ok], table.real_3800.to_numpy()[ok],
                    table.ai_3800.to_numpy()[ok], table.d_homolog.to_numpy()[ok],
                    f"Plant proteins without a structure (n = {int(ok.sum()):,}): PDB without the 3,825 recent plant "
                    "sequences, adding the real ones (b, d) or 3,800 AI-selected proteins (c, e)")
    plot_comparison("690", len(r690), len(cf690), table.base_690.to_numpy()[ok], table.real_690.to_numpy()[ok],
                    table.ai_690.to_numpy()[ok], table.d_homolog.to_numpy()[ok],
                    f"Plant proteins without a structure (n = {int(ok.sum()):,}): PDB without the 690 most recent plant "
                    "sequences, adding the real ones (b, d) or 690 AI-selected proteins (c, e)")
    summary = {"checks": checks, "overlaps": overlaps, "n_eval": int(len(table)), "near_identical_lt_0.05": near_identical,
               "pdb3825_first_release": [str(recent.first_release.min()), str(recent.first_release.max())],
               "recent690_first_release": [str(r690.first_release.min()), str(r690.first_release.max())]}
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
    print(metrics[metrics.group.isna()].to_string(index=False), flush=True)
    print(json.dumps(summary, indent=2), flush=True)
    (OUT / "CF_DONE").write_text("done\n")
    log("CF_DONE")


if __name__ == "__main__":
    main()
