#!/usr/bin/env python3
"""Apply the corrected average homolog distance (d_avg; the candidate's own zero distance excluded) to the source data of
the SI and of main-text Fig. 3c (author request, 6 Nov 2026). Input: davg_self_excluded.py (cluster) outputs.

d_avg(i) = sum over the other candidates j of d(i, j) / (N - 1), N = candidates in the candidate pool at that point.
  FigS11 candidate_pool_round0.tsv.gz       d_average_homolog <- d_avg of round 0 (by pool_position)
  FigS12 *_rounds_*.tsv.gz                  round 0 <- d_avg (by pool_position); last round: x N_r / (N_r - 1)
  FigS12 selected_accessions_*.tsv          avg_dist: round 0 <- d_avg (by seq_id); later rounds x N_r / (N_r - 1)
  FigS13 evaluation_distances.parquet       d_average_homolog <- d_avg (by seq_id)
  FigS13 davg_change.tsv.gz (new)           d_avg of the 175,537 evaluated candidates with the full candidate pool and
                                            after the proteins of sel3800 or sel690 leave it, and the changes (Fig. S13)
  FigS15_S18 cohort and pairs tables        d_homolog (mean distance to all 182,286 plant sequences, self included)
                                            -> d_average_homolog (candidate pool, self excluded; by seq_id/accession)
  updated_maintext_figure/data/Fig3c_pool_minus_sel3800.tsv.gz   d_average_homolog <- d_avg (round 0, by seq_id)
  Dataset S15 (both copies)                 + d_mean_plant: mean distance to all plant sequences other than identical
                                            ones (Fig. S21J), with its tests in comparison_statistics.tsv
For the last rounds of FigS12 the round's own sums are not stored; scaling by N_r / (N_r - 1) removes the zero term
exactly, leaving only the float32 rounding of the self distance (at most 0.013 before division, < 1e-7 after).
Usage: python apply_davg_self_excluded.py <folder with candidate_pool_davg.tsv.gz and sets_mean_plant.tsv> <backup dir>
"""
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

SRC, BAK = Path(sys.argv[1]), Path(sys.argv[2])
B = Path("/path/to/project")
S17 = B / "SI_2026/datasets/Dataset_S17_source_data_SI_figures"
S15 = [B / "SI_2026/pdb3825_vs_selections", B / "SI_2026/datasets/Dataset_S15_compared_protein_sets"]
FIG3C = B / "updated_maintext_figure/data/Fig3c_pool_minus_sel3800.tsv.gz"
N = 179596
BAK.mkdir(parents=True, exist_ok=True)
report = {}


def backup(path):
    dst = BAK / path.relative_to(B).as_posix().replace("/", "__")
    if not dst.exists():
        shutil.copy2(path, dst)


cp = pd.read_csv(SRC / "candidate_pool_davg.tsv.gz", sep="\t", dtype={"seq_id": str})
assert len(cp) == N
by_pos = cp.set_index("pool_position").d_avg
old_by_pos = cp.set_index("pool_position").d_avg_with_self_old
by_seq = cp.set_index("seq_id").d_avg
acc_rows = cp[["seq_id", "accessions"]].assign(acc=cp.accessions.str.split(";")).explode("acc")
by_acc = acc_rows.set_index("acc").seq_id.map(by_seq)
assert by_acc.index.is_unique

# ---------------------------------------------------------------- FigS11
p = S17 / "FigS11/candidate_pool_round0.tsv.gz"
backup(p)
t = pd.read_csv(p, sep="\t")
check = float((t.d_average_homolog - t.pool_position.map(old_by_pos)).abs().max())
assert check < 1e-4, check
t["d_average_homolog"] = t.pool_position.map(by_pos).round(5)
t.to_csv(p, sep="\t", index=False, compression="gzip")
report["FigS11"] = {"rows": len(t), "old_vs_cluster_max_abs": check, "median": float(t.d_average_homolog.median())}

# ---------------------------------------------------------------- FigS12 rounds
for name in ("sel3800_rounds_0_and_37.tsv.gz", "sel690_rounds_0_and_68.tsv.gz"):
    p = S17 / "FigS12" / name
    backup(p)
    t = pd.read_csv(p, sep="\t")
    sizes = t.groupby("round").size().to_dict()
    new = t.d_average_homolog.astype(float).copy()
    r0 = t["round"].eq(0)
    check = float((t.loc[r0, "d_average_homolog"] - t.loc[r0, "pool_position"].map(old_by_pos)).abs().max())
    assert check < 1e-4, (name, check)
    new[r0] = t.loc[r0, "pool_position"].map(by_pos)
    for rnd, n_r in sizes.items():
        if rnd != 0:
            m = t["round"].eq(rnd)
            new[m] = t.loc[m, "d_average_homolog"] * n_r / (n_r - 1)
    t["d_average_homolog"] = new.round(5)
    t.to_csv(p, sep="\t", index=False, compression="gzip")
    report[f"FigS12 {name}"] = {"rows per round": {int(k): int(v) for k, v in sizes.items()}, "round0_old_vs_cluster": check}

# ---------------------------------------------------------------- FigS12 selected proteins
for name, per_round in (("selected_accessions_sel3800.tsv", 100), ("selected_accessions_sel690.tsv", 10)):
    p = S17 / "FigS12" / name
    backup(p)
    t = pd.read_csv(p, sep="\t", dtype={"seq_id": str})
    n_r = N - t["round"] * per_round
    r0 = t["round"].eq(0)
    check = float((t.loc[r0, "avg_dist"] - t.loc[r0, "seq_id"].map(cp.set_index("seq_id").d_avg_with_self_old)).abs().max())
    assert check < 1e-4, (name, check)
    new = t.avg_dist * n_r / (n_r - 1)
    new[r0] = t.loc[r0, "seq_id"].map(by_seq)
    t["avg_dist"] = new
    t.to_csv(p, sep="\t", index=False)
    report[f"FigS12 {name}"] = {"rows": len(t), "round0_old_vs_cluster": check}

# ---------------------------------------------------------------- FigS13
p = S17 / "FigS13/evaluation_distances.parquet"
backup(p)
ev = pd.read_parquet(p)
ev["seq_id"] = ev.seq_id.astype(str)
check = float((ev.d_average_homolog - ev.seq_id.map(cp.set_index("seq_id").d_avg_with_self_old)).abs().max())
assert check < 1e-4, check
ev["d_average_homolog"] = ev.seq_id.map(by_seq).astype(np.float32)
ev.to_parquet(p, index=False)
report["FigS13 evaluation_distances"] = {"rows": len(ev), "old_vs_cluster": check}

old = pd.read_csv(S17 / "FigS13/davg_eval.tsv.gz", sep="\t", dtype={"seq_id": str})
d_new = old.seq_id.map(by_seq).to_numpy(float)
s_other = d_new * (N - 1)
out = pd.DataFrame({"seq_id": old.seq_id, "d_avg": d_new})
for tag, k in (("sel3800", 3800), ("sel690", 690)):
    removed = N * old.d_base.to_numpy(float) - (N - k) * old[f"d_plus_{tag}"].to_numpy(float)
    after = (s_other - removed) / (N - 1 - k)
    out[f"d_avg_without_{tag}"] = after
    out[f"change_{tag}"] = after - d_new
    old_change = (old[f"d_plus_{tag}"] - old.d_base).to_numpy(float)
    report[f"FigS13 change {tag}"] = {"median": float(np.median(after - d_new)), "min": float((after - d_new).min()),
                                     "max": float((after - d_new).max()), "share_increased": float((after > d_new).mean()),
                                     "old_change_max_abs_diff": float(np.abs((after - d_new) - old_change).max())}
out.to_csv(S17 / "FigS13/davg_change.tsv.gz", sep="\t", index=False, compression="gzip")
moved = []
for name in ("davg_eval.tsv.gz", "davg_added.tsv", "added_radii.tsv.gz", "davg_summary.json"):
    f = S17 / "FigS13" / name
    if f.exists():
        shutil.move(str(f), str(BAK / f"FigS13__{name}"))
        moved.append(name)
report["FigS13 moved out of the dataset"] = moved
(S17 / "FigS13/davg_change_summary.json").write_text(json.dumps(
    {"candidate_pool": N, "evaluated": int(len(out)), "median_d_avg": float(np.median(d_new)),
     **{k: v for k, v in report.items() if k.startswith("FigS13 change")}}, indent=2))

# ---------------------------------------------------------------- FigS15-S18
for f in sorted((S17 / "FigS15_S18").glob("*_cohort.tsv")) + sorted((S17 / "FigS15_S18").glob("*_pairs.tsv")):
    backup(f)
    t = pd.read_csv(f, sep="\t", dtype={"seq_id": str, "db_accession": str}, low_memory=False)
    if "seq_id" in t:
        new = t.seq_id.map(by_seq)
    else:
        new = t.db_accession.map(by_acc)
    miss = int(new.isna().sum())
    assert miss == 0, (f.name, miss)
    diff = float((new - t.d_homolog).abs().max())
    cols = list(t.columns)
    t["d_homolog"] = new
    t = t.rename(columns={"d_homolog": "d_average_homolog"})
    t.to_csv(f, sep="\t", index=False)
    report[f"FigS15_S18 {f.name}"] = {"rows": len(t), "max_abs_change_vs_all_plant_with_self": diff}

# ---------------------------------------------------------------- main-text Fig. 3c data
backup(FIG3C)
t = pd.read_csv(FIG3C, sep="\t", dtype={"seq_id": str})
check = float((t.d_average_homolog - t.seq_id.map(cp.set_index("seq_id").d_avg_with_self_old)).abs().max())
assert check < 1e-4, check
t["d_average_homolog"] = t.seq_id.map(by_seq).round(5)
t.to_csv(FIG3C, sep="\t", index=False, compression="gzip")
report["main-text Fig3c data"] = {"rows": len(t), "old_vs_cluster": check}

# ---------------------------------------------------------------- Dataset S15: mean distance to all other plant sequences
mp = pd.read_csv(SRC / "sets_mean_plant.tsv", sep="\t", dtype={"query_id": str})
values = {(s, q): d for s, q, d in zip(mp.set, mp.query_id, mp.d_mean_plant)}
COL, LABEL, AFTER_COL, AFTER_MEASURE = "d_mean_plant", "Mean distance to all other plant sequences", "d_nearest_plant", \
    "Distance to the nearest other plant sequence"
for folder in S15:
    p = folder / "per_protein_comparison.tsv"
    backup(p)
    lines = p.read_text(encoding="utf-8").split("\n")
    end = lines[-1] == ""
    rows = [l.split("\t") for l in (lines[:-1] if end else lines)]
    head = rows[0]
    if COL in head:
        k = head.index(COL)
        rows = [r[:k] + r[k + 1:] for r in rows]
        head = rows[0]
    i_set, i_key, i_after = head.index("set"), head.index("key"), head.index(AFTER_COL)
    new_rows, found = [head[:i_after + 1] + [COL] + head[i_after + 1:]], 0
    for r in rows[1:]:
        v = values.get((r[i_set], r[i_key]))
        found += v is not None
        new_rows.append(r[:i_after + 1] + ["" if v is None else repr(float(v))] + r[i_after + 1:])
    p.write_text("\n".join("\t".join(r) for r in new_rows) + ("\n" if end else ""), encoding="utf-8")
    report[f"S15 {folder.name} per_protein"] = {"with d_mean_plant": found, "rows": len(rows) - 1}

pc = pd.read_csv(S15[0] / "per_protein_comparison.tsv", sep="\t", low_memory=False)
stats = []
for name, col in (("all proteins", None), ("Arabidopsis/rice/maize/soybean", "crop")):
    df = pc if col is None else pc[pc[col].notna()]
    for sel in ("sel3800", "sel690"):
        for ref in ("pdb3825_full", "pdb3825_construct", "random_matched", "random_pool"):
            x, y = df.loc[df.set.eq(sel), COL].dropna(), df.loc[df.set.eq(ref), COL].dropna()
            if len(x) < 10 or len(y) < 10:
                continue
            u = mannwhitneyu(x, y, alternative="two-sided")
            stats.append({"subset": name, "measure": LABEL, "AI set": sel, "reference": ref, "n_AI": len(x), "n_ref": len(y),
                          "AI median": x.median(), "ref median": y.median(),
                          "Cliff's delta (AI vs ref)": 2 * u.statistic / (len(x) * len(y)) - 1, "p": u.pvalue})
for folder in S15:
    p = folder / "comparison_statistics.tsv"
    backup(p)
    lines = p.read_text(encoding="utf-8").split("\n")
    end = lines[-1] == ""
    rows = [l.split("\t") for l in (lines[:-1] if end else lines)]
    head = rows[0]
    im, isub, isel, iref = head.index("measure"), head.index("subset"), head.index("AI set"), head.index("reference")
    rows = [r for r in rows if r[im] != LABEL]
    fmt = {}
    for s in stats:
        fmt[(s["subset"], s["AI set"], s["reference"])] = [
            str(s.get(h, "")) if isinstance(s.get(h, ""), (str, int, np.integer)) else repr(float(s[h])) for h in head]
    out_rows, placed = [head], 0
    for r in rows[1:]:
        out_rows.append(r)
        key = (r[isub], r[isel], r[iref])
        if r[im] == AFTER_MEASURE and key in fmt:
            out_rows.append(fmt.pop(key))
            placed += 1
    out_rows.extend(fmt.values())
    p.write_text("\n".join("\t".join(r) for r in out_rows) + ("\n" if end else ""), encoding="utf-8")
    report[f"S15 {folder.name} stats"] = {"added": len(stats), "placed after nearest-plant rows": placed}

print(json.dumps(report, indent=2, default=str))
show = pd.DataFrame(stats)
print(show[["subset", "AI set", "reference", "AI median", "ref median", "Cliff's delta (AI vs ref)", "p"]].to_string(
    index=False, float_format=lambda v: f"{v:.3g}"))
