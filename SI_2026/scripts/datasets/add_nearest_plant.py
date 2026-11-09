#!/usr/bin/env python3
"""Cutoff-free embedding neighbour measure for Section S7.5 (author request, 6 Nov 2026).

Adds d_nearest_plant, the Euclidean distance between the ESM-2 (650M) mean embedding of a protein and that of the
nearest other plant sequence of the AFDB v6 set (182,286 unique sequences; identical sequences not counted), from
neighbour_distances.py (cluster; sets_nearest_plant.tsv), to per_protein_comparison.tsv, and its Mann-Whitney tests
(Cliff's delta, as compare_sets.py) to comparison_statistics.tsv. Both copies are edited (SI_2026/pdb3825_vs_selections
and Dataset S15) at the text level, so that every existing value keeps its exact text. The former neighbour count
(n_within_1.5) is kept. Re-running is safe: an existing d_nearest_plant column / rows are replaced.

Usage (needs scipy): python add_nearest_plant.py <sets_nearest_plant.tsv>
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

B = Path("/path/to/project") if sys.platform.startswith("linux") else Path(r"/path/to/project")
COPIES = [B / "SI_2026/pdb3825_vs_selections", B / "SI_2026/datasets/Dataset_S15_compared_protein_sets"]
COL = "d_nearest_plant"
LABEL = "Distance to the nearest other plant sequence"
AFTER_COL = "n_within_2.0"
AFTER_MEASURE = "Plant proteins within ESM-2 distance 1.5"
SUBSETS = [("all proteins", None), ("Arabidopsis/rice/maize/soybean", "crop")]
REFS = ("pdb3825_full", "pdb3825_construct", "random_matched", "random_pool")


def cliffs_delta(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    u = mannwhitneyu(a, b, alternative="two-sided")
    return 2 * u.statistic / (len(a) * len(b)) - 1, u.pvalue


def add_column(path, values):
    lines = path.read_text(encoding="utf-8").split("\n")
    end = lines[-1] == ""
    rows = [l.split("\t") for l in (lines[:-1] if end else lines)]
    head = rows[0]
    if COL in head:                                              # re-run: drop the old column first
        k = head.index(COL)
        rows = [r[:k] + r[k + 1:] for r in rows]
        head = rows[0]
    i_set, i_key, i_after = head.index("set"), head.index("key"), head.index(AFTER_COL)
    out, found = [head[:i_after + 1] + [COL] + head[i_after + 1:]], 0
    for r in rows[1:]:
        v = values.get((r[i_set], r[i_key]))
        found += v is not None
        out.append(r[:i_after + 1] + ["" if v is None else repr(float(v))] + r[i_after + 1:])
    path.write_text("\n".join("\t".join(r) for r in out) + ("\n" if end else ""), encoding="utf-8")
    return found, len(rows) - 1


def stats_rows(pc):
    rows = []
    for name, col in SUBSETS:
        df = pc if col is None else pc[pc[col].notna()]
        for sel in ("sel3800", "sel690"):
            for ref in REFS:
                x = df.loc[df.set.eq(sel), COL].dropna()
                y = df.loc[df.set.eq(ref), COL].dropna()
                if len(x) < 10 or len(y) < 10:
                    continue
                d, p = cliffs_delta(x, y)
                rows.append({"subset": name, "measure": LABEL, "AI set": sel, "reference": ref, "n_AI": len(x),
                             "n_ref": len(y), "AI median": x.median(), "ref median": y.median(),
                             "Cliff's delta (AI vs ref)": d, "p": p})
    return rows


def add_stats(path, new_rows):
    lines = path.read_text(encoding="utf-8").split("\n")
    end = lines[-1] == ""
    rows = [l.split("\t") for l in (lines[:-1] if end else lines)]
    head = rows[0]
    im, isub, isel, iref = head.index("measure"), head.index("subset"), head.index("AI set"), head.index("reference")
    rows = [r for r in rows if r[im] != LABEL]                  # re-run
    fmt = {}
    for nr in new_rows:
        cells = []
        for h in head:
            v = nr.get(h, "")
            cells.append(str(v) if isinstance(v, (int, np.integer)) or isinstance(v, str) else repr(float(v)))
        fmt[(nr["subset"], nr["AI set"], nr["reference"])] = cells
    out, placed = [head], 0
    for r in rows[1:]:
        out.append(r)
        key = (r[isub], r[isel], r[iref])
        if r[im] == AFTER_MEASURE and key in fmt:
            out.append(fmt.pop(key))
            placed += 1
    out.extend(fmt.values())                                     # any without an anchor row go at the end
    path.write_text("\n".join("\t".join(r) for r in out) + ("\n" if end else ""), encoding="utf-8")
    return placed, len(new_rows)


def main():
    nbr = pd.read_csv(sys.argv[1], sep="\t", dtype={"query_id": str})
    values = {(s, q): d for s, q, d in zip(nbr.set, nbr.query_id, nbr.d_nearest_plant)}
    for folder in COPIES:
        found, n = add_column(folder / "per_protein_comparison.tsv", values)
        print(f"{folder.name}: per_protein_comparison.tsv, {COL} for {found:,} of {n:,} rows", flush=True)
    pc = pd.read_csv(COPIES[0] / "per_protein_comparison.tsv", sep="\t", low_memory=False)
    print(pc.groupby("set")[[COL, "n_within_1.5"]].count().to_string())
    new = stats_rows(pc)
    for folder in COPIES:
        placed, n = add_stats(folder / "comparison_statistics.tsv", new)
        print(f"{folder.name}: comparison_statistics.tsv, {n} rows added ({placed} after the neighbour-count rows)")
    show = pd.DataFrame(new)
    pd.set_option("display.width", 220)
    print(show[["subset", "AI set", "reference", "n_AI", "n_ref", "AI median", "ref median",
                "Cliff's delta (AI vs ref)", "p"]].to_string(index=False, float_format=lambda v: f"{v:.3g}"))


if __name__ == "__main__":
    main()
