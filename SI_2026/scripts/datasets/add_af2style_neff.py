#!/usr/bin/env python3
"""Add the AF2-style Neff of neff_v6/fetch_af2style_neff.py to Dataset S10 (merged_neff.parquet) as the column
AF2-style_Neff, placed after neff_per_len. Each protein takes the value of its assigned alignment (rep_id), as for neff.
Stops unless every alignment reproduced the Neff of the table (the fetch script recomputes it from the same rows) and
every existing column stays identical. Safe to rerun: an existing AF2-style_Neff column is replaced.
Usage: python add_af2style_neff.py <af2style_neff_by_cluster.tsv.gz> <merged_neff.parquet>"""
import os
import sys

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

COL = "AF2-style_Neff"
src, path = sys.argv[1], sys.argv[2]
a = pd.read_csv(src, sep="\t")
a = a[a.status.isin(["ok", "ok_cpu"])]
assert not a.rep_id.duplicated().any(), "duplicate alignments"
pf = pq.ParquetFile(path)
codec = pf.metadata.row_group(0).column(0).compression
old = pf.read().to_pandas()
if COL in old:
    old = old.drop(columns=COL)
print(f"proteins {len(old):,}; with neff {old.neff.notna().sum():,}; alignments in the new table {len(a):,}; codec {codec}")
value = old.rep_id.astype("Int64").map(a.set_index("rep_id").af2style_neff).astype("float64")
check = old.rep_id.astype("Int64").map(a.set_index("rep_id").neff_recomputed).astype("float64")
has = old.neff.notna()
assert value[has].notna().all(), "a protein with Neff has no AF2-style value"
assert value[~has].isna().all(), "a protein without Neff has an AF2-style value"
diff = (check[has] - old.neff[has]).abs().to_numpy()
assert diff.max() <= 0.001, f"Neff not reproduced: {diff.max()}"
print(f"Neff reproduced for all {has.sum():,} proteins; largest difference {diff.max():.4f}")
out = old.copy()
out.insert(out.columns.get_loc("neff_per_len") + 1, COL, value)
tmp = path + ".tmp"
pq.write_table(pa.Table.from_pandas(out, preserve_index=False), tmp, compression=codec)
back = pd.read_parquet(tmp)
for c in old.columns:
    assert back[c].equals(old[c]) or (back[c].fillna(-1).astype(str) == old[c].fillna(-1).astype(str)).all(), f"column {c} changed"
assert back[COL].notna().sum() == has.sum() and np.allclose(back[COL].dropna().to_numpy(), value.dropna().to_numpy())
print("existing columns unchanged; size", round(os.path.getsize(path) / 2 ** 20, 1), "->", round(os.path.getsize(tmp) / 2 ** 20, 1), "MiB")
os.replace(tmp, path)
print("columns:", list(back.columns))
