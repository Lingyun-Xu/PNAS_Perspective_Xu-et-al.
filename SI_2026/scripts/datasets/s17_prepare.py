"""Dataset S17 (source data of SI figures): convert the Fig. S12 round arrays to one TSV, and print the shape and
columns of every table in the folder (for the README)."""
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

D = Path("/path/to/project/SI_2026/datasets/Dataset_S17_source_data_SI_figures")

raw = D / "FigS12" / "raw"
if raw.exists():
    parts = []
    for r in ("000", "037"):
        mn, av, pk = (np.load(raw / f"round_{r}_{k}.npy") for k in ("min_dist", "avg_dist", "picked"))
        sel = np.zeros(len(mn), dtype=bool)
        sel[pk] = True
        parts.append(pd.DataFrame({"round": int(r), "pool_position": np.arange(len(mn)), "d_nearest_library": mn,
                                   "d_average_homolog": av, "selected_this_round": sel}))
    out = pd.concat(parts, ignore_index=True)
    out.to_csv(D / "FigS12" / "selection_rounds_0_and_37.tsv.gz", sep="\t", index=False, float_format="%.6g")
    shutil.rmtree(raw)
    print("FigS12 rounds:", out.groupby("round").size().to_dict(), "selected:", out.groupby("round").selected_this_round.sum().to_dict())

for f in sorted(D.rglob("*")):
    if f.is_dir():
        continue
    name = str(f.relative_to(D))
    try:
        if f.suffix == ".parquet":
            t = pd.read_parquet(f)
        elif f.name.endswith((".tsv", ".tsv.gz")):
            t = pd.read_csv(f, sep="\t", low_memory=False)
        elif f.suffix == ".csv":
            t = pd.read_csv(f, low_memory=False)
        elif f.suffix == ".json":
            print(f"{name}: {json.loads(f.read_text())}")
            continue
        else:
            print(f"{name}: {f.stat().st_size:,} bytes")
            continue
        print(f"{name}: {len(t):,} rows | {list(t.columns)}")
    except Exception as exc:
        print(f"{name}: unreadable ({exc})")
