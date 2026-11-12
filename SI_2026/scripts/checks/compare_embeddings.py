#!/usr/bin/env python3
"""Are the reference-set embeddings of the first version of Fig. 4 the same vectors as the v6 extractor gives?
  MaveDB           all 1,170: February 2025 extract.py embeddings (inputs/mavedb_old_emb.pt) vs emb_mavedb
  thermostability  400 random sequences: February 2025 embeddings (inputs/thermostability_library.pt) vs emb_thermo_check
  mdCATH           Diego's December 2024 tensor (inputs/mdCATH_embeddings_2024.pt, no ids) vs emb_mdcath: nearest row
Distances are Euclidean, as everywhere in the selection; for scale, the median nearest-template distance of the
candidate pool is about 1.6."""
import json
import sys
from pathlib import Path

ROOT = Path("/path/to/cluster/redo_PDB/ESM_embedding/template_v6")
sys.path.insert(0, str(ROOT))
import numpy as np
import torch

import compute_distances as cd

W = ROOT / "campaign/fig4_v6"
res = {}


def norm_id(x):
    x = str(x).lstrip(">").split()[0]
    head, _, rest = x.partition("|")
    return rest if head in ("mega23", "fireprot", "meltome") else x


def paired(tag, old_file, new_dir):
    old = torch.load(W / old_file, map_location="cpu", weights_only=False)
    oid = {norm_id(i): k for k, i in enumerate(old["ids"])}
    oemb = old["emb"].float().numpy()
    nemb, nidx = cd.load_emb([str(W / new_dir)])
    nidx = {norm_id(i): k for i, k in nidx.items()}
    common = [i for i in nidx if i in oid and not i.startswith("nan")]       # 'nan' labels are ambiguous
    d = np.linalg.norm(oemb[[oid[i] for i in common]] - nemb[[nidx[i] for i in common]], axis=1)
    res[tag] = {"new": len(nidx), "matched_by_id": len(common), "median_distance": float(np.median(d)),
                "p99_distance": float(np.percentile(d, 99)), "max_distance": float(d.max()),
                "n_distance_gt_0.01": int((d > 0.01).sum())}
    worst = np.argsort(-d)[:5]
    res[tag]["worst"] = [(common[j], float(d[j])) for j in worst]
    print(tag, json.dumps(res[tag]), flush=True)


paired("mavedb", "inputs/mavedb_old_emb.pt", "emb_mavedb")
if (W / "inputs/thermo_old_emb.pt").exists():
    paired("thermostability_sample", "inputs/thermo_old_emb.pt", "emb_thermo_check")

old = torch.load(W / "inputs/mdCATH_embeddings_2024.pt", map_location="cpu", weights_only=False).float().numpy()
MD_DIR, MD_FA = (("emb_mdcath_atom", "inputs/mdcath_5398_s20atom.fasta") if (W / "emb_mdcath_atom").exists()
                 else ("emb_mdcath", "inputs/mdcath_5398_domains.fasta"))
nemb, nidx = cd.load_emb([str(W / MD_DIR)])
doms = [l.split()[0][1:] for l in open(W / MD_FA) if l.startswith(">")]
seqlen = {}
name = None
for line in open(W / MD_FA):
    if line.startswith(">"):
        name = line[1:].strip()
    else:
        seqlen[name] = len(line.strip())
pos = {d: k for k, d in enumerate(doms)}
on, nn = (old * old).sum(1), (nemb * nemb).sum(1)
d2 = nn[:, None] + on[None, :] - 2 * nemb @ old.T
d = np.sqrt(np.maximum(d2, 0))
best, bestd = d.argmin(1), d.min(1)
names = sorted(nidx, key=nidx.get)
same_pos = np.array([best[nidx[n]] == pos[n] for n in names])
res["mdcath"] = {"v6_domains": len(names), "rows_2024": int(len(old)),
                 "median_nearest_2024_row": float(np.median(bestd)), "p90": float(np.percentile(bestd, 90)),
                 "nearest_row_is_same_list_position": int(same_pos.sum()),
                 "median_distance_same_position": float(np.median([d[nidx[n], pos[n]] for n in names])),
                 "nearest_lt_0.01": int((bestd < 0.01).sum())}
print("mdcath", json.dumps(res["mdcath"]), flush=True)
(W / "results").mkdir(exist_ok=True)
(W / "results/embedding_checks.json").write_text(json.dumps(res, indent=1))
print("CHECKS_DONE", flush=True)
