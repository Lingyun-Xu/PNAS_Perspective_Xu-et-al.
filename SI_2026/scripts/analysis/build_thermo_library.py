#!/usr/bin/env python3
"""Thermostability template pool for Fig. 4 on AFDB v6.
patch     : FASTA records of the 2025 library (FireProtDB wild types, mega-scale Dataset 2/3, Meltome Atlas) without an
            unambiguous February 2025 embedding -> inputs/thermo_patch.fasta (embedded with the v6 extractor). The
            seven ids 'nan', 'nan_2' ... 'nan_7' each name one FireProtDB and one Meltome record with different
            sequences; one file per label survived in 2025, so both records of each pair are re-embedded.
assemble  : inputs/thermostability_library_2025.pt = the February 2025 embeddings (without the 'nan' labels) + the
            patch, one row per FASTA record (781,320), ids '<source>|<id>' and the source of each row.
python build_thermo_library.py patch|assemble"""
import re
import sys
from pathlib import Path

ROOT = Path("/path/to/cluster/redo_PDB/ESM_embedding/template_v6")
sys.path.insert(0, str(ROOT))
import numpy as np
import pandas as pd
import torch

import compute_distances as cd

W = ROOT / "campaign/fig4_v6"
TAG = {"mega-scale Dataset 2/3": "mega23", "FireProtDB": "fireprot", "Meltome Atlas": "meltome"}
NAN = re.compile(r"^nan(_\d+)?$")
norm = lambda x: str(x).lstrip(">").split()[0]

rec = pd.read_csv(W / "inputs/thermostability_fasta_records.tsv.gz", sep="\t", keep_default_na=False)
rec["nid"] = rec.id.map(norm)
rec["key"] = rec.source.map(TAG) + "|" + rec.nid
assert rec.key.is_unique, rec.key[rec.key.duplicated()].tolist()[:5]
old = torch.load(W / "inputs/thermo_old_emb.pt", map_location="cpu", weights_only=False)
oid = [norm(i) for i in old["ids"]]
keep = [k for k, i in enumerate(oid) if not NAN.match(i)]
by_nid = rec[~rec.nid.map(lambda s: bool(NAN.match(s)))].set_index("nid")
assert by_nid.index.is_unique
kept_ids = [oid[k] for k in keep]
unknown = [i for i in kept_ids if i not in by_nid.index]
covered = set(kept_ids) - set(unknown)
missing = rec[~rec.nid.isin(covered)]
print(f"2025 embeddings {len(oid):,}; 'nan' labels dropped {len(oid) - len(keep)}; labels not in the FASTA records "
      f"{len(unknown)} {unknown[:5]}; FASTA records without an embedding {len(missing)}:")
print(missing[["source", "nid"]].assign(length=missing.sequence.str.len()).to_string(index=False))

if sys.argv[1] == "patch":
    with open(W / "inputs/thermo_patch.fasta", "w") as fh:
        for k, s in zip(missing.key, missing.sequence):
            fh.write(f">{k}\n{s}\n")
    print("written inputs/thermo_patch.fasta", len(missing))
else:
    pemb, pidx = cd.load_emb([str(W / "emb_thermo_patch")])
    assert set(pidx) == set(missing.key), (len(pidx), len(missing))
    emb = np.concatenate([old["emb"].float().numpy()[keep], pemb[[pidx[k] for k in missing.key]]])
    src_of = dict(zip(by_nid.index, by_nid.source))
    keys = [src_of[i] for i in kept_ids]
    ids = [TAG[s] + "|" + i for s, i in zip(keys, kept_ids)] + list(missing.key)
    src = keys + list(missing.source)
    assert len(ids) == len(rec) == len(set(ids)), (len(ids), len(rec))
    seqs = dict(zip(rec.key, rec.sequence))
    print(f"library rows {len(ids):,}; unique sequences {len({seqs[i] for i in ids}):,}; by source "
          f"{pd.Series(src).value_counts().to_dict()}")
    torch.save({"ids": ids, "source": src, "emb": torch.from_numpy(emb.astype(np.float32))},
               W / "inputs/thermostability_library_2025.pt")
    print("written inputs/thermostability_library_2025.pt")
