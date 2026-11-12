#!/usr/bin/env python3
"""Protein sequences of all 1,841,285 rows of the mega-scale Dataset 1 (Tsuboyama et al. 2023), in the form the 2025
thermostability library used for Dataset 2/3 (aa_seq: the designed or natural domain without the SAGG linker padding).
Each row's dna_seq is translated (the padded construct, = aa_seq_full in Dataset 2/3) and its padding removed:
  1. construct also in Dataset 2/3: its aa_seq (exact);
  2. else, construct family (name before '.pdb') in Dataset 2/3 and the construct starts and ends with that family's
     most common left and right pads: those pads removed;
  3. else the generic rule: the longest prefix and suffix that are prefixes of SAGGSAGG... removed (on Dataset 2/3 this
     reproduces aa_seq for 82% of constructs; most misses differ by one terminal residue).
Writes <out>/mega_dataset1_rows.tsv.gz (row, name, sequence, rule) and <out>/mega_dataset1_new.fasta (the unique
sequences not already in the 2025 library FASTA records), and prints the counts.
python build_mega_d1.py <thermostability_fasta_records.tsv.gz> <out dir>"""
import sys
from collections import Counter
from pathlib import Path

import pandas as pd

D = "/path/to/project/mega_scale_stability_dataset/"
REC, OUT = sys.argv[1], Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)
CODON = {}
bases = "TCAG"
aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
k = 0
for a in bases:
    for b in bases:
        for c in bases:
            CODON[a + b + c] = aas[k]
            k += 1
PAD = "SAGG" * 12


def tr(dna):
    dna = str(dna).upper()
    return "".join(CODON.get(dna[i:i + 3], "X") for i in range(0, len(dna) - len(dna) % 3, 3))


def family(name):
    return str(name).split(".pdb")[0]


def generic(f):
    i = 0
    while i < len(f) and i < len(PAD) and f[i] == PAD[i]:
        i += 1
    core = f[i:]
    kb = max([k for k in range(1, min(len(core), len(PAD)) + 1) if core[-k:] == PAD[:k]], default=0)
    return core[:len(core) - kb] if kb else core


d23 = pd.read_csv(D + "Tsuboyama2023_Dataset2_Dataset3_20230416.csv", usecols=["name", "aa_seq_full", "aa_seq"],
                  low_memory=False)
full2seq = dict(zip(d23.aa_seq_full, d23.aa_seq))
pads = {}
for n, f, s in zip(d23.name, d23.aa_seq_full.astype(str), d23.aa_seq.astype(str)):
    i = f.find(s)
    if i >= 0:
        pads.setdefault(family(n), Counter())[(f[:i], f[i + len(s):])] += 1
fam_pad = {fam: c.most_common(1)[0][0] for fam, c in pads.items()}

d1 = pd.read_csv(D + "Tsuboyama2023_Dataset1_20230416.csv", usecols=["name", "dna_seq"], low_memory=False)
assert len(d1) == 1841285, len(d1)
full = d1.dna_seq.map(tr)
seq, rule = [], []
for n, f in zip(d1.name, full):
    s = full2seq.get(f)
    if s is not None:
        seq.append(s); rule.append("in Dataset 2/3")
        continue
    lp, rp = fam_pad.get(family(n), (None, None))
    if lp is not None and f.startswith(lp) and f.endswith(rp) and len(f) > len(lp) + len(rp):
        seq.append(f[len(lp):len(f) - len(rp)] if rp else f[len(lp):]); rule.append("family pads")
        continue
    seq.append(generic(f)); rule.append("generic rule")
rows = pd.DataFrame({"row": range(len(d1)), "name": d1.name, "sequence": seq, "rule": rule})
rows["has_stop_or_X"] = rows.sequence.str.contains(r"[*X]")
rows.to_csv(OUT / "mega_dataset1_rows.tsv.gz", sep="\t", index=False, compression="gzip")
print("rows by rule:", rows.rule.value_counts().to_dict())
print("rows with a stop codon or untranslatable codon in the domain:", int(rows.has_stop_or_X.sum()))
ok = rows[~rows.has_stop_or_X & rows.sequence.str.len().gt(0)]
lib = pd.read_csv(REC, sep="\t")
have = set(lib.sequence)
uniq = ok.drop_duplicates("sequence")
new = uniq[~uniq.sequence.isin(have)]
print(f"unique Dataset 1 domain sequences {uniq.sequence.nunique():,}; already in the 2025 library {int(uniq.sequence.isin(have).sum()):,}; "
      f"new {len(new):,} (by rule: {new.rule.value_counts().to_dict()})")
print("Dataset 2/3 library sequences also in Dataset 1:",
      int(lib[lib.source.eq('mega-scale Dataset 2/3')].sequence.isin(set(uniq.sequence)).sum()), "of 740,697")
print("new sequence lengths: median", new.sequence.str.len().median(), "max", new.sequence.str.len().max(),
      "total residues", int(new.sequence.str.len().sum()))
with open(OUT / "mega_dataset1_new.fasta", "w") as fh:
    for i, s in enumerate(new.sequence):
        fh.write(f">megaD1_{i:07d}\n{s}\n")
new.assign(id=[f"megaD1_{i:07d}" for i in range(len(new))])[["id", "name", "rule", "sequence"]].to_csv(
    OUT / "mega_dataset1_new_ids.tsv.gz", sep="\t", index=False, compression="gzip")
print("written", OUT)
