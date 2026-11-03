"""Dataset S16: source data of the main-text figures (Figs. 1-3; Fig. 4a reuses Fig. 3, Fig. 4b has no data).
Inputs were located in the 2024-2026 notebooks (workstation /path/to/home), their copies (MD5-identical) and the
analysis server (moved from the cluster; pulled to a local working folder (s16_raw) with MD5 checks). Tensors (.pt) are converted to TSV; the
plotted counts of Fig. 2 are recomputed with the notebook logic and checked against the printed or labelled values."""
import ast
import gzip
import json
import re
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import torch

B = Path("/path/to/project")
H = Path("/path/to/home")
RAW = Path("/path/to/scratch/s16_raw")
OUT = B / "SI_2026/datasets/Dataset_S16_source_data_main_text_figures"
OUT.mkdir(parents=True, exist_ok=True)


def sub(name):
    p = OUT / name
    p.mkdir(parents=True, exist_ok=True)
    return p


def copy(src, folder, gz=False):
    src = Path(src)
    dst = sub(folder) / (src.name + (".gz" if gz else ""))
    if gz:
        with open(src, "rb") as fi, gzip.open(dst, "wb", compresslevel=6) as fo:
            shutil.copyfileobj(fi, fo)
    else:
        shutil.copy2(src, dst)
    print(f"  {folder}/{dst.name} ({dst.stat().st_size:,} B)")


def tensor(path):
    return torch.load(path, map_location="cpu", weights_only=False).detach().cpu().numpy().astype(np.float64).ravel()


# ---------------------------------------------------------------- Fig. 1a: global pLDDT by group (AFDB v4)
print("Fig. 1a")
d = pd.read_parquet(B / "SI_2026/datasets/Dataset_S4_AFDB_v4_48_proteomes/AF2_48_model_oganisms.parquet",
                    columns=["pdbx_PDB_id_code", "db_accession", "organism_scientific", "plddt"])
gm = pd.read_csv(B / "alphafoldDB/AF2_unique_species_find_what_they_are_shortenedagain_082524.csv", encoding="ISO-8859-1")
gm["group"] = gm["superkindom"].astype(str).str.replace(r"[^A-Za-z]", "", regex=True)
grp = dict(zip(gm.cleaned_species.str.strip(), gm.group))
d["group"] = d.organism_scientific.str.strip().map(grp)
d["global_plddt"] = pd.to_numeric(d.plddt, errors="coerce")
print("  unmapped organisms:", int(d.group.isna().sum()), "| per group:", d.group.value_counts().to_dict())
d.rename(columns={"pdbx_PDB_id_code": "afdb_entry_id"})[["afdb_entry_id", "db_accession", "organism_scientific", "group", "global_plddt"]] \
    .to_csv(sub("Fig1a") / "Fig1a_global_plddt_by_group.tsv.gz", sep="\t", index=False)
gm[["cleaned_species", "species_NIH", "group"]].to_csv(sub("Fig1a") / "group_of_each_proteome.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- Fig. 1b and 1h-k: MSA depth (Neff/length) and pLDDT
print("Fig. 1b, 1h-k")
copy(B / "uniclust30_MSA/MSA_depth_with_plddt_34921.csv", "Fig1b_Fig1h-k")

# ---------------------------------------------------------------- Fig. 1c: distance to nearest PDB (May 2024) vs pLDDT
print("Fig. 1c")
for name in ("AF2plant_exclude_crystalizedentries", "ECOLI", "HUMAN", "YEAST"):
    copy(B / f"2021_now_published_plant_PDB/redo_PDB/recal_distance/{name}_to_redoPDB_distance.parquet", "Fig1c")
for sp in ("ECOLI", "HUMAN", "YEAST"):
    copy(B / f"AF2_clustering/Embeddings_distance/AF2_{sp}.csv", "Fig1c")

# ---------------------------------------------------------------- Fig. 1d-g: nearest-PDB distance, February 2025 tensors
print("Fig. 1d-g")
parts = []
for sp in ("Plants", "Ecoli", "Yeast", "Human"):
    v = tensor(H / f"AF2_clustering/min_pairwise_dist_{sp}_021125.pt")
    parts.append(pd.DataFrame({"species": sp, "d_nearest_pdb": v}))
    print(f"  {sp}: {len(v):,} values, median {np.median(v):.3f}")
pd.concat(parts).to_csv(sub("Fig1d-g") / "nearest_PDB_distance_Feb2025.tsv.gz", sep="\t", index=False, float_format="%.6g")

# ---------------------------------------------------------------- Fig. 2
print("Fig. 2")
counts = []
# a: UniProtKB entries, typed into Fig2_Pie_chart.ipynb cell 116 (snapshot 2 July 2026)
a = {"Bacteria": 74117227, "Fungi": 17581490, "Viridiplantae": 17344395, "Other Vertebrata": 14959291, "Other": 12374573,
     "Mammalia": 5844956, "Insecta": 4715693, "Archaea": 1767834, "Viruses": 1104680}
counts += [("2a", "UniProtKB entries, 2 July 2026", g, n) for g, n in a.items()]
# b: PDB unique protein sequences, July 2026
raw_b = json.load(open(RAW / "fig2b/unique_protein_sequences.json"))
if isinstance(raw_b, dict):   # either {"records": [...]} or {sequence: record}
    raw_b = raw_b.get("records", list(raw_b.values()))
recs = pd.DataFrame(raw_b)
print("  2b columns:", list(recs.columns))
copy(RAW / "fig2b/unique_protein_sequences.json", "Fig2b", gz=True)
copy(B / "Fig2/20260705/pdb_superkingdom_still_unmapped_hand_curated.csv", "Fig2b")
copy(RAW / "fig2b/nih_species_cache.json", "Fig2b")
copy(RAW / "fig2b/fetch_pdb_unique_sequences.py", "Fig2b")
has_org = recs.source_organism.notna() & recs.source_organism.astype(str).str.strip().ne("")
auto = recs.loc[has_org, "superkingdom"].value_counts(dropna=False)
cur = pd.read_csv(B / "Fig2/20260705/pdb_superkingdom_still_unmapped_hand_curated.csv")
cmap = dict(zip(cur.species, cur.proposed_superkingdom.str.replace("_", " ")))
b = recs.loc[has_org].copy()
b["group"] = b.superkingdom.where(b.superkingdom.notna(), b.species.map(cmap))
print(f"  2b: {len(recs):,} sequences; {int(has_org.sum()):,} with an organism; automatic label {int(recs.loc[has_org, 'superkingdom'].notna().sum()):,};"
      f" + hand-curated file {int((recs.loc[has_org, 'superkingdom'].isna() & b.group.notna()).sum()):,}; still unassigned {int(b.group.isna().sum()):,}")
for g, n in b.group.value_counts(dropna=False).items():
    counts.append(("2b", "PDB unique sequences with organism, automatic + hand-curated file (unassigned rows kept as NA)", str(g), int(n)))
# c: thermostability measurements (FireProtDB + mega-scale + Meltome Atlas)
fm = pd.read_csv(H / "GEO_dataset/firemega_complete_species.csv", low_memory=False)
me = pd.read_csv(H / "GEO_dataset/meltome_altas_with_organisms_species_superkindoms.csv", low_memory=False)
copy(H / "GEO_dataset/firemega_complete_species.csv", "Fig2c", gz=True)
copy(H / "GEO_dataset/meltome_altas_with_organisms_species_superkindoms.csv", "Fig2c", gz=True)
c = pd.concat([fm.superkingdom, me.superkingdom]).replace({"Nematoda": "Others"}).value_counts()
print(f"  2c: {len(fm):,} + {len(me):,} = {len(fm) + len(me):,} rows; " + "; ".join(f"{k} {v:,}" for k, v in c.items()))
counts += [("2c", "thermostability measurements (rows)", g, int(n)) for g, n in c.items()]
# d: MaveDB score sets, groups from the lists of Fig2_Pie_chart.ipynb cell 58
nb = json.load(open(H / "GEO_dataset/Fig2_Pie_chart.ipynb"))
lists = {n.targets[0].id: set(ast.literal_eval(n.value)) for n in ast.parse("".join(nb["cells"][58]["source"])).body
         if isinstance(n, ast.Assign) and isinstance(n.value, ast.List)}
order = ["Bacteria", "Fungi", "Insecta", "Other_Mammalia", "Homo", "Nematoda", "Viridiplantae", "Other_Vertebrata",
         "Archaea", "Viruses", "Others", "synthetic"]
lab = {"Other_Mammalia": "Other Mammalia", "Other_Vertebrata": "Other Vertebrata", "synthetic": "Synthetic"}
mave = pd.read_csv(H / "GEO_dataset/mavedb-organisms_species_and_numbers.csv")
mave["group_listed"] = mave.Species.map(lambda s: next((lab.get(k, k) for k in order if s in lists[k]), None))
mave["pie_group"] = mave.group_listed.replace({"Insecta": "Others", "Nematoda": "Others"})
mave.to_csv(sub("Fig2d") / "mavedb_organisms_score_sets_groups.tsv", sep="\t", index=False)
dd = mave.groupby("pie_group")["Associated Score Sets"].sum().sort_values(ascending=False)
print("  2d: " + "; ".join(f"{k} {v} ({100 * v / dd.sum():.2f}%)" for k, v in dd.items()))
counts += [("2d", "MaveDB score sets (August 2024 export)", g, int(n)) for g, n in dd.items()]
# e: mdCATH domains' PDB entities, groups re-merged with the later species table (as in the plotted cell)
md = pd.read_parquet(H / "MDcath/MD_PDB_entries_with_superkingdom.parquet")
sp = pd.read_csv(H / "MDcath/MD_unique_species_with_superkingdom.csv")
copy(H / "MDcath/MD_PDB_entries_with_superkingdom.parquet", "Fig2e")
copy(H / "MDcath/MD_unique_species_with_superkingdom.csv", "Fig2e")
e = md.drop(columns=["superkingdom"]).merge(sp[["cleaned_species", "superkingdom"]], left_on="organism_name",
                                            right_on="cleaned_species", how="left").superkingdom.value_counts()
print(f"  2e: {int(e.sum()):,} rows; " + "; ".join(f"{k} {v:,}" for k, v in e.items()))
counts += [("2e", "mdCATH PDB entities (rows)", g, int(n)) for g, n in e.items()]
pd.DataFrame(counts, columns=["panel", "unit", "group", "count"]).to_csv(sub("Fig2") / "Fig2_plotted_counts.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- Fig. 3 (and Fig. 4a): selection rounds
print("Fig. 3")
runs = [("a-e_PDB_without_2021-2025_plant_structures", RAW / "fig3_pdb_excl", (0, 17)),
        ("b_all_PDB_2025-01-24", RAW / "fig3_pdb_all", (0,)),
        ("f-j_thermostability", RAW / "fig3_thermo", (0, 1, 20)),
        ("k-o_MaveDB", RAW / "fig3_mavedb", (0, 1, 40)),
        ("p-t_mdCATH", H / "MDcath/Diego_code_and_distance_data/AF2plants_vs_mdCATH", (0, 1, 20))]
for label, root, rounds in runs:
    frames = []
    for r in rounds:
        mn = tensor(root / f"round_{r}/min_pairwise_dist.pt")
        av = tensor(root / f"round_{r}/average_internal_dist.pt")
        assert len(mn) == len(av), (label, r)
        frames.append(pd.DataFrame({"round": r, "pool_position": np.arange(len(mn)), "d_nearest_reference": mn,
                                    "d_average_homolog": av}))
        print(f"  {label} round {r}: {len(mn):,} plant proteins")
    pd.concat(frames).to_csv(sub("Fig3") / f"Fig3{label}.tsv.gz", sep="\t", index=False, float_format="%.6g")
# the analysis server run of panels a and c equals the copies used earlier (round 17 = 'round_18' file; round 0 average)
w = B / "2026/AF2_species_vs_250124PDB_exclude_4yearplants250115"
for vm_file, win_file in ((RAW / "fig3_pdb_excl/round_17/min_pairwise_dist.pt", w / "min_pairwise_dist_round_18_PlantvsPDB.pt"),
                          (RAW / "fig3_pdb_excl/round_0/average_internal_dist.pt", w / "average_internal_dist_plant_round_0.pt")):
    if win_file.exists():
        x, y = tensor(vm_file), tensor(win_file)
        print(f"  check {win_file.name}: same length {len(x) == len(y)}, identical {len(x) == len(y) and np.allclose(x, y)}")
print("done")
