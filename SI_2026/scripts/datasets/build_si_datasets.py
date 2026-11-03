#!/usr/bin/env python3
"""Build SI Datasets S1-S3 (Excel) from existing result tables; print summary numbers for the SI text.

S1: original AI-guided selection (Feb 2025; AFDB v4 plants vs PDB of 24 Jan 2025 without 2021-2025 plant
    structures; 18 rounds x 50; the first 850 = main-text Fig. 3c), replayed from the saved per-round indices.
S2: sel3800 (AFDB v6; 38 rounds x 100), with annotation, Pfam, AIUPred disorder and GenEra age.
S3: sel690 (AFDB v6; 69 rounds x 10), same columns.
"""
from pathlib import Path

import pandas as pd

ROOT = Path(r"/path/to/project")
CAMP = ROOT / "template_esm" / "campaign" / "figures"
OUT = ROOT / "SI_2026" / "datasets"


def readme(ws_rows):
    return pd.DataFrame({"README": ws_rows})


def write_xlsx(path, data, notes):
    with pd.ExcelWriter(path, engine="openpyxl") as xw:
        readme(notes).to_excel(xw, sheet_name="README", index=False)
        data.to_excel(xw, sheet_name="data", index=False)
    print("wrote", path, len(data), "rows")


def first_existing(*paths):
    for p in paths:
        if Path(p).exists():
            return Path(p)
    return None


def selection_table(tag):
    sel = pd.read_csv(CAMP / f"genera_{tag}" / "selected_with_taxonomy.tsv", sep="\t")
    sel = sel.rename(columns={"representative_accession": "uniprot_accession", "round": "selection_round",
                              "plddt_global": "afdb_v6_global_plddt"})
    sel["selection_round"] = sel["selection_round"] + 1
    sel["afdb_id"] = "AF-" + sel["uniprot_accession"].astype(str) + "-F1"
    ann = first_existing(CAMP / f"annotation_{tag}" / f"selection_{tag}_annotated.tsv")
    if ann is not None:
        a = pd.read_csv(ann, sep="\t")
        acc = "accession" if "accession" in a.columns else ("representative_accession" if "representative_accession" in a.columns else a.columns[0])
        keep = [c for c in a.columns if c != acc and (("name" in c.lower() and "organism" not in c.lower()) or c in (
            "uncharacterized", "poorly_characterized", "has_interpro", "has_function_cc", "reviewed", "annotation_score"))]
        a = a[[acc] + keep].drop_duplicates(acc).rename(columns={acc: "uniprot_accession"})
        sel = sel.merge(a, on="uniprot_accession", how="left")
    pf = first_existing(CAMP / f"pfam_{tag}" / f"selection_{tag}_pfam.tsv")
    if pf is not None:
        p = pd.read_csv(pf, sep="\t")[["seq_id", "n_domain_hits", "n_unique_pfam", "best_pfam", "best_name"]]
        p = p.rename(columns={"n_domain_hits": "pfam_domain_hits", "n_unique_pfam": "pfam_unique_families",
                              "best_pfam": "pfam_best_accession", "best_name": "pfam_best_family"})
        sel = sel.merge(p.drop_duplicates("seq_id"), on="seq_id", how="left")
    ages = first_existing(CAMP / f"genera_ages_{tag}" / f"sel{tag}_genera_ages.csv")
    if ages is not None:
        g = pd.read_csv(ages)[["db_accession", "phylostratum", "taxonomic_representativeness", "age_bin", "neff", "frac_disordered"]]
        g = g.rename(columns={"db_accession": "uniprot_accession", "age_bin": "genera_age_bin",
                              "phylostratum": "genera_phylostratum", "neff": "neff_uniclust30",
                              "frac_disordered": "aiupred_disordered_fraction"})
        sel = sel.merge(g.drop_duplicates("uniprot_accession"), on="uniprot_accession", how="left")
    return sel


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    s1 = pd.read_csv(OUT / "original_selection_900_replayed.tsv", sep="\t")
    s1 = s1.rename(columns={"db_accession": "uniprot_accession", "plddt": "afdb_v4_global_plddt",
                            "pdbx_PDB_id_code": "afdb_id", "organism_scientific": "organism"})
    write_xlsx(OUT / "Dataset_S1_original_selection_850.xlsx", s1, [
        "Dataset S1. Original AI-guided selection of plant proteins (main-text Fig. 3a-e).",
        "Candidates: 161,460 AFDB v4 models of A. thaliana, G. max, O. sativa subsp. japonica and Z. mays (ESM-2 650M mean embeddings).",
        "Reference: ESM-2 embeddings of PDB sequences of 24 Jan 2025 excluding plant structures released 2021-01-01 to 2025-01-24 (346,853 unique masked sequences).",
        "Score per round: (L2-normalized nearest-reference distance x 100) / (L2-normalized mean distance to all remaining plant proteins x pLDDT); top 50 per round moved into the reference.",
        "18 rounds were run (900 proteins); the first 17 rounds (in_850_set = TRUE, 850 proteins) correspond to Fig. 3c/e.",
        "Rows were reconstructed from the saved per-round indices and validated against the saved per-round pLDDT vectors (exact match for all 18 rounds).",
        "Columns: selection_order, round, rank_in_round, score, nearest_pdb_distance (Euclidean, ESM-2), average_homolog_distance, UniProt accession, AFDB v4 pLDDT, organism, gene name, AFDB id, in_850_set.",
    ])
    s850 = s1[s1.in_850_set]
    print("S1 850: organisms", s850.organism.value_counts().to_dict())
    print("S1 850: median pLDDT %.1f, nearest %.2f, homolog %.2f" % (
        s850.afdb_v4_global_plddt.median(), s850.nearest_pdb_distance.median(), s850.average_homolog_distance.median()))
    for tag, rounds in (("3800", "38 rounds x 100"), ("690", "69 rounds x 10")):
        t = selection_table(tag)
        write_xlsx(OUT / f"Dataset_S{2 if tag == '3800' else 3}_sel{tag}.xlsx", t, [
            f"Dataset S{2 if tag == '3800' else 3}. Revised AI-guided selection sel{tag} ({rounds}; Section S16).",
            "Pool: 182,501 Viridiplantae proteins of the merged AFDB v6 set without a SIFTS-mapped PDB chain (179,596 unique sequences; one accession per sequence, lowest pLDDT).",
            "Reference library: ESM-2 650M embeddings of all 186,823 unique PDB protein sequences released to 30 Sep 2026, plus previously selected proteins.",
            "Score per round: (L2-normalized nearest-library distance x 100) / (L2-normalized mean distance to the remaining pool x pLDDT).",
            "Columns: AFDB v6 sequence id, UniProt accession, selection round, organism, genus, species, AFDB v6 global pLDDT, AFDB id;",
            "UniProt annotation (protein name, uncharacterized / poorly characterized flags, InterPro, function comment, reviewed, annotation score; Oct 2026);",
            "GenEra phylostratum, taxonomic representativeness and age bin (all_life / eukaryotes / kingdom_only / family_or_genus / flagged); Neff (OpenProteinSet UniClust30); AIUPred disordered fraction.",
        ])
        print(tag, "rows", len(t), "cols", list(t.columns))
    acc850 = set(s850.uniprot_accession)
    for tag in ("3800", "690"):
        t = pd.read_csv(CAMP / f"genera_{tag}" / "selected_with_taxonomy.tsv", sep="\t")
        print(f"overlap original-850 with sel{tag}:", len(acc850 & set(t.representative_accession)))


if __name__ == "__main__":
    main()
