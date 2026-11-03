#!/usr/bin/env python3
"""Analysis server: per-protein tables behind PNAS concise SI Figs. S4 and S10 that exist only as raw inputs or summaries.
Uses the same loading code as the figure scripts. Output: /path/to/server_data/si_pnas/datasets/ (originals untouched).

  esmfold_vs_afdb_v6_per_protein.tsv   Fig. S4: AFDB v6 global pLDDT and ESMFold pLDDT (ESM Metagenomic Atlas, full
                                       precision, 0-100) for the 92,385 identical sequences plotted
  published_gene_ages_afdb_v6.tsv      Fig. S10: published GenEra gene ages (Barrera-Redondo et al., 8 species) joined
                                       to AFDB v6 proteins by UniProt accession (all joined rows; the figure uses the
                                       age bins all_life, eukaryotes, kingdom_only, family, species)
"""
import sys
from pathlib import Path

import pandas as pd

OUT = Path("/path/to/server_data/si_pnas/datasets")
OUT.mkdir(parents=True, exist_ok=True)


def esmfold():
    sys.path.insert(0, "/path/to/server_data/esm_atlas_plddt/scripts")
    import plot_precise_esm as pe
    hashes = pd.read_csv(pe.ROOT / "hashes.tsv", sep="\t",
                         dtype={"sequence_hash": str, "db_accession": str, "tax_group7": str})
    scores = pe.load_means()
    finite = scores.mean_plddt.dropna()
    scale = 100.0 if len(finite) and float(finite.max()) <= 1.5 else 1.0
    scores["esmfold_plddt"] = scores.mean_plddt * scale
    keep = ["sequence_hash", "esmfold_plddt"] + [c for c in ("ptm", "source") if c in scores.columns]
    df = hashes.merge(scores[keep].rename(columns={"source": "atlas_source"}), on="sequence_hash", how="left")
    df = df[df.esmfold_plddt.notna() & df.plddt_global.notna()]
    df = df.rename(columns={"plddt_global": "afdb_v6_plddt_global"})
    df.to_csv(OUT / "esmfold_vs_afdb_v6_per_protein.tsv", sep="\t", index=False)
    print("esmfold rows", len(df), df.tax_group7.value_counts().to_dict(), flush=True)


def published_ages():
    sys.path.insert(0, "/path/to/server_data/conservation/scripts")
    import plot_published_suite as ps
    parent, rank, names = ps.mga.load_taxdump(ps.TAXDUMP)
    meta = ps.load_meta()
    ages = ps.assign_names(ps.barrera_genes(ps.load_sgd()), parent, rank, names)
    parts = []
    for proteome, *_rest in ps.BARRERA:
        hit = ps.join_meta(ages[ages.proteome.eq(proteome)], meta, proteome)
        parts.append(hit.sort_values("rank").drop_duplicates("db_accession"))
    df = pd.concat(parts, ignore_index=True)
    cols = ["proteome", "species_label", "source", "gene", "acc", "db_accession", "query_taxid", "phylostratum", "rank",
            "taxonomic_representativeness", "age_bin_lumped", "age_bin", "age_taxid", "plddt_global", "neff",
            "frac_disordered", "tax_group7"]
    df = df[[c for c in cols if c in df.columns]]
    df.to_csv(OUT / "published_gene_ages_afdb_v6.tsv", sep="\t", index=False)
    print("published ages rows", len(df), df.proteome.value_counts().to_dict(), flush=True)


if __name__ == "__main__":
    for job in sys.argv[1:] or ["esmfold", "published_ages"]:
        {"esmfold": esmfold, "published_ages": published_ages}[job]()
        print("DONE", job, flush=True)
