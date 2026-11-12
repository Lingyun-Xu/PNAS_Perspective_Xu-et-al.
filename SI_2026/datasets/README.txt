Datasets of the Supporting Information (PNAS manuscript 2026-25352)
===================================================================================

One file or folder per dataset, numbered as in the SI (Datasets S1-S17). All tables are copies of the files of the original analysis;
files copied between machines were checked by MD5 after transfer. Tab-separated unless .parquet, .csv or .xlsx.

Dataset_S1_original_selection_850.xlsx   Original AI-guided selection (main-text Fig. 3): 850 proteins, with the 50 of
                                          round 18 (900 rows; sheet README explains the columns).
  Dataset_S1_extras/original_selection_900_replayed.tsv   the same 900 rows as a TSV.
Dataset_S2_sel3800.xlsx                   sel3800, with annotation, Pfam domains, disorder, Neff, AF2-style Neff and gene age.
Dataset_S3_sel690.xlsx                    sel690, same columns.

Dataset_S4_AFDB_v4_48_proteomes/AF2_48_model_oganisms.parquet
    AlphaFold DB version 4, the 48 model-organism proteomes of the original study: 564,446 models of 561,559 UniProt
    accessions. Columns: AAseq_one_letter_code, plddt (global pLDDT, _ma_qa_metric_global.metric_value), db_accession,
    db_code, gene_name, ncbi_taxonomy_id, organism_scientific, EMBI_EBI_sequence_checksum, template_ids,
    pdbx_PDB_id_code (AFDB entry identifier). Extracted by reading the gzipped mmCIF files (*.cif.gz) inside the AFDB
    v4 proteome archives as text. Sequences spanning several mmCIF lines keep a leading ";" (strip before use).
Dataset_S5_AFDB_v4_plant_proteomes/AF2_Plants.csv
    The four plant proteomes of the same extraction (Arabidopsis, soybean, rice, maize): 166,181 models.

Dataset_S6_AFDB_v6_47_proteomes/AF2_48_proteomes_v6.parquet (+ _summary.csv)
    AlphaFold DB version 6 (downloaded 5-6 October 2026), the 47 proteomes available in v6: 558,194 models.
Dataset_S7_AFDB_v6_SwissProt/AF2_swissprot_v6.parquet (+ _summary.csv)
    The AFDB v6 Swiss-Prot archive: 550,122 models.
    Columns of S6 and S7: afdb_entry_id, pdbx_PDB_id_code, model_version, db_accession, db_code, gene_name,
    ncbi_taxonomy_id, organism_scientific, seq_db_sequence_checksum, seq_db_isoform, uniprot_start, uniprot_end,
    uniprot_sequence_version_date, AAseq_one_letter_code, length, plddt_global, plddt_mean_residue,
    frac_plddt_very_low / low / confident / very_high, plddt_per_residue, template_ids, template_chains, n_templates,
    proteome_tar, parse_error, msa_raw_depth, is_swissprot, superkingdom, tax_domain, tax_group6, tax_group_fine,
    tax_kingdom, tax_phylum, tax_class, tax_group7. In S7, superkingdom is filled only for entries of the proteome
    species and tax_group7 is empty; tax_group6 is the lineage-based group of every entry (the seven groups of the
    figures add the parasitic-worm species to tax_group6 at plotting time).

Dataset_S8_PDB_RCSB_API_2025-01-24/fetch_API_for_PDB20250124_seq_and_posi.parquet
    PDB sequences fetched through the RCSB API, 24 January 2025 (PDB snapshot of the original selection): 906,025
    entities. Columns: pdb_id, entity_id, sequence, instance_id, unobserved_residue_xyz (positions of residues without
    coordinates, used to mask unobserved residues).
Dataset_S9_PDB_RCSB_API_2024-12-31/fetch_API_for_PDB20241231_seq_and_posi_123024.parquet
    The same for 31 December 2024: 944,532 entities.

Dataset_S10_alignment_depth/merged_neff.parquet
    991,829 proteins of the merged AFDB v6 set (47 proteomes + Swiss-Prot): db_accession, afdb_entry_id, length,
    tax_group7, source, plddt_global, msa_raw_depth (AFDB MSA depth), neff and neff_per_len (OpenProteinSet UniClust30
    alignments; 80% identity), AF2-style_Neff (the same alignments, as defined in the AlphaFold2 paper: 80% identity
    over the region non-gap in either row, effective count per column, median over columns;
    scripts/fetch_af2style_neff.py), n_sequences, query_length, representative / rep_id / how / pident / qcov (how
    each protein received its alignment), has_msa.
Dataset_S11_ESMFold_confidence/esmfold_vs_afdb_v6_per_protein.tsv (+ correlations)
    92,385 proteins of the 47 proteomes with an ESM Metagenomic Atlas model: db_accession, afdb_entry_id, tax_group7,
    proteome_tar, afdb_v6_plddt_global, seq_len, sequence_hash, esmfold_plddt (0-100), ptm, atlas_source.
Dataset_S12_disorder/merged_per_protein.parquet (+ D1_binned.csv, D2_per_proteome_medians.csv,
    disorder_by_superkingdom.csv, regressions.csv, adopt_subset_merged.csv, summary.txt)
    558,194 entries of the 47 v6 proteomes with AIUPred disorder (residue disordered if score > 0.5): frac_disordered,
    mean_disorder, n_idr_30, longest_idr, frac_ordered_residues, plddt_all, plddt_ordered, plddt_disordered,
    plddt_global, msa_raw_depth, n_templates, length.
    disorder_vs_template_correlations.tsv, difference_from_plants_matched.tsv (with the values for each plant
    proteome) and difference_from_plants_sensitivity.tsv: the comparison of disorder with template distance and MSA
    depth (comprehensive SI only), from disorder_template_analysis.py with the template distances of Dataset S17
    (FigS06_S07/template_distances.parquet).
    disorder_by_group_descriptive.tsv, disorder_group_tests.tsv, disorder_among_plant_dunn.tsv and
    disorder_pLDDT_correlation_tests.tsv: the statistics of the disorder table of the concise PNAS SI (Table S19;
    seven groups, 554,943 canonical entries), from disorder_group_stats.py.
    Not used in the SI: regressions.csv (an earlier regression analysis, replaced by the comparison above),
    adopt_subset_merged.csv (ADOPT disorder for a stratified subset of 18,160 proteins) and the ADOPT part of
    summary.txt; kept for a possible further check with a second predictor.
Dataset_S13_homolog_class/panel_strata.parquet (+ panel_bin_summary.csv)
    554,943 proteins of the 47 proteomes: panel_bin (homolog class), best_outgroup_evalue (MMseqs2), plddt_global, neff,
    AF2-style_Neff, frac_disordered, tax_group7, proteome, family_flag, polyploid_proteome.
Dataset_S14_published_gene_ages/published_gene_ages_afdb_v6.tsv (+ suite_bin_summary.csv)
    210,457 proteins of 9 species: published gene ages joined to AFDB v6 by UniProt accession: phylostratum, rank,
    taxonomic_representativeness, age_bin, plddt_global, neff, AF2-style_Neff, frac_disordered. The GenEra maps of
    eight species (Barrera-Redondo et al. 2023; 190,434 proteins) and, as source 'Lotharukpong 2023', the human map
    of the phylomapr R package (Homo_sapiens.PhyloMap; GenEra, DIAMOND sensitive mode; numbered phylostrata only;
    20,023 proteins): rank k is named after the k-th of the 28 GenEra phylostrata of the human NCBI lineage
    (cellular organisms ... Homo sapiens; Craniata, Teleostomi and the genus omitted, as GenEra does; the same in
    all six NCBI Taxonomy releases checked, January 2019 to October 2026) and binned with the same code as the
    other species (Fig. S8; scripts figS8_published_ages_9species.py and table_s19_human.py). The human map has no
    taxonomic representativeness values (column empty).
    genera_balanced_sample_ages.csv (+ genera_balanced_sample_summary.csv): our GenEra ages (v1.4.3; NCBI nr of
    1 October 2026) of the balanced sample of 50 random proteins per proteome (seed 20261009), 2,306 proteins with a
    GenEra result: proteome, gene, organism, group, phylostratum, rank, taxonomic representativeness, age bin, pLDDT,
    Neff, AF2-style Neff, disordered fraction.
Dataset_S15_compared_protein_sets/per_protein_comparison.tsv (+ comparison_statistics.tsv)
    18,556 rows, one per protein of each compared set (sel3800 3,800; sel690 690; pdb3825 full-length 2,644 and
    constructs 3,825; random pool 3,800; random with the sel3800 species mix 3,797), with every measure of Section S7.
    Embedding measures against the 182,286 unique plant sequences of the AFDB v6 set: d_nearest_plant (pairwise distance
    to the nearest other plant sequence; Fig. S21I and Table S36), d_mean_plant (mean distance to all other plant
    sequences; Fig. S21J), both with identical sequences not counted; and, from earlier versions, d_homolog_mean (the
    same mean including the identical sequence) and n_within_1.0 / 1.5 / 2.0 (plant sequences within these distances).

Dataset_S16_source_data_main_text_figures/   Source data of main-text Figs. 1-3 (Fig. 4a reuses Fig. 3; Fig. 4b is a
    schematic), one folder per figure, located in the original notebooks (README.txt inside lists the files, columns,
    panels and notebook cells). Fig. 2b: the hand assignment of 9,451 of the 176,322 PDB sequences was not saved.

Dataset_S17_source_data_SI_figures/   Source data of SI Figs. S1-S21, one folder per figure (README.txt inside lists
    the files, columns and panels; figures drawn from Datasets S10-S15 point to them).
