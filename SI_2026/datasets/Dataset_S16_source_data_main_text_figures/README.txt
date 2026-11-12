Dataset S16. Source data of the main-text figures
=================================================

Current main text (rebuttal, November 2026) -> source data (Table S36 of SI_PNAS_rebuttal_11062026 gives every file):
  Fig. 1 (AFDB v6)  drawn from Datasets S10, S12 and S17, no folder here: a, Dataset S10 merged_neff.parquet
                    (plddt_global, tax_group7); b and d-g, Dataset S17 FigS06_S07/template_distances.parquet
                    (d_template_min, the closest-template distance, models with n_templates > 0); c and h-k, Dataset S10
                    (AF2-style_Neff / query_length, with plddt_global); the species of b-k (four plant proteomes, E. coli,
                    yeast, human) from the proteome column of Dataset S12 merged_per_protein.parquet
  Fig. 2a           Fig2a/   UniProtKB release 2026_03
  Fig. 2b           Fig2b/   Fig2b_pdb_unique_sequences.tsv.gz and Fig2b_source_organism_mapping.tsv (September 2026)
  Fig. 2c           Fig2c/   all 1,841,285 mega-scale rows (Fig2c_thermostability_counts.tsv and the files below it)
  Fig. 2d, 2e       Fig2d/, Fig2e/
  Fig. 3            Dataset S17 FigS13/evaluation_distances.parquet (base_3800, real_3800, ai_3800, d_average_homolog),
                    no folder here
  Fig. 4 (AFDB v6)  Fig4_AFDB_v6/  Fig4_v6_source_data.tsv.gz (a-o), Fig4_v6_selected_candidates.tsv.gz,
                    Fig4_v6_summary.tsv (November 2026; candidate pool of main-text Fig. 3)
  Fig. 4 (first version, AFDB v4 plant sequences)  Fig3/  Fig3f-j_thermostability.tsv.gz (a-e),
                    Fig3k-o_MaveDB.tsv.gz (f-j), Fig3p-t_mdCATH.tsv.gz (k-o): the thermostability, MaveDB and
                    mdCATH rows of the first version of Fig. 3
  Fig. 5            schematic, no data
The folders Fig1a, Fig1b_Fig1h-k, Fig1c and Fig1d-g, the panel a and b counts of Fig2/Fig2_plotted_counts.tsv,
Fig2b/unique_protein_sequences.json.gz and Fig3/Fig3a-e*, Fig3b* hold the source data of the first version of the
figures (AFDB v4; July 2026 pies; first selection) and are described below as they were.

Fig2a/   uniprotkb_taxonomy_counts.json  UniProtKB entries (Swiss-Prot + TrEMBL) per taxonomy query of the UniProt REST
         API, release 2026_03 (2 September 2026), retrieved 1 November 2026: counts, query URLs and the plotted
         categories ("figure": Bacteria, Archaea, Viruses, Fungi, Viridiplantae, Insecta, Mammalia (H. sapiens included),
         Other Vertebrata, Other); 150,006,383 entries, Viridiplantae 17,434,591 (11.6%)
Fig2b/   (September 2026 panel)
         Fig2b_pdb_unique_sequences.tsv.gz  the 184,499 of the 186,823 unique protein sequences of the PDB released up to
                                            30 September 2026 (the September 2026 PDB template pool of the SI) that have a
                                            source organism: embed_id (sequence identifier of the template pool),
                                            source_organism (the most frequent of its chains), taxid (NCBI Taxonomy),
                                            exact_name_match, category (group); Viridiplantae 8,371 (4.5%)
         Fig2b_source_organism_mapping.tsv  the 12,126 organism names: taxid, exact_name_match, category

Notes on the first version of the figures
-----------------------------------------
Each folder holds the table(s) one main-text figure was drawn from, located in the 2024-2026 notebooks of the original
analysis and copied from there or from MD5-identical copies (MD5 checked after every transfer). PyTorch tensors (.pt) were converted to
tab-separated tables. Distances are Euclidean distances between ESM-2 (650M) embeddings; pLDDT is the AFDB v4 global
pLDDT.

Fig1a/   global pLDDT density and box plot by group (notebook plddt_vs_species.ipynb, cells 1, 29, 42, 43)
  Fig1a_global_plddt_by_group.tsv.gz  564,446 AFDB v4 models of the 48 proteomes: afdb_entry_id, db_accession,
                                      organism_scientific, group, global_plddt (Animalia 204,479; Viridiplantae 166,181;
                                      Fungi 77,026; Bacteria 61,727; Protist 53,260; Archaea 1,773). Built from Dataset S4
                                      and the group table below, as in cell 29.
  group_of_each_proteome.tsv          the group of each of the 48 proteomes (hand assignment from NCBI Taxonomy)

Fig1b_Fig1h-k/   global pLDDT versus MSA depth (b) and MSA-depth distributions with Jensen-Shannon distances (h-k)
  MSA_depth_with_plddt_34921.csv  34,921 rows: Neff, Neff/len (the "MSA depth" of the figure), Total sequences, Source
                                  (PDB chain or UniClust30 alignment), db_accession, organism_scientific, plddt
                                  (notebook Refer_back_to folder_names_and_calculate_MSA_depth.ipynb, cells 33, 44, 47-50;
                                  organisms other than human, E. coli and yeast are the "plant" group). Section S4.1 of the
                                  SI explains why this table gives a weak correlation; Table S10 lists its statistics.

Fig1c/   global pLDDT versus distance to the nearest PDB sequence (May 2024 PDB; notebook
         Plot_Plants_that_have_PDB.ipynb, cells 5, 7, 18-20, 27, 28)
  AF2plant_exclude_crystalizedentries_to_redoPDB_distance.parquet  161,074 plant proteins: Plants_protein_name,
                                  nearest_PDB_protein_name, distance_to_nearest_PDB
  ECOLI_ / HUMAN_ / YEAST_to_redoPDB_distance.parquet              4,321 / 23,322 / 5,973 proteins, same columns
  AF2_ECOLI.csv, AF2_HUMAN.csv, AF2_YEAST.csv                     AFDB v4 tables with the global pLDDT, joined on accession
  -> plant pLDDT: Dataset S5 (AF2_Plants.csv). Pooled n = 194,690; the regression line was fitted on a random subsample
     of 34,921 (seed 42).

Fig1d-g/   distributions of the nearest-PDB distance and Jensen-Shannon distances, plants versus E. coli, yeast and human
           (February 2025 recomputation; notebook Plots_now_for_paper.ipynb, cells 18, 20, 32, 35, 40)
  nearest_PDB_distance_Feb2025.tsv.gz  species, d_nearest_pdb: Plants 161,460, Ecoli 4,321, Yeast 5,973, Human 23,322
                                       (from min_pairwise_dist_*_021125.pt; values in the order of the embedding sets).
                                       Reference: the complete PDB of 24 January 2025 (the plant values are identical
                                       to round 0 of Fig3b_all_PDB_2025-01-24.tsv.gz).
                                       The 1.9-threshold percentages of the main text come from these values.

Fig2/   composition pies (percentages printed on the figure; counts below)
  Fig2_plotted_counts.tsv  the plotted counts of every panel:
     a  UniProtKB entries, snapshot of 2 July 2026, typed from the UniProt statistics page into Fig2_Pie_chart.ipynb,
        cell 116 (total 149,810,139; Viridiplantae 17,344,395)
     b  PDB unique protein sequences, July 2026 (see Fig2b)
     c  thermostability measurements (rows): FireProtDB + mega-scale dataset (1,102,020) + Meltome Atlas (40,381) =
        1,142,401; Nematoda shown within Others. NOTE: the mega-scale part holds only the first 1,048,575 of the
        1,841,285 rows of Tsuboyama et al. Dataset 1 (spreadsheet row limit; FireProtDB is the other 53,445 rows).
        RECOMPUTED 6 Nov 2026 with all 1,841,285 rows (1,935,111 rows in total; Viridiplantae 2.85% instead of
        4.30%): see Fig2c/Fig2c_thermostability_counts.tsv
     d  MaveDB score sets (August 2024 export, 2,588), groups from the lists of Fig2_Pie_chart.ipynb cell 58; Insecta and
        Nematoda shown within Others (36)
     e  mdCATH: PDB entities of the 5,398 domains (6,361 rows), groups from MD_unique_species_with_superkingdom.csv;
        Nematoda shown within Others
  Fig2b/  unique_protein_sequences.json.gz (177,092 unique PDB protein sequences: sequence, pdb_ids, source_organism,
          species (NCBI description), superkingdom), pdb_superkingdom_still_unmapped_hand_curated.csv (hand assignment
          of 201 species descriptions, 2,144 sequences), nih_species_cache.json (organism -> NCBI description),
          fetch_pdb_unique_sequences.py (the fetch and automatic assignment).
          NOTE: 176,322 sequences have a source organism. 164,727 were assigned automatically and 2,144 through the
          hand-curated table; the hand assignment of the other 9,451 was not saved, so the printed percentages (also in
          Fig2_plotted_counts.tsv) cannot be recomputed exactly. With the recorded assignments, Viridiplantae is 8,000
          sequences (4.54% of 176,322; 4.6% on the figure).
  Fig2c/  firemega_complete_species.csv.gz (FireProtDB + mega-scale rows with species and group; the first version),
          meltome_altas_with_organisms_species_superkindoms.csv.gz (Meltome Atlas rows with species and group)
          Recomputed panel (updated_maintext_figure/scripts/fig2c_all_rows.py; images Fig2c_thermostability.png and
          _labeled.png in updated_maintext_figure/Fig2):
          Fig2c_thermostability_counts.tsv  rows per group for FireProtDB, the mega-scale dataset (all rows) and the
                                            Meltome Atlas, with the shares of the first version for comparison
          Fig2c_mega_entry_groups.tsv.gz    each mega-scale entry (PDB entry, UniProt or GenBank accession, or design
                                            family): group, how it was assigned (published mapping, RCSB source
                                            organism, design pattern, curated), rows in all 1,841,285 and in the first
                                            1,048,575. The first 1,048,575 rows keep their published groups, except the
                                            1,648 rows of 1E68 and 3E21, whose names had been read as numbers ('1', '3')
                                            and counted as Synthetic. 4,099 rows (0.21%) without an identifiable
                                            organism (mostly Pfam-family sequences) are drawn within Other
          fig2c_named_entries.tsv           the ten named natural peptides assigned by hand, with organism and source
  Fig2d/  mavedb_organisms_score_sets_groups.tsv (99 organisms: score sets, NCBI description, group)
  Fig2e/  MD_PDB_entries_with_superkingdom.parquet (6,361 PDB entities of the mdCATH domains; its own superkingdom
          column is outdated: re-join on organism_name), MD_unique_species_with_superkingdom.csv (1,260 species ->
          group)

Fig3/   free-energy maps of the plant proteins (nearest-reference distance versus average homolog distance) during the
        AI-guided selections, one file per reference set; each file holds the plotted rounds:
  Fig3a-e_PDB_without_2021-2025_plant_structures.tsv.gz  rounds 0 (a) and 17 (c; after 17 rounds x 50 = 850 proteins,
        Dataset S1); reference = PDB of 24 January 2025 without the plant structures released since 2021
  Fig3b_all_PDB_2025-01-24.tsv.gz                         round 0 (b); reference = the complete PDB of 24 January 2025
  Fig3f-j_thermostability.tsv.gz                           rounds 0, 1 and 20 (f, g, h; 500 proteins per round)
  Fig3k-o_MaveDB.tsv.gz                                    rounds 0, 1 and 40 (k, l, m; 50 per round)
  Fig3p-t_mdCATH.tsv.gz                                    rounds 0, 1 and 20 (p, q, r; 500 per round)
  Columns: round, pool_position (position in that round's pool of remaining plant proteins; not linked to
  accessions), d_nearest_reference (distance to the nearest sequence of the reference set, including the proteins
  moved into it), d_average_homolog (mean distance to the plant proteins remaining in the pool). The difference
  panels (d, e, i, j, n, o, s, t) are differences between these maps. Rounds were matched to the published panels by
  file names and the figure legend (image index = round + 1); the 2026 plotting code of panels a-o was not found.
  The round-17 and round-0 tensors of panels a and c are identical to the earlier copies of these files.

Fig. 4  a reuses the maps of Fig. 3a-e (source data above); b is a schematic without data.

Fig4_AFDB_v6/   main-text Fig. 4 recomputed on the AFDB v6 candidate pool of main-text Fig. 3 (179,596 unique sequences of
                the AFDB v6 Viridiplantae proteins without a SIFTS-mapped PDB chain; SI Section S3). Selection as in the
                first version: score = nearest template distance / average homolog distance (no pLDDT term); per round
                the top 500 (thermostability, mdCATH) or 50 (MaveDB) candidates move into the template pool; panels show
                rounds 0, 1 and 20 (40 for MaveDB). Distances: Euclidean between ESM-2 (650M) mean embeddings (layer
                33, sequences truncated at 1,022 residues). Template pools:
                  thermostability  1,778,010 unique sequences: FireProtDB wild types (242), every protein sequence of
                                   the 1,841,285 rows of mega-scale Dataset 1 (DNA translated, SAGG linker padding
                                   removed; 1,734,297) with the 3,138 sequences found only in Datasets 2 and 3, and the
                                   Meltome Atlas proteins (40,381)
                  thermostability (Datasets 2 and 3 only)  the pool of the first version: FireProtDB, mega-scale
                                   Datasets 2 and 3 (740,697) and Meltome Atlas, 781,273 unique sequences (comparison)
                  MaveDB           1,170 unique target sequences (export of 14 November 2024)
                  mdCATH           5,398 domains (ATOM-record sequences of the CATH v4.2.0 S20 non-redundant set)
  Fig4_v6_source_data.tsv.gz  row, round (0, 1, last), candidate_slot (position in candidates of the run), seq_id,
                              nearest_template_distance (x), average_homolog_distance (y; mean distance to the other
                              candidates still in the pool, own zero distance excluded), for every candidate still in
                              the pool in that round
  Fig4_v6_selected_candidates.tsv.gz  the selected candidates: row, round, rank_in_round, seq_id,
                              representative_accession, organism_scientific, plddt (AFDB v6 global), score, d_min,
                              d_avg (at selection)
  Fig4_v6_summary.tsv         per panel: candidates, median and mean nearest template distance, shares within 1.9 and
                              1.0, median homolog distance, points outside the plotted range, share brought closer
                              than in round 0
