Dataset S17. Source data of the SI figures (concise PNAS SI, Figs. S1-S21)
==========================================================================

Each folder holds the table(s) one SI figure was drawn from. Figures drawn from another dataset of the SI point to it
instead of duplicating it. Unless stated otherwise, pLDDT is the AFDB v6 global pLDDT, groups (tax_group7) are the seven
taxonomic groups of Section S1.2, and distances are Euclidean distances between ESM-2 (650M) mean embeddings.
Plotting code: scripts/analysis/si_pnas_server.py (Figs. S1, S3-S5, S9, S10), si_pnas_cluster.py (S6, S7, S15-S18),
figS11_S13_candidate_pool.py (S11-S13) and si_pnas_local.py (S2, S8 composite, S14, S20, S21). Files are tab-separated
unless they are .parquet.

Fig. S1  Global pLDDT density by group (AFDB v6 merged set)
  -> Dataset S10 (Dataset_S10_alignment_depth/merged_neff.parquet): columns plddt_global, tax_group7; 991,829 proteins.

Fig. S2  AFDB v4 versus v6 global pLDDT
  FigS02/plddt_v4_vs_v6_entries.tsv.gz  503,005 matched entries: entry, db_accession, organism_v4, tax_group7, plddt_v4,
                                         plddt_v6, diff (v6 - v4), same_sequence, source_v6 (panels A-C)
  FigS02/plddt_v4_vs_v6_summary.tsv     15 subsets: n, Pearson r, Spearman rho, Lin's CCC, mean difference with 95% CI,
                                         limits of agreement, share changing > 5 units, share crossing 70 (panel D; Table S6)

Fig. S3  Mapped pLDDT versus lDDT-Calpha (AFDB v6 models versus post-2018 experimental structures)
  FigS03/lddt_primary.tsv  12,222 eligible model-structure pairs (status = ok for the 8,181 scored): db_accession,
                           n_ca, n_lddt, lddt_ca, plddt_mapped, rmsd_ca, plddt_global, tax_group7, proteome, PDB, CHAIN,
                           n_obs, emethod

Fig. S4  AF2 (AFDB v6) versus ESMFold pLDDT
  -> Dataset S11 (Dataset_S11_ESMFold_confidence/esmfold_vs_afdb_v6_per_protein.tsv): afdb_v6_plddt_global,
     esmfold_plddt, tax_group7; 92,385 proteins.

Fig. S5  Global pLDDT versus alignment depth
  -> Dataset S10 (merged_neff.parquet): msa_raw_depth (A, AFDB MSA depth), neff (B, plotted as log10(neff + 1)),
     neff_per_len (C), plddt_global, tax_group7.

Figs. S6 and S7  Template distance and nearest-PDB distance
  FigS06_S07/template_distances.parquet  991,829 proteins of the merged AFDB v6 set: plddt_global, tax_group7, in_pdb (the
     protein has its own PDB structure), n_templates, d_template_min (closest template used by the model), d_nearest_asof
     (nearest PDB chain released on or before 3 February 2021, own structures removed), d_nearest_now (September 2026
     PDB), with the identifiers of the closest template and nearest chains.
     Fig. S6: A, in_pdb = false and d_template_min defined (906,513); B, in_pdb = true (39,524); C, both (946,037).
     Fig. S7: d_template_min versus d_nearest_asof (A-C) or d_nearest_now (D-F), same subsets.
     Tables S16-S18: in_pdb = false; d_template_min for the 906,513 models with it defined, d_nearest_asof and
     d_nearest_now for all 951,630 (45,090 built without templates and 27 whose templates could not be mapped have
     no d_template_min). Computed with dtemplate_stats.py (scripts/analysis; results in statistics/dtemplate_*.tsv).

Fig. S8  Intrinsic disorder (AIUPred)
  -> Dataset S12 (Dataset_S12_disorder): merged_per_protein.parquet (558,194 entries; B and C), D1_binned.csv (A),
     D2_per_proteome_medians.csv (D), disorder_by_superkingdom.csv.

Fig. S9  Homolog class within the 47 proteomes
  -> Dataset S13 (Dataset_S13_homolog_class): panel_strata.parquet (554,943 proteins), panel_bin_summary.csv.

Fig. S10  Published gene ages joined to AFDB v6
  -> Dataset S14 (Dataset_S14_published_gene_ages): published_gene_ages_afdb_v6.tsv (190,434 proteins), suite_bin_summary.csv.

Figs. S11-S13 share the definitions of the selection (select_plants.py):
  candidate pool     the 179,596 unique sequences of the 182,501 AFDB v6 plant accessions without a SIFTS-mapped PDB chain,
                     in the pool order of the selection (pool_position); in a selection round, the candidates still in
                     the pool
  reference library  the 186,823 unique protein sequences of the September 2026 PDB (complete PDB), or the 182,998
                     without the 3,825 plant sequences first released after 3 February 2021 (pdb3825), plus added proteins
  nearest distance   raw minimum ESM-2 distance from a candidate to the reference library
  average homolog distance   mean ESM-2 distance from a candidate to all candidates of the pool (in a selection round,
                     to those still in the pool), including the zero distance to itself
  Figures: scripts/analysis/figS11_S13_candidate_pool.py <this dataset folder> <output folder>. Pool positions are not
  linked to accessions, except in FigS13 (seq_id) and for the selected proteins (selected_accessions_*.tsv).

Fig. S11  Candidate pool against the PDB with and without pdb3825
  FigS11/candidate_pool_round0.tsv.gz  the 179,596 candidates: d_nearest_complete (to the complete PDB; A),
                                       d_nearest_without_pdb3825 (B), d_average_homolog, nearest_is_pdb3825 (black dots
                                       in A: the closest library sequence is a pdb3825 sequence). Same values as round 0
                                       of the selections (FigS12); from counterfactual_selection.py (mins_round0.npz,
                                       pool_sums_round0.npy).

Fig. S12  Selection rounds (sel3800, A-C; sel690, D-F)
  FigS12/sel3800_rounds_0_and_37.tsv.gz  every candidate still in the pool at round 0 (179,596) and round 37, the last
                                         (175,896): pool_position (position in that round's pool), d_nearest_library (to
                                         the complete PDB plus the proteins selected in earlier rounds),
                                         d_average_homolog, selected_this_round
  FigS12/sel690_rounds_0_and_68.tsv.gz   the same for sel690, round 0 (179,596) and round 68, the last (178,916)
  FigS12/selected_accessions_sel3800.tsv  all 3,800 selected sequences (3,801 accessions): round, rank_in_round, score,
                                         min_dist and avg_dist when selected, pLDDT, organism
  FigS12/selected_accessions_sel690.tsv  the same for the 690 sequences of sel690
  FigS12/limits.json                     histogram range of the density differences (Figs. S11 and S12)

Fig. S13  Real versus AI-selected additions (3,800 comparison)
  FigS13/evaluation_distances.parquet  175,537 candidates selected in no arm: pool_position, seq_id, base_3800 (nearest
                                       distance to the PDB without pdb3825; A), real_3800 (with the 3,825 real
                                       sequences added = complete PDB; B), ai_3800 (with the 3,800 proteins selected on
                                       the baseline library added; C), the same for the 690 comparison,
                                       d_average_homolog (y of all panels)
  FigS13/real_vs_ai_metrics.tsv        the comparison statistics by species group (Table S28)
  FigS13/summary.json                  reproduction checks and set overlaps

Fig. S14  UMAP and HDBSCAN maps
  FigS14/hdbscan_on_umap*.tsv  umap1, umap2 and HDBSCAN cluster (-1 = noise) for pdb3825 (A, K97 file), sel3800
                               (B, K109), sel690 optimum (C, K5) and finer view (D, K10), with the protein identifiers

Figs. S15-S18  HHsearch and Foldseek top hits before and after additions
  FigS15_S18/{hhsearch,foldseek}_cohort.tsv   per protein: top hit and its ESM-2 distance in each of the four searches
                                              (pdbonly = September 2026 PDB; minus3825 = without pdb3825; plus3800,
                                              plus690 = with the selected proteins added); empty = no hit
  FigS15_S18/{search}_{set}_pairs.tsv         the proteins of each comparison: d_before, d_after, top_is_added, top_after,
                                              d_homolog (S15 hhsearch_sel3800, S16 hhsearch_pdb3825, S17 foldseek_sel3800,
                                              S18 foldseek_pdb3825; the sel690 files back Table S29)
  FigS15_S18/summary.tsv, new_coverage.tsv    the values of Tables S29 and S30 and the first-hit counts

Fig. S19  Alignment depth of the compared sets
  FigS19/neff_sets.tsv        14,731 set members (sel3800, sel690, full-length pdb3825, two random sets): neff, log10_neff
                              (empty when the protein has no alignment)
  FigS19/neff_sets_tests.tsv  the values of Table S31

Fig. S20  Gene age of the selections
  FigS20/sel3800_genera_ages.csv, sel690_genera_ages.csv  GenEra phylostratum, rank, taxonomic representativeness,
                                                           age_bin, plddt_global, neff, frac_disordered
  -> pdb3825 ages and the other set measures: Dataset S15.

Fig. S21  Per-crop comparison of the sets
  -> Dataset S15 (Dataset_S15_compared_protein_sets/per_protein_comparison.tsv): one row per protein of each set (18,556).
