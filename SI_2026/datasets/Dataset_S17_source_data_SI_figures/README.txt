Dataset S17. Source data of the SI figures (concise PNAS SI, Figs. S1-S21)
==========================================================================

Folders are named after the figure numbers of the concise PNAS SI. In the SI for the rebuttal
(SI_PNAS_rebuttal_11062026; Fig. S1 of the concise SI became main-text Fig. 1a and was removed) they hold:
  FigS02      -> Fig. S1  (AFDB v4 versus v6)
  FigS03      -> Fig. S2  (pLDDT versus lDDT-Ca)
  FigS06_S07  -> Figs. S5 and S6 (closest-template and nearest-PDB distances), and main-text Fig. 1b, 1d-g
  FigS11      -> Fig. S10A-C (candidate pool, round 0)
  FigS12      -> Fig. S10D-I (selection rounds; Figs. S10 and S11 of the first rebuttal draft were combined on
                 7 Nov 2026: candidate_pool_and_rounds_AI.png, figS11_S13_candidate_pool.py option s10ai)
  FigS13      -> Fig. S11 (change of the average homolog distance), and main-text Fig. 3 (evaluation_distances.parquet:
                 base_3800, real_3800, ai_3800 and d_average_homolog of the 175,537 evaluated candidates)
  FigS14      -> Fig. S12 (A, B: UMAP-HDBSCAN maps of pdb3825 and sel3800; C, D: joint UMAP coloured by the
                 nearest-neighbour distance; figS12_umap_diversity.py, 7 Nov 2026)
  FigS15_S18  -> Figs. S13-S16 (HHsearch and Foldseek)
  FigS19      -> Fig. S17 (alignment depth of the compared sets)
  FigS20      -> Fig. S18 (gene ages of the selections)
The other SI figures of the rebuttal SI are drawn from Datasets S10-S15 (Table S36 of that SI). The figure numbers
below are those of the concise PNAS SI.

Each folder holds the table(s) one SI figure was drawn from. Figures drawn from another dataset of the SI point to it
instead of duplicating it. Unless stated otherwise, pLDDT is the AFDB v6 global pLDDT, groups (tax_group7) are the seven
taxonomic groups of Section S1.2, and distances are Euclidean distances between ESM-2 (650M) mean embeddings.
Each distance has a name (concise SI Table S1): closest-template distance (d_temp, to the closest template used by the
model), nearest-PDB distance (at model building, or September 2026), nearest template-pool distance (d_min, to the
template pool of the selection), average homolog distance (d_avg, mean distance from a candidate to the other
candidates of the candidate pool, itself not counted), pairwise distance (between two proteins) and top-hit distance
(to the top HHsearch or Foldseek hit).
Protein sets: candidate pool (179,596 unique AFDB v6 plant sequences without a PDB chain; its members are candidates);
September 2026 PDB template pool (186,823 unique PDB sequences); baseline template pool (the same without pdb3825;
182,998); 690-comparison template pool (the same without the 690 most recently released plant sequences; 186,133);
selections, sel3800 and sel690 (chosen against the September 2026 PDB template pool).
d_avg values: corrected on 6 Nov 2026 to exclude the candidate's own zero distance (davg_self_excluded.py on the
cluster, apply_davg_self_excluded.py). The selection code had divided by the number of candidates including the
candidate itself; this rescales all values of a round by the same factor (N - 1)/N and changes no ranking.
Plotting code: SI_2026/scripts/si_pnas_server.py (Figs. S1, S3-S5, S9, S10), si_pnas_cluster.py (S6, S7),
figS11_S13_candidate_pool.py (S11, S12), figS13_average_homolog_distance.py (S13), figS15_S18_search_density.py
(S15-S18) and si_pnas_local.py (S2, S8 composite, S14, S20, S21). Files are tab-separated unless they are .parquet.

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
     AF2-style_Neff (C, plotted as log10(AF2-style_Neff + 1)), neff_per_len (D), plddt_global, tax_group7.

Figs. S6 and S7  Template distance and nearest-PDB distance
  FigS06_S07/template_distances.parquet  991,829 proteins of the merged AFDB v6 set: plddt_global, tax_group7, in_pdb (the
     protein has its own PDB structure), n_templates, d_template_min (closest-template distance d_temp: to the closest
     template used by the model), d_nearest_asof (nearest-PDB distance at model building: to the nearest PDB chain
     released on or before 3 February 2021, the newest template date of the models; own structures removed),
     d_nearest_now (nearest-PDB distance, September 2026), with the identifiers of the closest template and nearest
     chains.
     Fig. S6: A, in_pdb = false and d_template_min defined (906,513); B, in_pdb = true (39,524); C, both (946,037).
     Fig. S7: d_template_min versus d_nearest_asof (A-C) or d_nearest_now (D-F), same subsets.
     Tables S16-S18: in_pdb = false; d_template_min for the 906,513 models with it defined, d_nearest_asof and
     d_nearest_now for all 951,630 (45,090 built without templates and 27 whose templates could not be mapped have
     no d_template_min). Computed with dtemplate_stats.py (SI_2026/scripts; results in SI_2026/statistics/dtemplate_*.tsv).

Fig. S8  Intrinsic disorder (AIUPred)
  -> Dataset S12 (Dataset_S12_disorder): merged_per_protein.parquet (558,194 entries; B and C), D1_binned.csv (A),
     D2_per_proteome_medians.csv (D), disorder_by_superkingdom.csv.

Fig. S9  Homolog class within the 47 proteomes
  -> Dataset S13 (Dataset_S13_homolog_class): panel_strata.parquet (554,943 proteins; C plots AF2-style_Neff),
     panel_bin_summary.csv.

Fig. S10  Published gene ages joined to AFDB v6
  -> Dataset S14 (Dataset_S14_published_gene_ages): published_gene_ages_afdb_v6.tsv (190,434 proteins), suite_bin_summary.csv.

Figs. S11-S13 share the definitions of the selection (select_plants.py):
  candidate pool     the 179,596 unique sequences of the 182,501 AFDB v6 plant accessions without a SIFTS-mapped PDB
                     chain, in the order of the selection (pool_position); in a selection round, the candidates still
                     in the candidate pool
  template pool      the September 2026 PDB template pool (186,823 unique protein sequences) or the baseline template
                     pool (182,998, without pdb3825, the 3,825 plant sequences first released after 3 February 2021),
                     plus added proteins
  d_nearest_*        nearest template-pool distance d_min: raw minimum distance from a candidate to the template pool
  d_average_homolog  d_avg: mean distance from a candidate to the other candidates of the candidate pool (in a selection
                     round, to the other candidates still in it); the candidate's own zero distance is not counted
                     (round 0: sum over the other candidates / 179,595; later rounds of FigS12: the stored round
                     values x N_r / (N_r - 1), N_r = candidates in that round)
  Figures: SI_2026/scripts/figS11_S13_candidate_pool.py <this dataset folder> <output folder> s11 s12. Pool positions
  are not linked to accessions, except in FigS13 (seq_id) and for the selected proteins (selected_accessions_*.tsv).

Fig. S11  Candidate pool against the September 2026 PDB template pool and the baseline template pool
  FigS11/candidate_pool_round0.tsv.gz  the 179,596 candidates: d_nearest_complete (to the September 2026 PDB template
                                       pool; A), d_nearest_without_pdb3825 (to the baseline template pool; B),
                                       d_average_homolog, nearest_is_pdb3825 (black dots in A: the nearest template is
                                       a pdb3825 sequence). Same values as round 0
                                       of the selections (FigS12); from counterfactual_selection.py (mins_round0.npz,
                                       pool_sums_round0.npy).

Fig. S12  Selection rounds (sel3800, A-C; sel690, D-F)
  FigS12/sel3800_rounds_0_and_37.tsv.gz  every candidate still in the candidate pool at round 0 (179,596) and round 37,
                                         the last (175,896): pool_position (position in that round's candidate pool),
                                         d_nearest_library (to the September 2026 PDB template pool plus the proteins
                                         selected in earlier rounds),
                                         d_average_homolog, selected_this_round
  FigS12/sel690_rounds_0_and_68.tsv.gz   the same for sel690, round 0 (179,596) and round 68, the last (178,916)
  FigS12/selected_accessions_sel3800.tsv  all 3,800 selected sequences (3,801 accessions): round, rank_in_round, score,
                                         min_dist (d_min) and avg_dist (d_avg, without the self distance) when
                                         selected, pLDDT, organism
  FigS12/selected_accessions_sel690.tsv  the same for the 690 sequences of sel690
  FigS12/limits.json                     histogram range of the density differences (Figs. S11 and S12)

Fig. S13  Change of the average homolog distance when the candidates of sel3800 or sel690 leave the candidate pool
  FigS13/davg_change.tsv.gz  the 175,537 candidates selected in no arm: seq_id, d_avg (candidate pool of 179,596),
                             d_avg_without_sel3800 (after the 3,800 candidates of sel3800 leave it),
                             d_avg_without_sel690, change_sel3800, change_sel690 (plotted). Adding pdb3825 (PDB
                             sequences, not candidates) changes no d_avg. Written by apply_davg_self_excluded.py from
                             davg_change.py and davg_self_excluded.py (cluster); figure by
                             figS13_average_homolog_distance.py
  FigS13/davg_change_summary.json  sizes, medians and ranges of the changes
  The real-versus-AI comparison (formerly Fig. S13, now in the main text; Table S30) uses the files below.
  FigS13/evaluation_distances.parquet  175,537 candidates selected in no arm: pool_position, seq_id, base_3800 (d_min to
                                       the baseline template pool), real_3800 (with the 3,825 pdb3825 sequences added =
                                       the September 2026 PDB template pool), ai_3800 (with the 3,800 candidates
                                       selected on the baseline template pool added), base_690 / real_690 / ai_690 (the
                                       same with the 690-comparison template pool), d_average_homolog
  FigS13/real_vs_ai_metrics.tsv        the comparison statistics by species group (Table S30)
  FigS13/summary.json                  reproduction checks and set overlaps

Fig. S14  UMAP and HDBSCAN maps
  FigS14/hdbscan_on_umap*.tsv  umap1, umap2 and HDBSCAN cluster (-1 = noise) for pdb3825 (A, K97 file), sel3800
                               (B, K109), sel690 optimum (C, K5) and finer view (D, K10), with the protein identifiers.
                               In the rebuttal SI (Fig. S12) only the K97 and K109 maps are shown (A, B)
  FigS14/joint_umap_nn1_pdb3825_sel3800.tsv   rebuttal SI Fig. S12C, D: one UMAP of both sets (embeddings standardized
                               over the 7,625 proteins of pdb3825 and sel3800, PCA50, UMAP n_neighbors 30, min_dist 0.45,
                               random_state 0; plot_dr_diversity_umap_opt.py): cohort, umap1, umap2, nn1_joint_z (Euclidean
                               distance to the nearest other protein of either set in the standardized 1,280-dimensional
                               space), HDBSCAN labels of that script (not shown), seq_id (sel3800) or embed_id (pdb3825)
  FigS14/diversity_metrics_sel3800_pdb3825.tsv, diversity_metric_winners.tsv   the 15 diversity measures of rebuttal
                               SI Section S3.3 for each set and which set is more diverse on each
                               (compare_dr_diversity_3800_pdb3825.py)

Figs. S15-S18  HHsearch and Foldseek top hits before and after additions
  FigS15_S18/{hhsearch,foldseek}_cohort.tsv   per candidate: top hit and its top-hit distance in each of the four
                                              searches (pdbonly = September 2026 PDB template pool; minus3825 = baseline
                                              template pool; plus3800, plus690 = September 2026 PDB template pool with the
                                              selection added); empty = no hit; d_average_homolog (d_avg)
  FigS15_S18/{search}_{set}_pairs.tsv         the candidates of each comparison: d_before, d_after, top_is_added,
                                              top_after, d_average_homolog (S15 hhsearch_sel3800, S16 hhsearch_pdb3825,
                                              S17 foldseek_sel3800, S18 foldseek_pdb3825; the sel690 files back Table S31).
                                              d_average_homolog replaces the former d_homolog (mean distance to all
                                              182,286 plant sequences, itself included; differences of at most 0.012)
  FigS15_S18/summary.tsv, new_coverage.tsv    the values of Tables S31 and S32 and the first-hit counts

Fig. S19  Alignment depth of the compared sets
  FigS19/neff_sets.tsv        14,731 set members (sel3800, sel690, full-length pdb3825, two random sets): AF2-style_Neff and
                              log10_AF2-style_Neff (plotted), neff and log10_neff (empty when the protein has no alignment)
  FigS19/neff_sets_tests.tsv  the values of Table S33 (AF2-style Neff), and the same for Neff

Fig. S20  Gene age of the selections
  FigS20/sel3800_genera_ages.csv, sel690_genera_ages.csv  GenEra phylostratum, rank, taxonomic representativeness,
                                                           age_bin, plddt_global, neff, AF2-style_Neff (A, B),
                                                           frac_disordered
  -> pdb3825 ages and the other set measures: Dataset S15.

Fig. S21  Comparison of the sets (all proteins of each set; before 6 Nov 2026, the four crop species only)
  -> Dataset S15 (Dataset_S15_compared_protein_sets/per_protein_comparison.tsv): one row per protein of each set (18,556).
     Panel I: d_nearest_plant, the pairwise distance to the nearest other plant sequence among the 182,286 unique plant
     sequences of the AFDB v6 set (identical sequences not counted; added 6 Nov 2026 from neighbour_distances.py, with
     its tests in comparison_statistics.tsv); panel J: d_mean_plant, the mean distance to all other plant sequences
     (identical sequences not counted; davg_self_excluded.py). Not the average homolog distance of the selection, which
     is defined only within the candidate pool. Kept from earlier versions: n_within_1.0 / 1.5 / 2.0 (plant sequences
     within these pairwise distances) and d_homolog_mean (the mean including the identical sequence).
