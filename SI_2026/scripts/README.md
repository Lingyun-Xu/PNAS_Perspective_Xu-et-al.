# Scripts

The Python scripts that built the datasets, drew the figures and computed the statistics of the Supporting
Information, as they were run (on a workstation, an analysis server and a computing cluster). File paths at the top of
each script are placeholders for the folders of the original analysis (`/path/to/project`, `/path/to/home`,
`/path/to/server_data`, `/path/to/cluster`, `/path/to/scratch`); set them to the local folders to run a script. The
scripts of a folder import each other by name, so keep each folder together.

Figure, table and section numbers below are those of the current SI. Some script names and the folders of Dataset S17
keep the numbering of the first SI version (for example `figS13_average_homolog_distance.py` draws Fig. S11 from folder
FigS13); `datasets/Dataset_S17_source_data_SI_figures/README.txt` maps every folder to the current figures.

Some scripts import modules of the original pipelines that are not in this repository: `stats_vm`,
`plot_merged_7groups`, `plot_plddt_density`, `plot_precise_esm`, `plot_published_suite`, `plot_template_distances`,
`select_plants`, `compute_distances`, `plant_library_distances` and `family_metrics`.

## analysis/ (figures and statistics)

| Script | Purpose |
|---|---|
| `si_pnas_server.py` | SI figures drawn on the analysis server (Figs. S2-S4, S8 and S9, and the AFDB v6 group comparison of the first SI version, now main-text Fig. 1a). Needs plot_merged_7groups, plot_plddt_density, plot_precise_esm and plot_published_suite (analysis server pipeline, not included) |
| `si_pnas_cluster.py` | SI figures drawn on the cluster (Figs. S5, S6 and S13-S16). Needs plot_template_distances and select_plants (cluster pipeline, not included) |
| `si_pnas_local.py` | SI figures drawn on the workstation (Figs. S1, S7 composite, S18, S19 and S20, and the earlier four-panel version of Fig. S12) |
| `annotation_pfam_names.py` | InterPro names of the Pfam families of Fig. S20 |
| `plot_neff_sets.py` | Fig. S17 and Table S30, alignment depth of the compared protein sets |
| `fetch_af2style_neff.py` | AF2-style Neff of the OpenProteinSet alignments (per column, median across columns; needs torch and a GPU) |
| `af2style_neff_stats.py` | statistics of AF2-style Neff by group and plant proteome. Needs stats_vm (not included) |
| `two_peaks_server.py` | column coverage of the alignments, and the three OpenProteinSet alignments of the same PDB chains (used by neff_two_peaks.py) |
| `neff_two_peaks.py` | why Neff and AF2-style Neff have two peaks (SI Section S1.4); tables in `statistics/neff_two_peaks/` |
| `maintext_fig1_species_tests.py` | Kruskal-Wallis and Dunn's tests of main-text Fig. 1 (four plant proteomes versus *E. coli*, *S. cerevisiae* and *H. sapiens*: closest-template distance, AF2-style Neff per residue, global pLDDT); table in `statistics/main_text_fig1/` |
| `figS11_S13_candidate_pool.py` | Fig. S10 (candidate pool and selection rounds; option `s10ai` draws the combined figure) and the 3,800-comparison figure of the first SI version, now main-text Fig. 3, from Dataset S17 (folders FigS11-FigS13) |
| `figS13_average_homolog_distance.py` | Fig. S11, change of the average homolog distance when the proteins of sel3800 or sel690 leave the candidate pool (Dataset S17, folder FigS13) |
| `figS12_umap_diversity.py` | Fig. S12: HDBSCAN maps of pdb3825 and sel3800 (A, B) and their joint UMAP coloured by the distance to the nearest other protein (C, D), from Dataset S17 (folder FigS14). Uses redraw_hdbscan_maps.py |
| `redraw_hdbscan_maps.py` | HDBSCAN-on-UMAP maps of Fig. S12 from the saved cluster tables (used by figS12_umap_diversity.py and si_pnas_local.py) |
| `figS15_S18_search_density.py` | Figs. S13-S16 (HHsearch and Foldseek top-hit changes) drawn on the workstation from Dataset S17 (folder FigS15_S18); local version of the cluster code in si_pnas_cluster.py |
| `nolabel_hook.py` | re-exports the figures of a plotting script as label-free panels |
| `exact_p.py` | exact p values of correlation tests that underflow double precision (used by fig1_stats.py and si_pnas_cluster.py) |
| `fig1_stats.py` | statistics of main-text Fig. 1b and 1h-k of the first version, with the code of the 2024 notebook |
| `fig1hk_v6.py` | main-text Fig. 1h-k recomputed on the revision alignments (Dataset S10) |
| `orig_fig1b_why.py` | why the correlation of the original Fig. 1b was weak |
| `fig2d_counts.py` | main-text Fig. 2d counts rebuilt from the organism table |
| `conservation_compare.py` | conservation and gene age, plants versus the other groups |
| `disorder_group_stats.py` | disorder by group and plant proteome. Needs stats_vm (not included) |
| `disorder_template_analysis.py` | disorder versus template distance as explanations of low pLDDT |
| `dtemplate_stats.py` | template-distance statistics (Figs. S5 and S6 and their tables). Needs stats_vm |
| `dtemplate_stats_now.py` | the template-distance tables with the September 2026 PDB. Needs stats_vm |
| `esm_stats_92385.py` | ESMFold statistics on all 92,385 proteins of Dataset S11. Needs stats_vm |
| `esm_corr_92385.py` | ESMFold correlations by group and plant proteome (Dataset S11). Needs stats_vm |
| `plant_plddt_factors.py` | joint model of what determines the pLDDT of plant proteins |
| `plddt_gap_decomposition.py` | decomposition of the pLDDT gap between plants and the other groups |
| `counterfactual_selection.py` | real recent plant PDB structures versus AI-selected plant proteins (cluster; main-text Fig. 3, Table S27). Needs compute_distances, plant_library_distances and select_plants (not included) |
| `davg_self_excluded.py` | average homolog distance without the candidate's own zero distance, for the candidate pool and for the compared sets (cluster). Needs plant_library_distances, select_plants and family_metrics (not included) |
| `davg_change.py` | average homolog distance of the candidate pool after sel3800, sel690 or pdb3825 are added to the template pool (cluster; Fig. S11). Needs plant_library_distances and select_plants (not included) |
| `neighbour_distances.py` | cutoff-free neighbour measures: candidates within each pairwise distance of the added proteins, and the distance from each protein of the compared sets to the nearest other plant sequence (cluster; Section S3.7, Fig. S19). Needs plant_library_distances, select_plants and family_metrics (not included) |

## datasets/ (building the SI datasets)

| Script | Purpose |
|---|---|
| `build_si_datasets.py` | Datasets S1-S3 (Excel) from the result tables |
| `build_s16.py` | Dataset S16, source data of the main-text figures (reads PyTorch tensors; needs torch) |
| `add_af2style_neff.py` | adds AF2-style Neff to Dataset S10 |
| `add_af2style_neff_other_datasets.py` | adds AF2-style Neff to Datasets S2, S3, S13, S14 and S17 (FigS20) |
| `s17_prepare.py` | Dataset S17: the round arrays of Fig. S10D-I (folder FigS12) to TSV, and a summary of every table |
| `s11_s13_source_data.py` | Dataset S17 folders FigS11-FigS13 |
| `apply_davg_self_excluded.py` | applies the corrected average homolog distance (the candidate itself not counted) to Dataset S17 (folders FigS11-FigS13 and FigS15_S18) and writes FigS13/davg_change.tsv.gz |
| `add_nearest_plant.py` | adds the distance to the nearest other plant sequence, and its tests, to Dataset S15 |
| `si_datasets_server.py` | per-protein tables behind Figs. S3 and S9 (analysis server). Needs plot_precise_esm and plot_published_suite |

## checks/ (one-off verifications)

| Script | Purpose |
|---|---|
| `cohort_check.py` | HHsearch and Foldseek cohort counts (Figs. S13-S16) |
| `cohort_check2.py` | Foldseek cohort counts for the selections versus pdb3825 |
| `cohort_check3.py` | Foldseek proteins with and without hits |
| `fig1dg_ref.py` | the Fig. 1d-g plant distances of the first version equal round 0 of its Fig. 3b |
| `fig1hk_source.py` | Jensen-Shannon distances of the original Fig. 1h-k table |
| `fig2_counts_post.py` | tidies the Fig. 2 plotted-count table of Dataset S16 |
| `fig2_counts_post2.py` | second tidy step of the Fig. 2 plotted-count table |
| `fig2b_check.py` | Fig. 2b counts versus the printed percentages |
| `genera_verify.py` | recomputes the gene-age numbers before their use in the SI |
| `genera_verify2.py` | further gene-age checks (published maps, selection composition) |
| `jsd_grid.py` | Jensen-Shannon distance of Fig. 1h-k on the v6 alignments |
| `md_units.py` | units of the mdCATH table (Fig. 2e) |
| `mega_impact.py` | how Fig. 2c changes if all mega-scale rows are counted |
| `nb_stats_cells.py` | lists the statistics cells of the 2024 notebook |
| `orig_fig1b_check.py` | per-protein correlations of the original Fig. 1b table |
| `panel_n.py` | homolog-class counts (Dataset S13) |
| `reg_check.py` | AIUPred regressions of Dataset S12 |
| `s11_s13_harmonize_check.py` | checks before redrawing Figs. S10 and S11 on the selection's definitions |
| `stats_compare.py` | compares recomputed and earlier statistics tables |
| `t11_check.py` | alignment-depth correlations by group (Dataset S10) |
| `t11_check4.py` | Neff and Neff per residue correlations of plant proteins |
| `v4_plant_gap.py` | sequence lengths of the AFDB v4 plant proteins |
| `v4_unique_check.py` | unique sequences per species in AFDB v4 |
