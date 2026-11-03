# Scripts

The Python scripts that built the datasets, drew the figures and computed the statistics of the Supporting
Information, as they were run (on a workstation, an analysis server and a computing cluster). File paths at the top of
each script are placeholders for the folders of the original analysis (`/path/to/project`, `/path/to/home`,
`/path/to/server_data`, `/path/to/cluster`, `/path/to/scratch`); set them to the local folders to run a script. The
scripts of a folder import each other by name, so keep each folder together.

Some scripts import modules of the original pipelines that are not in this repository: `stats_vm`,
`plot_merged_7groups`, `plot_plddt_density`, `plot_precise_esm`, `plot_published_suite`, `plot_template_distances`,
`select_plants`, `compute_distances` and `plant_library_distances`.

## analysis/ (figures and statistics)

| Script | Purpose |
|---|---|
| `si_pnas_server.py` | SI figures drawn on the analysis server (Figs. S1, S3-S5, S9, S10). Needs plot_merged_7groups, plot_plddt_density, plot_precise_esm and plot_published_suite (analysis server pipeline, not included) |
| `si_pnas_cluster.py` | SI figures drawn on the cluster (Figs. S6, S7, S15-S18). Needs plot_template_distances and select_plants (cluster pipeline, not included) |
| `si_pnas_local.py` | SI figures drawn on the workstation (Figs. S2, S8 composite, S14, S20, S21) |
| `figS11_S13_candidate_pool.py` | Figs. S11-S13 from Dataset S17 |
| `redraw_hdbscan_maps.py` | HDBSCAN-on-UMAP maps of Fig. S14 from the saved cluster tables (used by si_pnas_local.py) |
| `nolabel_hook.py` | re-exports the figures of a plotting script as label-free panels |
| `exact_p.py` | exact p values of correlation tests that underflow double precision (used by fig1_stats.py and si_pnas_cluster.py) |
| `fig1_stats.py` | statistics of main-text Fig. 1b and 1h-k, with the code of the 2024 notebook |
| `fig1hk_v6.py` | main-text Fig. 1h-k recomputed on the revision alignments (Dataset S10) |
| `orig_fig1b_why.py` | why the correlation of the original Fig. 1b was weak |
| `fig2d_counts.py` | main-text Fig. 2d counts rebuilt from the organism table |
| `conservation_compare.py` | conservation and gene age, plants versus the other groups |
| `disorder_group_stats.py` | disorder by group and plant proteome. Needs stats_vm (not included) |
| `disorder_template_analysis.py` | disorder versus template distance as explanations of low pLDDT |
| `dtemplate_stats.py` | template-distance statistics (Figs. S6-S7 and their tables). Needs stats_vm |
| `dtemplate_stats_now.py` | the template-distance tables with the September 2026 PDB. Needs stats_vm |
| `esm_stats_92385.py` | ESMFold statistics on all 92,385 proteins of Dataset S11. Needs stats_vm |
| `esm_corr_92385.py` | ESMFold correlations by group and plant proteome (Dataset S11). Needs stats_vm |
| `plant_plddt_factors.py` | joint model of what determines the pLDDT of plant proteins |
| `plddt_gap_decomposition.py` | decomposition of the pLDDT gap between plants and the other groups |
| `counterfactual_selection.py` | real recent plant PDB structures versus AI-selected plant proteins (cluster). Needs compute_distances, plant_library_distances and select_plants (not included) |

## datasets/ (building the SI datasets)

| Script | Purpose |
|---|---|
| `build_si_datasets.py` | Datasets S1-S3 (Excel) from the result tables |
| `build_s16.py` | Dataset S16, source data of the main-text figures (reads PyTorch tensors; needs torch) |
| `s17_prepare.py` | Dataset S17: Fig. S12 round arrays to TSV, and a summary of every table |
| `s11_s13_source_data.py` | Dataset S17 folders FigS11-FigS13 |
| `si_datasets_server.py` | per-protein tables behind Figs. S4 and S10 (analysis server). Needs plot_precise_esm and plot_published_suite |

## checks/ (one-off verifications)

| Script | Purpose |
|---|---|
| `cohort_check.py` | HHsearch and Foldseek cohort counts (Figs. S15-S18) |
| `cohort_check2.py` | Foldseek cohort counts for the selections versus pdb3825 |
| `cohort_check3.py` | Foldseek proteins with and without hits |
| `fig1dg_ref.py` | the Fig. 1d-g plant distances equal round 0 of Fig. 3b |
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
| `s11_s13_harmonize_check.py` | checks before redrawing Figs. S11-S13 on the selection's definitions |
| `stats_compare.py` | compares recomputed and earlier statistics tables |
| `t11_check.py` | alignment-depth correlations by group (Dataset S10) |
| `t11_check4.py` | Neff and Neff per residue correlations of plant proteins |
| `v4_plant_gap.py` | sequence lengths of the AFDB v4 plant proteins |
| `v4_unique_check.py` | unique sequences per species in AFDB v4 |
