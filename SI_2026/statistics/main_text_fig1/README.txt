Statistics of main-text Fig. 1 (AFDB v6), behind Tables S2, S7 and S12 of SI_PNAS_rebuttal_11062026
====================================================================================================

Copied from updated_maintext_figure/statistics (written by updated_maintext_figure/scripts/maintext_statistics.py and
fig1_panels.py), except fig1_species_kruskal_dunn.tsv (SI_2026/scripts/maintext_fig1_species_tests.py).

fig1a_plddt_by_group.tsv, fig1a_plants_vs_groups.tsv   Fig. 1a: global pLDDT by group; plants vs each group (Table S2)
fig1_kruskal.tsv                                       Kruskal-Wallis across the seven groups (1a; MSA depth; template
                                                       distance)
fig1bc_species_fits.tsv                                least-squares lines and correlations of the species panels:
                                                       rows "1c species *" = main-text Fig. 1b (closest-template
                                                       distance; Table S12), rows "1b AF2-style Neff length *" =
                                                       main-text Fig. 1c (AF2-style Neff per residue; Table S7)
fig1dg_template_distance_*.tsv, jsd_g.tsv              main-text Fig. 1d-g: closest-template distance by species,
                                                       plants vs each species, Jensen-Shannon distance, crossings
                                                       (Table S12)
fig1hk_af2style_neff_length_*.tsv,
jsd_k_af2_style_neff_length.tsv                        main-text Fig. 1h-k: AF2-style Neff per residue (Table S7)
fig1_species_kruskal_dunn.tsv                          Kruskal-Wallis across the four species and Dunn's tests of plants
                                                       vs each species (closest-template distance; AF2-style Neff per
                                                       residue; global pLDDT)
Panel names in the file names are those of the script (1c = closest-template distance, 1b = alignment depth); in the
current main text the closest-template distance is Fig. 1b and the alignment depth Fig. 1c.
