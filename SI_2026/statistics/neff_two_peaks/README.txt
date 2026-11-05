Tables behind the explanation of the two peaks of Neff and AF2-style Neff (SI Section S4.1).

Written by two_peaks_server.py (reads the alignments):
  uniclust_sample_coverage.tsv  rows and column coverage of 8,000 random uniclust30_unfiltered alignments
  pdb_rows.tsv                  rows of the three OpenProteinSet alignments of a random 5,000 PDB chains of the original analysis
  pdb_af2style.tsv              AF2-style Neff, Neff, rows and median column coverage of each alignment, for 600 of these chains

Written by neff_two_peaks.py (from these three tables and Datasets S10, S12, S13 and S16; its header describes each table):
  by_assignment.tsv, by_alignment.tsv, rows_per_alignment.tsv, rows_vs_depth.tsv, saturation.tsv, saturation_correlation.tsv,
  homolog_class.tsv, length_bands.tsv, original_data.tsv, same_proteins.tsv, same_chains.tsv, same_chains_rows.tsv
