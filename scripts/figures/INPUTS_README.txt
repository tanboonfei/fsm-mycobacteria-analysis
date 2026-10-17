inputs/ (everything the builders read except the 307 MB uni824_cluster.tsv used by Fig. S2, set MYCO_UNI_DIR)
S1_genome_metadata.tsv   first 862 rows of Additional_file_2.xlsx sheet S1 (workbook values; M. tilburgii source_material is
                         'cell culture (shell vial)', as in the workbook; the earlier script TSV said 'axenic culture')
myco_bac120.rooted.treefile   identical to scripts/inputs/myco_bac120.rooted.treefile
relax_best_fits.tsv      written by scripts/extract_relax_best_fits.py from analyses/relax_stem/relax/runs (Fig. 3a)
tsv/  exports named by CURRENT table numbers (old number in brackets):
  S4b_decayer_prodigal (S10b), S5h_pseudogenes_MAG_panel (S5b), S5i_per_genome_decay_screen_grouped (extended C_pseudogene results),
  S7_functional_modules (S17; superseded for Fig. 4a by fig4a_values.tsv), S9_codon_usage (S12), S9e_within_pop_pNpS (S13),
  S11d_mag_global_grid_both, S11i_sandpiper_prevalence_abundance_by_habitat (S12i), S13_FSM_enriched_genes (S18),
  S19_discriminating_families (corrected, 2,477 rows)
Optional (produced by other agents): fig4a_values.tsv (Fig. 4a), cog_corrections.tsv (Fig. S3).
Run: python3 scripts/build_figure1_phylogeny.py etc. from the folder above scripts/. Fig. 4 and S3 need no change to run;
set MYCO_FIG_OUT / edit figstyle.SUPP only for preview builds.
