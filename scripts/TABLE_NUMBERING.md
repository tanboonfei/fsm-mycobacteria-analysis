# Supplementary tables: old to new numbering

28 tables became 15, numbered in order of first citation in the main text; a sixteenth (protist tests) was added later, and the per-genome Bakta command-line table (old S28, later S2) was then removed, leaving 15 (S1–S15). Merged tables keep their content as lettered sub-tables.

| New | Contents | Old |
|---|---|---|
| S1 | Per-genome metadata | S1 |
| S2 | Binning and fragmentation controls | S4 |
| S3 | Phylogenetically corrected comparisons and robustness; S3f trait associations | S7, plus S8 → S3f |
| S4 | Reduction regimes within and beyond the family | S10, plus S25 → S4c, S25b → S4d |
| S5 | Pseudogene load | S6, plus S5 → S5g, S5b → S5h |
| S6 | Gene-family loss | S15, plus S14 → S6h, S16 → S6i, S16b → S6j |
| S7 | Functional capacity | S17, plus S21a–d → S7b–e, S22 → S7f, S22b → S7g |
| S8 | RELAX | S11 |
| S9 | Codon usage and within-population diversity | S12, plus S13 → S9e, S13b → S9f |
| S10 | Rhodopsin and retinal supply | S20, plus S23 → S10j |
| S11 | Habitat association | S3 |
| S12 | Small genomes outside FSM | S9 |
| S13 | FSM-enriched genes | S18 |
| S14 | Clade summary statistics | S2 |
| S15 | Tests of association with protists | new (cited as S29 in the build sources) |
| removed | Literature-search record | S26; a summary paragraph is now in Additional file 1 |
| removed | Value-to-source map | S27; its main-text citation now points to S2b, and the notes now point to S1 |
| removed | Pairwise ANI | S24; the median within-FSM ANI (80.1%) is stated in the Results, and the pairwise values ship with the analysis code (scripts/inputs/tsv/S24_ani.tsv), as Additional file 1 says |
| removed | Clade-discriminating gene families (Fig. S2) | S19; the list ships with the analysis code as scripts/inputs/tsv/S19_discriminating_families.tsv, and the Fig. S2 and Fig. S3 legends point to it |
| removed | Bakta provenance (per-genome command lines) | S28; removed because every genome used the same Bakta version, database and settings (stated in the Methods); the per-genome log is in scripts/analyses/bakta_provenance_per_genome.tsv |

## Sub-tables added in the eighth revision round

These exist only in the current numbering: S3h–k (placement AU tests, orthologue-tree concordance, unsupervised genome-size shifts), S5i (orthologue-guided decay screen), S7h–k (osmolyte, transport, nitrogen, metal and vitamin genes; transporter substrates and TPP riboswitches), S8m–q (RELAX control clades, stem tests, alignment cleaning, parametric bootstrap, optimiser reproducibility), S10k–l (opsin lineages and loci), S11f–i (source studies, study-controlled habitat tests, opsin by light regime, Sandpiper prevalence by habitat) and S2g (isolate-minus-MAG offsets within species). Table S1 culture_status was re-derived from NCBI records in this round (scripts/analyses/extended/G_culture). In the build sources these appear under the original host numbers (S7h–k, S6i, S17h–k, S11m–q, S20k–l, S3f–i, S4g).
