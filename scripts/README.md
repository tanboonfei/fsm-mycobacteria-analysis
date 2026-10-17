# Analysis code

The code behind the reported analyses, with the inputs it reads. This is a selected
set: exploratory work, the annotation pipeline drivers and the data-preparation and
export utilities are not included; the figure builders are in `figures/`. Per-genome Bakta provenance (version, database, run time) is in `analyses/bakta_provenance_per_genome.tsv`; it was Table S2 before that table was removed.

Run everything from this directory. Absolute paths were replaced by placeholders (`/path/to/project` and similar); set them to your own locations first.

| Script | Produces |
|---|---|
| `paths.py` | shared paths and the Table S1 loader; not run directly |
| `binning_artefact_control.py` | Table S2, the isolate-versus-MAG binning control |
| `bakta_isolates_analyse.py` | Table S5, pseudogene load under one annotator |
| `s8_ingroup_models.R` | Table S3f, phylogenetically corrected trait models |
| `s8_gls_corpagel_crosscheck.R` | the GLS cross-check quoted in Table S3f |
| `genome_size_bm_null.R` | parametric bootstrap under Brownian motion for the single-origin genome-size shift (Table S3c) |
| `opsin_simulation_null.R` | Table S11b, the Mk simulation null for opsin presence |

Rerun results are written to `run_output/`, which is created on first use, so that
they can be compared against the workbook without overwriting it.

## Inputs

`inputs/` holds what the scripts read that is not in Additional file 2.

| Input | Read by |
|---|---|
| `tsv/` | all scripts; a per-table text export of `Additional_file_2.xlsx`, one file per table |
| `myco_bac120.rooted.treefile` | the four R scripts |

`bakta_isolates_analyse.py` additionally reads the 260 Bakta output directories and
`analysis_gaps/assembly_stats.tsv` from the project tree, which are too large to
bundle. Its results are in Table S5.

## Reproduction notes

- `inputs/tsv/S1_genome_metadata.tsv` now matches Table S1 of the current Additional file 2 (corrected culture status and provenance, the M. madagascariense opsin call, closure-based completeness). The R scripts and `binning_artefact_control.py` were written against the earlier S1; rerun with the current input, `binning_artefact_control.py` reproduces Table S2 (see `analyses/extended/G_culture/` for the reruns).
  `bakta_isolates_analyse.py` differs from the workbook only in row order and in
  displayed precision.
- HyPhy RELAX optimisation is stochastic; see `ENVIRONMENT.md`.
- Where an older exported table under `inputs/tsv/` disagrees with `Additional_file_2.xlsx`, the workbook holds the final values; the analyses behind the later changes are under `analyses/extended/` and `analyses/relax_stem/`.

## Environment

`../ENVIRONMENT.md` lists interpreter, library and tool versions, and records four
reproduction notes to read before a rerun. `environment.yml` and
`requirements.txt` are alongside it.

## Table numbering in `inputs/tsv/`

The file names in `inputs/tsv/` and the comments inside the scripts follow the original 28-table numbering (and in places older working numbers). `inputs/tsv/` is a per-table export of that original workbook. In the revised Additional file 2 the tables were merged and renumbered to S1–S14 in order of first citation; Table S15 (protist-association tests) was added afterwards, giving S1–S15. `TABLE_NUMBERING.md` maps old numbers to new. Two tables were removed from Additional file 2 and are provided here only:
- `S19_discriminating_families.tsv`: clade-discriminating families, Fig. S2.
- `S24_ani.tsv`: pairwise ANI. The median within-FSM ANI is 80.1% (66 pairs).

## Revision analyses and tree files

- `analyses/`: scripts and outputs of the further analyses:  the phylogenetic robustness and AU tests, size-stratified RELAX, the 178-genome stem losses, the opsin nulls and within-FSM opsin models, the opsin-tree splits, the Blh search, the πN/πS sensitivity panel and the ANI clustering of unassigned FSM genomes. `analyses/RESULTS.md` summarises them; it uses the original table numbering. The scripts were run from a sibling directory of `Publication/`, so their relative input paths assume that layout.
- `analyses/opsin_tree/`: the 135-sequence opsin alignment and IQ-TREE output behind the opsin phylogeny (Table S10h), and the 105-sequence alignment used for the identities in Table S10e (`opsin_identity_105seq.aln`).
- `analyses/protist/`: the tests of association with protists behind Table S15: A/ sample metadata, size fractions and protist co-occurrence; B/ searches of public 16S records, EukProt and Lake Pavin protist single-cell genomes; C/ biosynthetic capacity, iron and host-interaction systems, eukaryotic-like domains and gene donors. `SUMMARY.md` gives the overview; each subdirectory has a its tables. Large databases and downloaded reads are not included.
- `relax_trees/`: the annotated Newick trees and HyPhy RELAX JSON outputs (Table S8d).
- `analyses/extended/`: the supplementary analyses (Tables S2g, S3h–k, S5i, S7h–k, S8m–q, S10k–l, S11f–i and the corrected culture status of Table S1). Each subdirectory (`A_habitat`, `B_opsin`, `C_pseudogene`, `D_relax`, `E_phylo`, `F_genes`, `G_culture`) holds its scripts, result tables and a its tables with methods, software versions and caveats; `figures/` holds the figure scripts used for the current Figures 1–4 and the corrected Table S1 they read. Alignments, per-locus calls, simulated alignments and HyPhy JSON files are too large for the package and are kept in `analyses/extended/`. Some RESULTS.md files name sub-tables by the numbers proposed during analysis; the installed numbers are those listed in `TABLE_NUMBERING.md`.
- `analyses/relax_stem/`: multi-start RELAX refits (Table S8r) and joint stem-loss probabilities (Tables S6l, S6m), with their scripts and tables.

## Figures and final workbook corrections

`figures/` holds the builders for the seven figures (`build_figure1_phylogeny.py` to `build_figureS3_gene_signature.py`, `figstyle.py`, `extract_relax_best_fits.py`), with the small inputs they read in `figures/inputs/` (see `figures/INPUTS_README.txt`; the 307 MB `uni824_cluster.tsv` used by Fig. S2 is set with `MYCO_UNI_DIR`). `analyses/table_recomputations/` records the corrections applied to Additional file 2 after an independent cell-by-cell check: one patch script per group of sheets, the change logs, the findings that could not be resolved, and the Fig. 4a rule and values. The older figure scripts under `analyses/extended/figures/` are superseded.

`bakta_isolates_analyse.py` hard-codes 1,133 pseudogenes for *M. leprae* (PGAP annotation available when the analysis was run; release not recorded). The current RefSeq annotation of GCF_000195855.1 (release 2026-02-06) lists 679.
