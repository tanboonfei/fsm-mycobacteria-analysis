# Archive for "A diverse uncultured freshwater clade with reduced genomes in the Mycobacteriaceae"

Selected analysis scripts with their inputs, the bac120 phylogeny, HyPhy RELAX output, the supplementary tables and the habitat (Sandpiper), protist-association and related outputs supporting the manuscript (ISME Communications submission).


All paths below are relative to the repository root.

## Contents

| Folder | Contents | Manuscript item |
|---|---|---|
| `phylogeny/` | bac120 alignment (`bac120.msa.faa`), IQ-TREE log/report, ML tree, consensus tree, rooted tree, splits | species phylogeny, Fig. 1, Table S3 |
| `relax/` | HyPhy RELAX alignments, branch-labelled Newick trees, `*.command.json`, full JSON output (common-reference, joint, synonymous-rate-variation runs in `new_runs/`, superseded trees in `superseded_relax_trees_old/`) | Fig. 3, Table S8 (a-l) |
| `outputs/relax_multistart_stem/` | multi-start RELAX fits (`relax/runs`: JSON, logs), driver scripts, stem-test outputs, result tables | Table S8r, S9r-related convergence tables, Fig. 3a |
| `outputs/relax_controls/` | RELAX control clades, bootstrap, stem tests, alignment cleaning, aBSREL/BUSTED-E, inputs, scripts | Table S8m-q |
| `outputs/habitat_sandpiper/` | source studies, Sandpiper habitat extraction (scripts, logs, `cache/`), light-regime opsin tests | Table S11f-i |
| `outputs/protist/` | protist screen: `A/` (Kraken2 reports, co-occurrence, size-fraction tables, Sandpiper profiles), `B/` (EukProt dataset summary and protein hits, candidate-protein tables, 16S screens, Lake Pavin single-cell read-mapping tables in `pavin_sag/`), `C/` (domain, module and host-gene tables) | Table S15 |
| `accessions/` | GTDB reference accession list (`gtdb_reference_accessions.txt`), SPIRE MAG table and HQ identifiers, ENA run lists for PRJEB54081 (Elbe), PRJEB94020 (Lake Pavin eukaryote-targeted) and PRJEB94149 (Lake Pavin single-cell genomes; `ENA_PRJEB94149_pavin_sag_runs.tsv`), other project run lists and queue metadata | Table S1, S15, Data availability |
| `scripts/` | selected analysis scripts (`README.md` inside), `inputs/` (per-table TSV exports of Additional file 2, rooted tree), `figures/` (figure scripts and their source data), `analyses/` (extended analyses, RELAX multi-start and stem-loss scripts, table recomputations, protist tests, per-genome Bakta provenance, opsin tree, reduced and size-stratified analyses) | Tables S1-S15, Figs. 1-4, S1-S3 |
| `supplementary_tables/` | Additional file 2 (Tables S1-S15, one worksheet per table) | Tables S1-S15 |
| `figures/` | final Figures 1-4 and, in `supplementary/`, Figures S1-S3 (PDF and PNG) | Figs. 1-4, S1-S3 |

Table numbering: `scripts/TABLE_NUMBERING.md` maps the current S1-S15 to the original 28-table numbering still used in file names under `scripts/inputs/tsv/`.

## Genomes and raw data not redistributed

MAGs come from SPIRE (https://spire.embl.de) and reference assemblies from NCBI/GTDB (identifiers in `accessions/` and Table S1). Raw reads (ENA PRJEB54081, PRJEB94020, PRJEB94149) are not redistributed; run accessions are listed in `accessions/`. BLAST/DIAMOND databases, the Kraken2 SILVA database, raw read files and files that exceed practical size were left out; the protist scripts show how they were built.

## Sandpiper extraction

Sandpiper queries were retrieved on 14 October 2026 (from file timestamps of `outputs/habitat_sandpiper/cache/`). Sandpiper's raw-metagenome cutoff is 3 March 2026 and it uses the GTDB R232 release (15 April 2026). Extraction scripts are in `outputs/habitat_sandpiper/scripts/`, logs in `logs/`, and the compressed caches `sandpiper_species_runs.tsv.gz`, `sandpiper_positive_runs.tsv.gz` and `sandpiper_all_runs_clade_abund.tsv.gz` (with fetch logs) in `cache/` and `logs/`.

## Software versions

See `ENVIRONMENT.md`, `environment.yml`, `requirements.txt`, `outputs/habitat_sandpiper/logs/versions.txt`, `outputs/habitat_sandpiper/logs/R_sessionInfo.txt`, `scripts/analyses/extended/B_opsin/VERSIONS.tsv` and `scripts/analyses/extended/F_genes/VERSIONS.txt`.

## Notes

- HyPhy modules directories in the original project were symbolic links; they were not copied (the archive contains no symlinks). Install HyPhy's standard `TemplateBatchFiles` libraries to rerun.
- Some RELAX alignments under `outputs/relax_controls/ci/` were symbolic links to the alignments in `relax/`; they are absent here.
- Local absolute paths in scripts, logs and HyPhy JSON files were replaced by placeholders: `/path/to/project` (project root), `/path/to/databases`, `/path/to/scratch` and `/path/to/conda`. Set these to your own locations before rerunning. Only path strings and internal tree-node labels were changed; no values were altered. The SHA-256 values recorded in Table S8 refer to these cleaned files.
- Left out for size (all regenerable with the scripts): HyPhy `.bf` fit files, minimap2 `.paf` read-mapping files for the Lake Pavin single-cell runs (the classified tables are kept), the whole-gene BUSTED-E JSON (`bustedE_all.json`; the per-stem files are kept), simulated bootstrap alignments (`boot/sim/`) and the full-Bacteria Sandpiper cache file.
- Not every intermediate file of the original working tree is included. The scripts are those behind the reported tables and figures; the genome sequences, sequence databases and raw reads are not redistributed, so a rerun needs them downloaded from the sources above.
