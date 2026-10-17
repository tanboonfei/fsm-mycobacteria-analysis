# Environment specification

Software used for the analyses reported in the manuscript. Versions in the first two
sections were captured from the machine the analyses were run on. Versions in the third
section are as recorded in the Methods and in the tool outputs themselves (Bakta writes its
version and database into every annotation summary; the per-genome provenance is in
`analyses/bakta_provenance_per_genome.tsv`); they were not reinstalled to re-verify.

Companion files: `environment.yml` (conda) and `requirements.txt` (pip) alongside this file.

## Interpreters and libraries

| Component | Version |
|---|---|
| Python | 3.13.7 |
| pandas | 2.2.3 |
| numpy | 2.4.3 |
| scipy | 1.16.2 |
| statsmodels | 0.14.6 |
| openpyxl | 3.1.5 |
| biopython | 1.83 |
| R | 4.5.2 |
| ape | 5.8.1 |
| phytools | 2.5.2 |
| phylolm | 2.6.5 |
| nlme | 3.1.169 |
| geiger | 2.0.11 |
| castor | 1.8.4 |
| picante | 1.8.2 |

## Command-line tools verified on this machine

| Tool | Version |
|---|---|
| HyPhy | 2.5.59(MP), AVX build; RELAX analysis version 3.1.1 |
| IQ-TREE | 2.0.7 |
| MAFFT | 7.526 |
| BLAST+ | 2.16.0+ |
| Prodigal | 2.6.3 |
| pandoc | 3.1.3 (document rendering only) |

## Tools recorded in the Methods and tool outputs

| Tool | Version |
|---|---|
| Bakta | 1.11.4, database 6.0 (per-genome provenance in `analyses/bakta_provenance_per_genome.tsv`) |
| GTDB-Tk | 2.6.1, release R226 |
| CheckM2 | 0.1.3 (as applied by SPIRE) |
| GUNC | 1.0.1 (as applied by SPIRE) |
| OrthoFinder | 2.5.5 |
| gapseq | 1.4.0 |
| MMseqs2 | 18.8cc5c |
| eggNOG-mapper | 2.1.12, eggNOG database 5.0.2 |
| DIAMOND | 2.1.13 |
| HMMER | 3.4 |
| fastANI | 1.34 |
| SingleM | 0.21.3, metapackage S6.5.0 on GTDB R232 |
| miniprot | 0.18-r281 (orthologue-guided decay screen, Table S5i) |
| IQ-TREE | 2.0.7 (constrained searches, C60 PMSF, AU tests, opsin tree); 3.0.1 (orthologue gene trees and concordance factors) |
| trimAl | 1.5 |
| MAFFT | 7.526 (L-INS-i; --add for RELAX control clades) |
| HyPhy | 2.5.59 (RELAX 3.1.1 with single and multiple starting points, aBSREL, BUSTED-E) |
| l1ou | 1.43 (GitHub khabbazian/l1ou, commit d738756) |
| PhylogeneticEM | 1.8.1 |
| GeneRax | 2.1.3 (likelihoods only; event counts not reported) |
| Infernal | 1.1.5 (TPP riboswitch, Rfam RF00059) |
| logistf / lme4 | 1.26.1 / 1.1.37 |
| Sandpiper | 2.1.1, GTDB R232, tables retrieved 14 October 2026 by file timestamp; Sandpiper raw-metagenome cutoff 3 March 2026; R232 released 15 April 2026 (913,071 runs) |

## Tools invoked inside Bakta 1.11.4

Recorded in every run log; they determine the annotation and are not separately installed.

| Tool | Version |
|---|---|
| Pyrodigal | 3.6.3 |
| tRNAscan-SE | 2.0.12 |
| Aragorn | 1.2.41 |
| Infernal cmscan | 1.1.5 |
| PilerCR | 1.6.0 |
| AMRFinderPlus | 4.0.23 |
| Pyhmmer | 0.11.1 |
| DIAMOND | 2.1.13 |
| blastn | 2.17.0 |

## Notes on reproduction

- `/usr/bin/hyphy` on this system is a `simd-dispatch` wrapper that collapses `"$@"` into a
  single argument and so cannot pass options through. Call the architecture binary directly:
  `/usr/lib/hyphy/bin/hyphy-avx LIBPATH=/usr/share/hyphy relax --alignment ... --tree ...`
- HyPhy RELAX optimisation is stochastic and, on this alignment, single-start fits stop at
  different optima: reruns of the FSM test gave LRT 0.7–12.6, 306–420 log-likelihood units below
  the best fit found. Five-start and warm-started refits (`analyses/relax_stem/relax_scripts`) give
  K = 1.04 and LRT 6.7–9.8 for FSM, still not stably optimised, and K 0.745–0.782 with LRT
  181–185 for the leprosy clade (Table S8r). Use `--starting-points 5` or more when rerunning.
- R's `read.delim` returns `""` rather than `NA` for empty cells. Metadata columns must be
  read with `na.strings = c("NA", "", "NaN", "None")` or missing biome and salinity calls
  enter analyses as absences.
- `castor::asr_mk_model` states must not be thresholded at p > 0.5; about 9% of families are
  unresolved and a 0.5 split assigns all of them to one side.
