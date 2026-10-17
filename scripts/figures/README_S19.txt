S19_discriminating_families.tsv (corrected copy)

Content: the 2,477 clade-discriminating MMseqs2 families (prevalence maximum minus minimum over FSM, RGM and SGM >= 0.5)
from the uniform Prodigal clustering of 822 ingroup genomes (134 FSM, 389 RGM, 299 SGM); 548 FSM-, 1,038 RGM- and 891 SGM-enriched.
Columns: family_rep, enriched_in (clade of highest prevalence), prev_FSM, prev_RGM, prev_SGM (fractions over 134/389/299 genomes),
gene_name, COG_category, description (eggNOG-mapper v2.1.12 on the representative; free text, not curated).

Changes against scripts/inputs/tsv/S19_discriminating_families.tsv: header + 2,477 data rows only. The stray footnote that was appended
as a final row (it broke parsers) is removed from the table and given here. Short rows were padded to 8 columns.

Table numbers in the footnote and in analysis_uniform_clustering/README.txt are in the OLD scheme. Current scheme (scripts/TABLE_NUMBERING.md):
  "Table S18" (FSM-enriched genes)        -> Table S13
  "Table S14" (ancestral stem losses)     -> Table S6h (S14 content), clade summary is now Table S14 (old S2)
  "Table S15" (ancestral state / gene-family loss) -> Table S6
  "Table S19" (this list)                 -> removed from Additional file 2; shipped as this file only
Footnote, corrected: "Prevalence values in this table are computed over all 134 FSM, 389 RGM and 299 SGM genomes and are the values plotted in
Fig. S3. They are not on the same denominator as Table S13, which uses 96 SPIRE FSM MAGs and pools RGM and SGM as the comparator."
Note: the README of analysis_uniform_clustering lists "Tables S19, S18, S14 and S15"; read as "this file, Tables S13, S6h/S14 and S6".
