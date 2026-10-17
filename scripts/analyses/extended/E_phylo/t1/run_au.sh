#!/bin/bash
# AU tests over the candidate set under three models. Arg: threads
T=${1:-2}
: > au_cand.nwk; : > au_cand_names.txt
add(){ [ -s "$2" ] && { tail -n1 "$2" >> au_cand.nwk; echo "$1" >> au_cand_names.txt; }; }
add ML_LG+G4_bac120 r1_lg.treefile
add ML_LG+F+R7 ml_lgfr7.treefile
add ML_PMSF_C60 pmsf_c60.treefile
add ML_PMSF_C20_published r3_pmsf_c20.treefile
add ML_fast30_removed r4_fast30.treefile
for h in $(cat hyps.txt); do add $h cs_$h.treefile; done
iqtree2 -s red.faa -m LG+G4 -te r1_lg.treefile -z au_cand.nwk -n 0 -zb 10000 -au -T $T --prefix au_lgg4 -redo -quiet
iqtree2 -s red.faa -m LG+F+R7 -te ml_lgfr7.treefile -z au_cand.nwk -n 0 -zb 10000 -au -T $T --prefix au_lgfr7 -redo -quiet
iqtree2 -s red.faa -m LG+C60+F+R4 -fs pmsf_c60.sitefreq -te pmsf_c60.treefile -z au_cand.nwk -n 0 -zb 10000 -au -T $T --prefix au_c60 -redo -quiet
