#!/bin/bash
# args: tag model [extra]   builds candidate set from available trees, runs AU (10,000 RELL)
tag=$1; model=$2; extra=$3
: > au_cand_$tag.nwk; : > au_cand_names_$tag.txt
add(){ [ -s "$2" ] && { tail -n1 "$2" >> au_cand_$tag.nwk; echo "$1" >> au_cand_names_$tag.txt; }; }
add ML_LG+G4_r1 r1_lg.treefile
[ -s pmsf_c60.iqtree ] && add ML_PMSF_C60 pmsf_c60.treefile
add ML_PMSF_C20_published r3_pmsf_c20.treefile
add ML_fast30_removed r4_fast30.treefile
for h in $(cat hyps.txt); do add $h cs_$h.treefile; done
iqtree2 -s red.faa -m $model $extra -te r1_lg.treefile -z au_cand_$tag.nwk -n 0 -zb 10000 -au -T 1 --prefix au_$tag -redo -quiet
