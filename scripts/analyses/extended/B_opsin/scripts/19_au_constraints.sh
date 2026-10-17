#!/bin/bash
# Constrained ML searches and AU tests (10,000 RELL) on the per-lineage mycobacterial sub-alignments
O=/path/to/project/analyses/extended/B_opsin/reconciliation
cd $O
for s in A_all A_FSM B_all B_FSM; do
  cons="sp95"; [[ $s == *_all ]] && cons="sp95 FSMmono"
  cp au_${s}_ml.treefile au2_${s}_trees.nwk
  for c in $cons; do
    iqtree2 -s $s.aln -m LG+F+R4 -g $s.constraint_$c.nwk -T 4 --seed 1 --prefix con_${s}_$c -redo -quiet
    cat con_${s}_$c.treefile >> au2_${s}_trees.nwk
  done
  cat $s.speciestopo.nwk >> au2_${s}_trees.nwk
  iqtree2 -s $s.aln -m LG+F+R4 -z au2_${s}_trees.nwk -n 0 -zb 10000 -au -T 4 --seed 1 --prefix au2_${s} -redo -quiet
  echo "== $s trees: ML, $cons, full species topology"; sed -n '/^Tree      logL/,/^$/p' au2_${s}.iqtree
done
