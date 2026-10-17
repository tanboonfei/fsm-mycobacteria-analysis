#!/bin/bash
# AU test per lineage: ML gene tree (sub-alignment) vs the species-tree topology pruned to the same genomes
O=/path/to/project/analyses/extended/B_opsin/reconciliation
cd $O
for s in A_all A_FSM B_all B_FSM; do
  iqtree2 -s $s.aln -m LG+F+R4 -T 4 --seed 1 --prefix au_${s}_ml -redo -quiet
  cat au_${s}_ml.treefile $s.speciestopo.nwk > au_${s}_trees.nwk
  iqtree2 -s $s.aln -m LG+F+R4 -z au_${s}_trees.nwk -n 0 -zb 10000 -au -T 4 --seed 1 --prefix au_${s}_test -redo -quiet
  echo "== $s"; sed -n '/USER TREES/,/deltaL/p;/^  1 /p;/^  2 /p' au_${s}_test.iqtree | tail -4
done
