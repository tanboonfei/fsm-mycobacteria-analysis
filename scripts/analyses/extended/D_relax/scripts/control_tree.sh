#!/bin/bash
D=/path/to/project/analyses/extended/D_relax; C=$1
cd $D/controls
iqtree2 -s $C.codon.fa -m GTR+F+G4 -g ${2:-$D/inputs/constraint.nwk} -T 1 -fast --prefix $C -redo -quiet > /dev/null 2>&1
rm -f $D/controls/$C.ckp.gz $D/controls/$C.parstree $D/controls/$C.mldist $D/controls/$C.bionj $D/controls/$C.uniqueseq.phy
Rscript $D/scripts/label_control.R $C.treefile $D/inputs/common_reference_FSM_tips_stem.nwk control_clades.json $C $C.lbl.nwk > $C.label.txt 2>&1
