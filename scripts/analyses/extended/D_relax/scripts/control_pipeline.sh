#!/bin/bash
# usage: control_pipeline.sh <clade name> <threads for RELAX>
set -e
D=/path/to/project/analyses/extended/D_relax; C=$1; T=$2
cd $D/controls
if [ ! -s $C.treefile ]; then
  iqtree2 -s $C.codon.fa -m GTR+F+G4 -g $D/inputs/constraint.nwk -T $T --prefix $C -redo -quiet > /dev/null 2>&1
  rm -f $C.ckp.gz $C.mldist $C.bionj $C.uniqueseq.phy $C.model.gz
fi
Rscript $D/scripts/label_control.R $C.treefile $D/inputs/common_reference_FSM_tips_stem.nwk control_clades.json $C $C.lbl.nwk > $C.label.txt
OMP_NUM_THREADS=$T /usr/lib/hyphy/bin/hyphy-avx LIBPATH=/usr/share/hyphy CPU=$T relax --alignment $D/controls/$C.codon.fa --tree $D/controls/$C.lbl.nwk \
   --test Test --reference Reference --code Universal --models Minimal --srv No --output $D/controls/$C.relax.json > $D/controls/$C.relax.log 2>&1
