#!/bin/bash
# T threads for the concatenated analysis
T=${1:-3}
# 1. concatenated ML tree (IQ-TREE 3), ModelFinder restricted to standard matrices, UFBoot + SH-aLRT
iqtree3 -s concat.faa -m MFP -mset LG,WAG,JTT,Q.pfam -B 1000 -alrt 1000 -T $T --prefix concat -redo > concat_iq.stdout 2>&1
# 2. gene trees (one per OG), ModelFinder per gene with the same matrix set
iqtree3 -S trim_kept -m MFP -mset LG,WAG,JTT,Q.pfam -T $T --prefix genes -redo > genes_iq.stdout 2>&1
# 3. concordance factors on the concatenated tree: gCF, likelihood-based sCF (sCFL, 100 quartets)
M=$(grep "Best-fit model:" concat.log | head -1 | awk '{print $3}')
iqtree3 -te concat.treefile -s concat.faa -m $M --gcf genes.treefile --scfl 100 -T $T --prefix concord -redo > concord.stdout 2>&1
echo T2_DONE
