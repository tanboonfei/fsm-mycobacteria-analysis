#!/bin/bash
# per-gene ML trees, ModelFinder restricted to LG/WAG/JTT/Q.pfam x {+G4} x {model freqs, +F}
mkdir -p gt
for f in trim_kept/*.trim.faa; do o=$(basename $f .trim.faa)
  [ -s gt/$o.treefile ] || iqtree3 -s $f -m MFP -mset LG,WAG,JTT,Q.pfam -mrate G4 -mfreq FU,F -T 1 --prefix gt/$o -redo -quiet > /dev/null 2>&1
done
cat gt/*.treefile > genes.treefile
grep -h "Best-fit model" gt/*.log | awk '{print $3}' | sort | uniq -c > gene_models.txt
echo GT_DONE
