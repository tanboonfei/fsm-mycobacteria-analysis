#!/bin/bash
# per-gene ML trees (fast search), ModelFinder restricted to LG/WAG/JTT/Q.pfam x +G4 x {model freqs,+F}. args: shard nshards
s=$1; n=$2; mkdir -p gt; i=0
for f in trim_kept/*.trim.faa; do o=$(basename $f .trim.faa); i=$((i+1))
  [ $(( (i-1) % n )) -eq $s ] || continue
  [ -s gt/$o.treefile ] && grep -q "fast" gt/$o.log 2>/dev/null && continue
  iqtree3 -s $f -m MFP -mset LG,WAG,JTT,Q.pfam -mrate G4 -mfreq FU,F --fast -T 1 --prefix gt/$o -redo -quiet > /dev/null 2>&1
done
echo SHARD_DONE $s
