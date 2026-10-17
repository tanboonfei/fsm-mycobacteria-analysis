#!/bin/bash
# Stream the first 10M R1 reads of each Lake Pavin protist SAG run (PRJEB94149) and map competitively to 862
# Mycobacteriaceae genomes with minimap2 -x sr -c -N 30 (base-level, secondaries). Keep alignments >=100 bp at >=90% identity.
cd "$(dirname "$0")"
tail -n +2 runs.tsv | tac | head -6 | while IFS=$'\t' read run samp alias title plat strat rc bc fb ftp; do
  [ -e paf/$run.paf ] && continue
  mkdir -p paf
  r1=$(echo $ftp | cut -d';' -f1)
  curl -s --retry 5 "https://$r1" | gzip -dc | head -n 40000000 | \
    minimap2 -x sr -c -N 30 -t 12 myco862.sr.mmi - 2>/dev/null | \
    awk -F'\t' '$11>=100 && $10/$11>=0.90' > paf/$run.paf
  echo "$run $(wc -l < paf/$run.paf)" > paf/$run.paf.done
  echo "$run done $(wc -l < paf/$run.paf)"
done
echo ALLDONE
