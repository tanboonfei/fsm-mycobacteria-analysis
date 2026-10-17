#!/bin/bash
# one OG: MAFFT L-INS-i then trimAl automated1
o=$1
mafft --localpair --maxiterate 1000 --thread 1 --quiet og_fa/$o.faa > aln/$o.aln
trimal -in aln/$o.aln -out trim/$o.trim.faa -automated1 > /dev/null 2>&1
