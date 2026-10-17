#!/bin/bash
# tBLASTn opsin screen of all 862 genomes with the Table S11d query panel (as in the manuscript Methods)
O=/path/to/project/analyses/extended/B_opsin
S=/path/to/scratch/gen
mkdir -p $S $O/screen/hits
one() { tip=$1; p=$2; f=$S/$tip.fna
  [ -s $O/screen/hits/$tip.tsv ] && return
  zcat -f $p > $f
  tblastn -query $O/screen/queries_S11d.faa -subject $f -evalue 1e-5 -seg no \
    -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore slen qlen" > $O/screen/hits/$tip.tsv.tmp
  mv $O/screen/hits/$tip.tsv.tmp $O/screen/hits/$tip.tsv
  rm -f $f; }
export -f one; export O S
cut -f1,2 /path/to/project/scripts/analyses/protist/A/genome_paths.tsv | parallel -j6 --colsep '\t' one {1} {2}
