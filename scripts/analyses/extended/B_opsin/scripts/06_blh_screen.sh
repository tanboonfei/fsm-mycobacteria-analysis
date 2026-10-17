#!/bin/bash
# genome-wide tBLASTn for Blh (28 mycobacterial Blh from opsin windows + Q9HNE6 + Q4PNI0) in all 862 genomes
O=/path/to/project/analyses/extended/B_opsin
S=/path/to/scratch/gen2
mkdir -p $S $O/synteny/blh_hits
cat $O/synteny/myco_blh_window.faa $O/synteny/blh_queries.faa > $O/synteny/blh_all_queries.faa
one() { tip=$1; p=$2; f=$S/$tip.fna
  [ -e $O/synteny/blh_hits/$tip.tsv ] && return
  zcat -f $p > $f
  tblastn -query $O/synteny/blh_all_queries.faa -subject $f -evalue 1e-10 -seg no \
    -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore slen qlen" > $O/synteny/blh_hits/$tip.tsv.tmp
  mv $O/synteny/blh_hits/$tip.tsv.tmp $O/synteny/blh_hits/$tip.tsv
  rm -f $f; }
export -f one; export O S
cut -f1,2 /path/to/project/scripts/analyses/protist/A/genome_paths.tsv | parallel -j2 --colsep '\t' one {1} {2}
