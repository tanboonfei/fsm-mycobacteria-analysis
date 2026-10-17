#!/bin/bash
# tBLASTn ThrS and TopA (M. smegmatis) against opsin-carrying genomes; keep the best HSP per query/genome
O=/path/to/project/analyses/extended/B_opsin/reconciliation/control
S=/path/to/scratch/gen3
mkdir -p $S $O/hits
one() { tip=$1; p=$2; f=$S/$tip.fna; [ -s $O/hits/$tip.tsv ] && return; zcat -f $p > $f
  tblastn -query $O/control_queries.faa -subject $f -evalue 1e-20 -seg no -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore sseq" > $O/hits/$tip.tsv; rm -f $f; }
export -f one; export O S
parallel -j4 --colsep '\t' one {1} {2} :::: $O/genomes.tsv
