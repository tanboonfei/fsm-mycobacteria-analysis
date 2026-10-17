#!/bin/bash
# Control: same sp95 AU test on single-copy housekeeping genes (ThrS, TopA; neither is a bac120 marker) for the same genome sets
O=/path/to/project/analyses/extended/B_opsin/reconciliation
C=$O/control; cd $C
python3 - <<'P'
import pandas as pd, glob, os
from Bio import SeqIO
C='/path/to/project/analyses/extended/B_opsin/reconciliation/control'
s1=pd.read_csv('/path/to/project/scripts/inputs/tsv/S1_genome_metadata.tsv',sep='\t'); s1['key']=s1.tip.str.replace(r'\.fa$','',regex=True)
tipof=dict(zip(s1.key,s1.tip)); qlen={'A0QWG2':684,'A0R5D9':936}
rows=[]
for f in glob.glob(C+'/hits/*.tsv'):
    g=os.path.basename(f)[:-4]
    if os.path.getsize(f)==0: continue
    d=pd.read_csv(f,sep='\t',header=None,names='q s pid len qs qe ss se ev bits sseq'.split()); d['g']=g; rows.append(d)
H=pd.concat(rows); H['q']=H.q.str.split('|').str[1]
for q,name in [('A0QWG2','thrS'),('A0R5D9','topA')]:
    b=H[H.q==q].sort_values('bits',ascending=False).drop_duplicates('g')
    b=b[(b.qe-b.qs+1)>=0.6*qlen[q]]
    with open(f'{C}/{name}.faa','w') as fh:
        for _,r in b.iterrows(): fh.write(f'>{tipof[r.g]}\n{r.sseq.replace("-","").replace("*","X")}\n')
    print(name,len(b),'genomes with >=60% query coverage')
P
for gname in thrS topA; do
  mafft --auto --thread 4 $gname.faa > $gname.aln 2>/dev/null
  for s in A_all A_FSM B_all B_FSM; do
    seqkit grep -n -f <(grep '>' $O/$s.aln | sed 's/>//') $gname.aln > ${gname}_$s.raw.aln 2>/dev/null
    trimal -in ${gname}_$s.raw.aln -out ${gname}_$s.aln -gt 0.1 -fasta
    n=$(grep -c '>' ${gname}_$s.aln)
    Rscript -e "suppressPackageStartupMessages(library(ape)); sp<-read.tree('/path/to/project/scripts/inputs/myco_bac120.rooted.treefile'); tips<-names(read.FASTA('${gname}_$s.aln',type='AA')); s<-keep.tip(sp,tips); sup<-suppressWarnings(as.numeric(s\$node.label)); sup[is.na(sup)]<-100; s\$edge.length<-rep(1,nrow(s\$edge)); s\$edge.length[s\$edge[,2] %in% (which(sup<95)+Ntip(s))]<-0; s2<-di2multi(s,tol=1e-8); s2\$node.label<-NULL; s2\$edge.length<-NULL; write.tree(unroot(s2),'${gname}_$s.sp95.nwk')"
    iqtree2 -s ${gname}_$s.aln -m LG+F+R4 -T 4 --seed 1 --prefix ${gname}_${s}_ml -redo -quiet
    iqtree2 -s ${gname}_$s.aln -m LG+F+R4 -g ${gname}_$s.sp95.nwk -T 4 --seed 1 --prefix ${gname}_${s}_con -redo -quiet
    cat ${gname}_${s}_ml.treefile ${gname}_${s}_con.treefile > ${gname}_${s}_trees.nwk
    iqtree2 -s ${gname}_$s.aln -m LG+F+R4 -z ${gname}_${s}_trees.nwk -n 0 -zb 10000 -au -T 4 --seed 1 --prefix ${gname}_${s}_au -redo -quiet
    L=$(python3 -c "from Bio import AlignIO; print(AlignIO.read('${gname}_$s.aln','fasta').get_alignment_length())")
    echo "== $gname $s taxa=$n columns=$L"; sed -n '/^Tree      logL/,/^$/p' ${gname}_${s}_au.iqtree | sed -n '3,4p'
  done
done
