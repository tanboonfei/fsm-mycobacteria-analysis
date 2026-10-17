"""Pick the Bakta CDS overlapping each tBLASTn opsin locus; write opsin proteins and a locus table."""
import pandas as pd, re
from Bio import SeqIO
O='/path/to/project/analyses/extended/B_opsin'
P='/path/to/project/scripts'
L=pd.read_csv(O+'/screen/opsin_loci.tsv',sep='\t')
off={};clen={};win={}
for r in SeqIO.parse(O+'/synteny/windows.fna','fasta'):
    f=r.description.split(); a,b=map(int,f[3].split('-')); off[r.id]=a; win[r.id]=r.seq; clen[r.id]=int(f[4].split('=')[1])
b=pd.read_csv(O+'/synteny/bakta_windows/windows.tsv',sep='\t',comment='#',header=None,
   names=['seq','type','start','stop','strand','tag','gene','product','xref'])
b=b[b.type=='cds']
faa={r.id:str(r.seq) for r in SeqIO.parse(O+'/synteny/bakta_windows/windows.faa','fasta')}
s1=pd.read_csv(P+'/inputs/tsv/S1_genome_metadata.tsv',sep='\t'); s1['key']=s1.tip.str.replace(r'\.fa$','',regex=True)
cl=dict(zip(s1.key,s1.clade.str.split().str[0])); sp=dict(zip(s1.key,s1.species)); org=dict(zip(s1.key,s1.origin))
rows=[]
for _,r in L.iterrows():
    lo=r.lo-off[r.locus]+1; hi=r.hi-off[r.locus]+1
    c=b[(b.seq==r.locus)&(b.stop>=lo)&(b.start<=hi)].copy()
    c['ov']=c.apply(lambda x:min(x.stop,hi)-max(x.start,lo),axis=1)
    c=c.sort_values('ov',ascending=False)
    if len(c):
        t=c.iloc[0]; seq=faa[t.tag]; tag,cs,ce,st,prod=t.tag,t.start,t.stop,t.strand,t['product']
    else:  # no CDS called: translate the tBLASTn hit span directly (stops -> X)
        nt=win[r.locus][lo-1:hi]; nt=nt[:len(nt)//3*3]
        if r.strand=='-': nt=nt.reverse_complement()
        seq=str(nt.translate()).replace('*','X'); tag,cs,ce,st,prod='none',lo,hi,r.strand,'no CDS called; hit span translated'
    key=r.tip.replace('.fa','')
    rows.append(dict(locus=r.locus,genome=key,clade=cl.get(key),species=sp.get(key),origin=org.get(key),contig=r.contig,
        contig_len=clen[r.locus],locus_start=r.lo,locus_end=r.hi,best_query=r.best_query.split('|')[1],best_bits=r.bits,
        bakta_tag=tag,cds_start_in_window=cs,cds_stop_in_window=ce,cds_strand=st,bakta_product=prod,aa_len=len(seq.rstrip('*')),
        seq=seq.rstrip('*')))
D=pd.DataFrame(rows)
D['dist_to_contig_end']=D.apply(lambda x:min(x.locus_start,x.contig_len-x.locus_end),axis=1)
D.to_csv(O+'/screen/opsin_loci_annotated.tsv',sep='\t',index=False)
with open(O+'/screen/mycobacterial_opsins.faa','w') as f:
    for _,x in D.iterrows(): f.write(f'>{x.locus}\n{x.seq}\n')
print(D.groupby('clade').agg(loci=('locus','size'),genomes=('genome','nunique')))
print(D.aa_len.describe()); print(D[D.aa_len<200][['locus','clade','aa_len','dist_to_contig_end','bakta_product']])
print(D.bakta_product.value_counts())
