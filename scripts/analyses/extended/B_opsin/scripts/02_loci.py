"""Merge tBLASTn hits into opsin loci and extract +-15 kb windows for uniform Bakta annotation."""
import pandas as pd, gzip, os, sys
from Bio import SeqIO
O='/path/to/project/analyses/extended/B_opsin'
P='/path/to/project/scripts'
gp=dict(l.rstrip('\n').split('\t')[:2] for l in open(P+'/analyses/protist/A/genome_paths.tsv'))
h=pd.read_csv(O+'/screen/tblastn_all_hits.tsv',sep='\t')
h=h[h.ev<=1e-5].copy()
h['lo']=h[['ss','se']].min(axis=1); h['hi']=h[['ss','se']].max(axis=1); h['strand']=(h.se>h.ss).map({True:'+',False:'-'})
loci=[]
for (tip,ctg),g in h.groupby(['tip','s']):
    g=g.sort_values('lo'); cur=None
    for _,r in g.iterrows():
        if cur and r.lo<=cur['hi']+300: cur['hi']=max(cur['hi'],r.hi); cur['bits']=max(cur['bits'],r.bits); cur['qs']=min(cur['qs'],r.qs); cur['qe']=max(cur['qe'],r.qe)
        else:
            if cur: loci.append(cur)
            cur=dict(tip=tip,contig=ctg,lo=r.lo,hi=r.hi,bits=r.bits,strand=r.strand,slen=r.slen,best_query=r.q,qs=r.qs,qe=r.qe)
    loci.append(cur)
L=pd.DataFrame(loci)
L['locus']=[f"{t.replace('.fa','')}_L{i+1}" for t,i in zip(L.tip,L.groupby('tip').cumcount())]
L.to_csv(O+'/screen/opsin_loci.tsv',sep='\t',index=False)
print(len(L),'loci in',L.tip.nunique(),'genomes'); print(L.groupby('tip').size().value_counts())
W=15000
with open(O+'/synteny/windows.fna','w') as out:
    for tip,g in L.groupby('tip'):
        need=set(g.contig)
        with gzip.open(gp[tip.replace('.fa','') if tip.replace('.fa','') in gp else tip],'rt') if gp.get(tip.replace('.fa',''),gp.get(tip,'')).endswith('.gz') else open(gp.get(tip.replace('.fa',''),gp.get(tip))) as fh:
            seqs={r.id:r.seq for r in SeqIO.parse(fh,'fasta') if r.id in need}
        for _,r in g.iterrows():
            s=seqs[r.contig]; a=max(0,r.lo-1-W); b=min(len(s),r.hi+W)
            out.write(f">{r.locus} {tip} {r.contig} {a+1}-{b} ctglen={len(s)}\n{str(s[a:b])}\n")
