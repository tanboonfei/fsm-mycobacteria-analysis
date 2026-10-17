"""Pairwise identity (shared ungapped columns of the L-INS-i alignment) of mycobacterial opsins to characterised references and to UniProt background."""
import pandas as pd, numpy as np
from Bio import AlignIO
O='/path/to/project/analyses/extended/B_opsin'
a={r.id:np.array(list(str(r.seq).upper())) for r in AlignIO.read(O+'/tree/opsins_all.linsi.aln','fasta')}
A=pd.read_csv(O+'/tree/tip_annotation.tsv',sep='\t').set_index('tip')
L=pd.read_csv(O+'/screen/opsin_loci_lineages.tsv',sep='\t'); L=L[L.locus.isin(a)]
L['grp']=L.apply(lambda x:f"{'FSM' if x.clade=='FSM' else 'nonFSM'}-{x.lineage}",axis=1)
def ident(x,y):
    m=(a[x]!='-')&(a[y]!='-'); n=m.sum()
    return (100*(a[x][m]==a[y][m]).mean(), n) if n>=80 else (np.nan,n)
refs=[t for t in A.index if A.loc[t,'set']!='mycobacterial']
rows=[]
for _,x in L.iterrows():
    best=None
    for r in refs:
        i,n=ident(x.locus,r)
        if not np.isnan(i):
            rows.append(dict(locus=x.locus,grp=x.grp,ref=r,ref_set=A.loc[r,'set'],ref_label=A.loc[r,'label'],ref_org=A.loc[r,'organism'],ref_order=A.loc[r,'order'],ref_phylum=A.loc[r,'phylum'],ident=i,cols=n))
I=pd.DataFrame(rows); I.to_csv(O+'/tree/identity_to_refs.tsv.gz',sep='\t',index=False)
ch=I[I.ref_set=='characterised']
print('## median identity to characterised references (top 6 per group)')
for g,gg in ch.groupby('grp'):
    s=gg.groupby(['ref','ref_label']).ident.agg(['median','min','max','size']).sort_values('median',ascending=False).head(6)
    print(g); print(s.round(1).to_string())
print('## best hit among all references per locus: order of best-hit reference')
bh=I.sort_values('ident',ascending=False).drop_duplicates('locus')
print(bh.groupby('grp').apply(lambda d: d.ref_order.value_counts().head(5).to_dict()))
print(bh.groupby('grp').ident.describe().round(1))
bh.to_csv(O+'/tree/best_reference_per_locus.tsv',sep='\t',index=False)
