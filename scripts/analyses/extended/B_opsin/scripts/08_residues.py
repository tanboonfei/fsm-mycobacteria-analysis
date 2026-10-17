"""Read BR positions 85/89/96/216 (P02945 precursor 98/102/109/229) from the L-INS-i alignment; map old-tree lineages."""
import pandas as pd
from Bio import AlignIO
O='/path/to/project/analyses/extended/B_opsin'
a={r.id:str(r.seq) for r in AlignIO.read(O+'/tree/opsins_all.linsi.aln','fasta')}
br=a['P02945']; cols={}; k=0
for i,c in enumerate(br):
    if c!='-':
        k+=1
        if k in (98,102,109,229): cols[k]=i
names={98:'pos85',102:'pos89',109:'pos96',229:'pos216'}
rows=[]
for t,s in a.items():
    d={'tip':t}
    for k,i in cols.items(): d[names[k]]=s[i].upper()
    d['motif']=d['pos85']+d['pos89']+d['pos96']; rows.append(d)
R=pd.DataFrame(rows)
A=pd.read_csv(O+'/tree/tip_annotation.tsv',sep='\t')
R=A.merge(R,on='tip')
L=pd.read_csv(O+'/screen/opsin_loci_annotated.tsv',sep='\t')
old=pd.read_csv(O+'/screen/old_tree_lineages.tsv',sep='\t')
old['genome']=old.old_id.str.split('__').str[0]; old['contig']=old.old_id.str.split('__').str[1].str.rsplit('_',n=1).str[0]
L=L.merge(old[['genome','contig','old_id','old_lineage']],on=['genome','contig'],how='left')
L=L.merge(R[['tip','pos85','pos89','pos96','pos216','motif']].rename(columns={'tip':'locus'}),on='locus',how='left')
L.drop(columns=['seq']).to_csv(O+'/screen/opsin_loci_residues.tsv',sep='\t',index=False)
R.to_csv(O+'/tree/tip_residues.tsv',sep='\t',index=False)
print(pd.crosstab(L.clade,L.motif.fillna('NA')))
print(pd.crosstab(L.old_lineage.fillna('not_in_old_tree'),L.motif.fillna('NA')))
print(L[L.old_lineage.isna()].groupby('clade').size())
