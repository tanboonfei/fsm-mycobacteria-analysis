"""Assemble the opsin tree input: mycobacterial opsins (>=150 aa) + characterised panel + clustered UniProt PF01036 background."""
import pandas as pd, re
from Bio import SeqIO
O='/path/to/project/analyses/extended/B_opsin'
D=pd.read_csv(O+'/screen/opsin_loci_annotated.tsv',sep='\t')
U=pd.read_csv(O+'/refs/raw/all_pf01036.tsv',sep='\t').set_index('Entry')
useq={r.id.split('|')[1]:str(r.seq) for f in ['all_pf01036.faa','characterised.faa'] for r in SeqIO.parse(f'{O}/refs/raw/{f}','fasta')}
C=pd.read_csv(O+'/refs/characterised.tsv',sep='\t')
def rank(l,r):
    m=re.search(r'([^,]+) \('+r+r'\)',str(l)); return m.group(1).strip() if m else ''
ann=[]; seqs={}
for _,x in D[D.aa_len>=150].iterrows():
    seqs[x.locus]=x.seq.replace('*','').replace('X','')
    ann.append(dict(tip=x.locus,set='mycobacterial',accession='',label=x.species,group=x.clade,organism=x.species,phylum='Actinomycetota',order='Mycobacteriales',family='Mycobacteriaceae'))
for _,c in C.iterrows():
    seqs[c.accession]=useq[c.accession]
    lin=U.loc[c.accession,'Taxonomic lineage'] if c.accession in U.index else ''
    org=U.loc[c.accession,'Organism'] if c.accession in U.index else ''
    ann.append(dict(tip=c.accession,set='characterised',accession=c.accession,label=c.label,group=c['class'],organism=org,phylum=rank(lin,'phylum'),order=rank(lin,'order'),family=rank(lin,'family')))
myco_uniprot={a for a in U.index if 'Mycobacteriaceae (family)' in str(U.loc[a,'Taxonomic lineage'])}
for k in ['actino','bact','arch']:
    for r in SeqIO.parse(f'{O}/refs/bg_{k}_reps.faa','fasta'):
        a=r.id
        if a in seqs or a in myco_uniprot: continue
        seqs[a]=str(r.seq); lin=U.loc[a,'Taxonomic lineage']
        ann.append(dict(tip=a,set='background_'+k,accession=a,label=U.loc[a,'Protein names'].split(' (')[0],group='',organism=U.loc[a,'Organism'],phylum=rank(lin,'phylum'),order=rank(lin,'order'),family=rank(lin,'family')))
A=pd.DataFrame(ann); A.to_csv(O+'/tree/tip_annotation.tsv',sep='\t',index=False)
with open(O+'/tree/opsins_all.faa','w') as f:
    for t in A.tip: f.write(f'>{t}\n{seqs[t]}\n')
print(A.set.value_counts(), len(A))
