"""Per-locus lineage (tree clade or best BLASTp hit for fragments) and per-genome presence table."""
import pandas as pd
O='/path/to/project/analyses/extended/B_opsin'
P='/path/to/project/scripts'
M=pd.read_csv(O+'/tree/mycobacterial_clades.tsv',sep='\t')
cl2lin={'M1':'A','M3':'A','M4':'A','M5':'A','M2':'B','M6':'C'}
sub={'M1':'A1 (92-seq core, FSM+RGM+SGM)','M3':'A2 (9 RGM)','M4':'A3 (5 FSM)','M5':'A-basal DTG singleton','M2':'B (DTG)','M6':'C (M. madagascariense)'}
loc={}
for _,m in M.iterrows():
    for t in m.members.split(','): loc[t]=(cl2lin[m.clade_id],sub[m.clade_id],'opsin ML tree')
bh=pd.read_csv(O+'/screen/short_fragments_besthits.tsv',sep='\t',header=None,names=['q','s','pid','len','ev','bits']).drop_duplicates('q')
for _,b in bh.iterrows(): loc[b.q]=(loc[b.s][0],loc[b.s][1],f'best BLASTp hit {b.s} ({b.pid:.0f}% over {b.len} aa)')
L=pd.read_csv(O+'/screen/opsin_loci_residues.tsv',sep='\t')
extra=pd.DataFrame([dict(locus='spire_mag_01108952_eggNOGfrag',genome='spire_mag_01108952',aa_len=61),dict(locus='spire_mag_01759821_eggNOGfrag',genome='spire_mag_01759821',aa_len=39)])
s1=pd.read_csv(P+'/inputs/tsv/S1_genome_metadata.tsv',sep='\t'); s1['key']=s1.tip.str.replace(r'\.fa$','',regex=True)
for c in ['clade','species','origin']: extra[c]=extra.genome.map(dict(zip(s1.key,s1[c].str.split().str[0] if c=='clade' else s1[c])))
L=pd.concat([L,extra],ignore_index=True)
L['lineage']=L.locus.map(lambda x:loc[x][0]); L['sublineage']=L.locus.map(lambda x:loc[x][1]); L['lineage_basis']=L.locus.map(lambda x:loc[x][2])
L['in_tree']=L.locus.isin(M.members.str.split(',').explode())
L.to_csv(O+'/screen/opsin_loci_lineages.tsv',sep='\t',index=False)
G=L.groupby('genome').lineage.apply(set).reset_index()
for k in 'ABC': G[k]=G.lineage.apply(lambda s:int(k in s))
G['any']=1
G[['genome','any','A','B','C']].to_csv(O+'/species/genome_lineages.tsv',sep='\t',index=False)
G['clade']=G.genome.map(dict(zip(s1.key,s1.clade.str.split().str[0])))
print(pd.crosstab(G.clade,[G.A,G.B]))
print(G.groupby('clade')[['A','B','C','any']].sum())
# S1 concordance
s1pos=set(s1.key[s1.rhodopsin=='present']); print('S1 positives',len(s1pos),'screen genomes',len(G),'missing from screen',s1pos-set(G.genome),'extra',set(G.genome)-s1pos)
print(pd.crosstab(L.clade,L.lineage))
