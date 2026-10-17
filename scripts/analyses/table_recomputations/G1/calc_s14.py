import pandas as pd, numpy as np
from scipy.stats import mannwhitneyu
s1=pd.read_excel('base.xlsx',sheet_name='S1',nrows=862)
for c in ['genome_size_Mb','coding_density_pct','acidic_residues_pct','completeness_pct']: s1[c]=pd.to_numeric(s1[c])
g=pd.read_excel('base.xlsx',sheet_name='S5',header=None,skiprows=58,nrows=503).iloc[:,:5]
g.columns=['acc','clade3','kind','sp','pg']; g['pg']=pd.to_numeric(g.pg)
sp=s1[s1.origin=='SPIRE MAG'].merge(g[g.kind=='MAG'][['acc','pg']],left_on='genome_id_or_accession',right_on='acc',how='left')
sp['c3']=np.where(sp.clade.str.startswith('FSM'),'FSM',np.where(sp.clade.str.startswith('RGM'),'RGM',np.where(sp.clade=='SGM','SGM','other')))
print(sp.c3.value_counts().to_dict(), sp.pg.notna().groupby(sp.c3).sum().to_dict())
cd='coding_density_pct_uniform_prodigal'
for met in ['genome_size_Mb','GC_pct','coding_density_pct',cd,'acidic_residues_pct','pg']:
    for a,b in [('FSM','RGM'),('FSM','SGM'),('RGM','SGM')]:
        x=sp[sp.c3==a][met].dropna();y=sp[sp.c3==b][met].dropna()
        print(met,a,b,len(x),len(y),round(x.median(),3),round(y.median(),3),'%.2e'%mannwhitneyu(x,y).pvalue)
# all-genome medians
s1['c3']=s1.clade.map(lambda c:'FSM' if c.startswith('FSM') else 'RGM' if c.startswith('RGM') else 'SGM' if c=='SGM' else 'Outgroup' if c.startswith('Outgroup') else 'Outside')
for c,gg in s1.groupby('c3'):
    print(c,len(gg),gg.genome_size_Mb.median(),gg.GC_pct.median(),gg.coding_density_pct.median(),gg.acidic_residues_pct.median(),gg[cd].median(),gg[cd].notna().sum(),
      [round(100*(gg[k]=='present').mean(),1) for k in ['rhodopsin','urease','ectoine_synthase']])
print(s1.rhodopsin.unique(),s1.urease.unique())
