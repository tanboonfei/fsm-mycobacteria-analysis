import pandas as pd, numpy as np
from scipy.stats import mannwhitneyu, spearmanr, binom
s1=pd.read_excel('base.xlsx',sheet_name='S1',nrows=862)
for c in ['genome_size_Mb','coding_density_pct','acidic_residues_pct','completeness_pct','contamination_pct']:
    s1[c]=pd.to_numeric(s1[c])
cdu='coding_density_pct_uniform_prodigal'
iso=s1.source_material.isin(['axenic culture','cell culture (shell vial)'])
mag=s1.source_material=='environmental metagenome'
sp=s1.origin=='SPIRE MAG'
print('iso counts',s1[iso].groupby('clade').size().to_dict(), s1[iso].groupby(['clade','source_material']).size().to_dict())
cl={'RGM':s1.clade.str.startswith('RGM'),'SGM':s1.clade=='SGM','FSM':s1.clade.str.startswith('FSM')}
def mw(a,b): return mannwhitneyu(a,b,alternative='two-sided').pvalue
print('--- deposited CD isolate - MAG')
for k in ['RGM','SGM']:
    a=s1[cl[k]&iso].coding_density_pct; 
    for lab,mm in [('allMAG',mag),('SPIRE',sp&mag)]:
        b=s1[cl[k]&mm].coding_density_pct
        print(k,lab,len(a),len(b),round(a.median(),3),round(b.median(),3),round(a.median()-b.median(),3),mw(a,b))
print('--- uniform CD isolate - MAG')
for k in ['RGM','SGM']:
    a=s1[cl[k]&iso][cdu].dropna()
    for lab,mm in [('allMAG',mag),('SPIRE',sp&mag)]:
        b=s1[cl[k]&mm][cdu].dropna()
        print(k,lab,len(a),len(b),round(a.median(),3),round(b.median(),3),round(a.median()-b.median(),3),mw(a,b))
print('--- GTDB vs SPIRE')
gt=mag&(s1.origin=='GTDB reference')
for k in ['SGM','RGM','FSM']:
  for col in ['coding_density_pct',cdu]:
    a=s1[cl[k]&gt][col].dropna();b=s1[cl[k]&sp][col].dropna()
    print(k,col,len(a),len(b),round(a.median(),3),round(b.median(),3),round(a.median()-b.median(),3),mw(a,b))
print('--- spearman FSM SPIRE')
f=s1[cl['FSM']&sp]
for comp in ['completeness_used_pct','completeness_pct','checkm2_this_study_completeness_pct']:
    print(comp,f[comp].notna().sum(),f[comp].median())
    for col in ['genome_size_Mb','GC_pct',cdu]:
        x=f[[comp,col]].apply(pd.to_numeric).dropna()
        r=spearmanr(x[comp],x[col]); print('  ',col,len(x),round(r.statistic,3),r.pvalue)
print('>=99', (f.completeness_pct>=99).sum(), f[f.completeness_pct>=99].genome_size_Mb.median(), f.genome_size_Mb.median())
print('>=99 used', (f.completeness_used_pct>=99).sum(), f[f.completeness_used_pct>=99].genome_size_Mb.median())
print('>=95',(f.completeness_pct>=95).sum(),(f.completeness_used_pct>=95).sum())
ff=s1[cl['FSM']]
print('FSM all134 used>=95',(ff.completeness_used_pct>=95).sum(),(ff.completeness_used_pct>=99).sum())
print('median FSM SPIRE comp_pct',f.completeness_pct.median(),'used',f.completeness_used_pct.median(),' all134 used',ff.completeness_used_pct.median())
n=59
for nn,k in [(59,54),(57,52)]:
    print(nn,k,binom.sf(k-1,nn,0.085),1245*binom.sf(k-1,nn,0.4))
# CDS
g=pd.read_csv('/path/to/project/mag_genomic_stats.tsv',sep='\t')
x=s1.merge(g,left_on='genome_id_or_accession',right_on='genome_id')
for k in ['FSM','RGM']:
    print(k,x[cl[k].reindex(x.index,fill_value=False)].n_genes.mean() if False else x[x.clade.str.startswith(k)].n_genes.agg(['mean','median','count']).to_dict())
print(5787.082474-3358.604167)
