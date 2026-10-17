"""Task 1: independent source studies and samples per clade and genome source set."""
import pandas as pd, numpy as np
O='/path/to/project/analyses/extended/A_habitat/'
d=pd.read_csv(O+'mag_master.tsv',sep='\t',dtype=str)
d['samples']=d.samples.fillna('')
def summ(g):
    samp=set(s for x in g.samples for s in x.split(';') if s)
    vc=g.study.value_counts()
    return pd.Series({'n_genomes':len(g),'n_species_units':g.species_key.nunique(),'n_studies':g.study.nunique(),
        'n_source_samples':len(samp),'largest_study':vc.index[0],'largest_study_n':int(vc.iloc[0]),
        'pct_in_largest_study':round(100*vc.iloc[0]/len(g),1),'n_studies_for_50pct':int(np.searchsorted(vc.cumsum().values,len(g)/2)+1),
        'genomes_per_study_median':float(vc.median())})
rows=[]
for (c,s),g in d.groupby(['clade3','source_set']): rows.append(pd.concat([pd.Series({'clade':c,'genome_set':s}),summ(g)]))
for c,g in d.groupby('clade3'): rows.append(pd.concat([pd.Series({'clade':c,'genome_set':'all MAG-labelled genomes'}),summ(g)]))
for c,g in d[d.source_set!='GTDB genome not flagged metagenome-derived'].groupby('clade3'): rows.append(pd.concat([pd.Series({'clade':c,'genome_set':'true MAGs (SPIRE + NCBI metagenome-derived)'}),summ(g)]))
# annotated subsets used in Fig. 4c / S12b
for lab,mask in [('S1 biome annotated',d.biome.notna()),('S1 salinity annotated',d.salinity_habitat.notna())]:
    for c,g in d[mask].groupby('clade3'): rows.append(pd.concat([pd.Series({'clade':c,'genome_set':lab}),summ(g)]))
for sal in ['low salinity','medium salinity','high salinity','very high salinity']:
    for c,g in d[d.salinity_habitat==sal].groupby('clade3'): rows.append(pd.concat([pd.Series({'clade':c,'genome_set':'salinity = '+sal}),summ(g)]))
t=pd.DataFrame(rows); t.to_csv(O+'T1_studies_per_clade.tsv',sep='\t',index=False)
pd.set_option('display.width',250); print(t.to_string())
# per-study listing (FSM, RGM, SGM)
st=d.groupby(['clade3','study']).agg(n_genomes=('t','size'),n_species=('species_key','nunique'),
    n_spire=('origin',lambda s:(s=='SPIRE MAG').sum()),bioprojects=('bioprojects',lambda s:';'.join(sorted({b for x in s.dropna() for b in x.split(';')}))[:300]),
    ncbi_bp=('ncbi_bioproject',lambda s:';'.join(sorted(set(s.dropna())))),
    biome=('biome',lambda s:','.join(f'{k}:{v}' for k,v in s.value_counts().items())),
    salinity=('salinity_habitat',lambda s:','.join(f'{k}:{v}' for k,v in s.value_counts().items())),
    example_source=('iso',lambda s:str(s.dropna().iloc[0])[:80] if s.notna().any() else ''),
    geo=('geo',lambda s:str(s.dropna().iloc[0])[:60] if s.notna().any() else '')).reset_index()
st.sort_values(['clade3','n_genomes'],ascending=[True,False]).to_csv(O+'T1b_study_list.tsv',sep='\t',index=False)
