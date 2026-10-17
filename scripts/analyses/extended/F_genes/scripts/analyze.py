"""Prevalence tables per group (eggNOG on all 862; KOfam HMM on the 296-genome panel),
FSM salinity contrasts, per-species FSM table."""
from common import *
import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
s1=load_meta(); T=targets()
E=pd.read_csv(SP+'calls_eggnog.tsv',sep='\t',index_col=0)
H=pd.read_csv(SP+'calls_hmm.tsv',sep='\t',index_col=0)
hp=set(H.index)
G=['FSM','FSM_low','FSM_medium','FSM_veryhigh','FSM_salNA','RGM_sister','RGM_other','SGM']
def cnt(df,gene,m):
    m=[x for x in m if x in df.index]
    if not m: return 'NA'
    v=df.loc[m,gene]
    if (v<0).all(): return 'NA'
    return f"{int(v.sum())}/{len(m)}"
def pct(s):
    if s=='NA': return ''
    a,b=map(int,s.split('/')); return f"{100*a/b:.0f}"
rows=[]
for t in T:
    g=t['gene']; r={'point':t['point'],'system':t['system'],'gene':g,'KOs':t['KOs']}
    for grp in G:
        m=members(s1,grp)
        r[f'{grp}_eggNOG']=cnt(E,g,m)
    r['FSM_hc90_eggNOG']=cnt(E,g,members(s1,'FSM',90))
    for grp in ['FSM','FSM_low','FSM_saline','RGM_sister']:
        r[f'{grp}_HMM']=cnt(H,g,[x for x in members(s1,grp) if x in hp])
    r['RGM_other_sub60_HMM']=cnt(H,g,[x for x in members(s1,'RGM_other') if x in hp])
    r['SGM_sub60_HMM']=cnt(H,g,[x for x in members(s1,'SGM') if x in hp])
    # salinity contrast within FSM (eggNOG, union with HMM where available)
    C=E[[g]].copy()
    if g in H.columns and (H[g]>=0).any(): C.loc[H.index,g]=np.maximum(E.loc[H.index,g],H[g])
    lo=members(s1,'FSM_low'); sa=members(s1,'FSM_saline')
    a=int(C.loc[lo,g].sum()); b=int(C.loc[sa,g].sum())
    r['FSM_low_combined']=f"{a}/{len(lo)}"; r['FSM_saline_combined']=f"{b}/{len(sa)}"
    r['fisher_low_vs_saline_p']=fisher_exact([[a,len(lo)-a],[b,len(sa)-b]])[1]
    fs=members(s1,'FSM'); si=members(s1,'RGM_sister')
    a2=int(C.loc[fs,g].sum()); b2=int(C.loc[si,g].sum())
    r['FSM_combined']=f"{a2}/{len(fs)}"; r['RGM_sister_combined']=f"{b2}/{len(si)}"
    r['fisher_FSM_vs_sister_p']=fisher_exact([[a2,len(fs)-a2],[b2,len(si)-b2]])[1]
    r['note']=t['note']
    rows.append(r)
df=pd.DataFrame(rows)
df['BH_low_vs_saline']=multipletests(df.fisher_low_vs_saline_p,method='fdr_bh')[1]
df['BH_FSM_vs_sister']=multipletests(df.fisher_FSM_vs_sister_p,method='fdr_bh')[1]
for c in ['fisher_low_vs_saline_p','BH_low_vs_saline','fisher_FSM_vs_sister_p','BH_FSM_vs_sister']: df[c]=df[c].map(lambda x:f"{x:.2g}")
df.to_csv(F+'TableS8_new_gene_prevalence_by_group.tsv',sep='\t',index=False)
# combined call matrix for FSM per genome (for per-species + per-genome tables)
C=E[[t['gene'] for t in T]].copy()
for t in T:
    g=t['gene']
    if (H[g]>=0).any(): C.loc[H.index,g]=np.maximum(E.loc[H.index,g],H[g])
fsm=s1[s1.group=='FSM'].set_index('gid')
out=C.loc[fsm.index].copy()
out.insert(0,'completeness_pct',fsm.comp); out.insert(0,'salinity_habitat',fsm.sal); out.insert(0,'biome',fsm.biome.fillna('NA')); out.insert(0,'species',fsm.species2); out.insert(0,'origin',fsm.origin)
out.to_csv(F+'FSM_per_genome_calls.tsv',sep='\t')
# per species (>=3 genomes) for informative genes
inf=['betP/betT','proP','mrpA','mrpD','nhaP','kup','kdpB','treS','treY','otsA','mscL','ptsG','oppA','pnuC','amtB','nirB','narK/nasA_MFS','cynS','ureC','feoB','efeU','efeB','irtA','mbtB','bfr','dps','znuA','zupT','sitABCD','thiC','thiD','thiE','thiG','thiM','tenA','panD','panE','egtB','egtD']
sp=out.groupby('species')
rows=[]
for s,dd in sp:
    if len(dd)<3: continue
    r={'species':s,'n':len(dd),'salinity':';'.join(f"{k}:{v}" for k,v in dd.salinity_habitat.value_counts().items()),'median_completeness':round(dd.completeness_pct.median(),1)}
    for g in inf: r[g]=f"{int(dd[g].sum())}/{len(dd)}"
    rows.append(r)
pd.DataFrame(rows).sort_values('n',ascending=False).to_csv(F+'FSM_per_species_selected.tsv',sep='\t',index=False)
print(df[['gene','FSM_eggNOG','FSM_HMM','RGM_sister_eggNOG','RGM_sister_HMM','FSM_low_combined','FSM_saline_combined','BH_low_vs_saline']].to_string())
