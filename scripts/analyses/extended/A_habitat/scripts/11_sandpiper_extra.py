"""Task 4 extras: (a) within-study paired FSM vs RGM prevalence per habitat (sign test);
(b) per-FSM-species habitat profile of Sandpiper detections."""
import pandas as pd, numpy as np
from scipy.stats import binomtest
O='/path/to/project/analyses/extended/A_habitat/'
R=pd.read_csv(O+'cache/sandpiper_all_runs_clade_abund.tsv.gz',sep='\t')
out=[]
for h,g in R.groupby('habitat'):
    st=g.groupby('study_key').agg(n=('run','size'),f=('FSM',lambda v:(v>0).mean()),r=('RGM',lambda v:(v>0).mean()))
    st=st[(st.f>0)|(st.r>0)]
    nf=int((st.f>st.r).sum()); nr=int((st.r>st.f).sum()); tie=int((st.f==st.r).sum())
    p=binomtest(nf,nf+nr).pvalue if nf+nr>0 else np.nan
    out.append(dict(habitat=h,studies_with_either=len(st),FSM_higher=nf,RGM_higher=nr,ties=tie,sign_test_p=p,
                    median_FSM_prev_pct=round(100*st.f.median(),1),median_RGM_prev_pct=round(100*st.r.median(),1)))
P=pd.DataFrame(out); P['sign_test_p']=P.sign_test_p.map(lambda x:f'{x:.2g}')
P.to_csv(O+'T4d_within_study_FSM_vs_RGM_prevalence.tsv',sep='\t',index=False); print(P.to_string(index=False))
A=pd.read_csv(O+'cache/sandpiper_species_runs.tsv.gz',sep='\t',usecols=['sample','species','clade','relative_abundance'])
x=A[A.clade=='FSM'].merge(R[['run','habitat','study_key']],left_on='sample',right_on='run')
fw=['lake/reservoir','river/stream','freshwater (unspecified)']
t=x.groupby('species').agg(runs=('run','size'),studies=('study_key','nunique'),median_pos=('relative_abundance','median'),max=('relative_abundance','max'))
comp=pd.crosstab(x.species,x.habitat,normalize='index').mul(100).round(1)
t['pct_freshwater_lake_river']=comp[[c for c in fw if c in comp]].sum(axis=1).round(1)
t['pct_estuary_marine_saline']=comp[[c for c in ['estuary/brackish','marine','saline/hypersaline lake, salt marsh'] if c in comp]].sum(axis=1).round(1)
t['pct_sediment']=comp.get('sediment',0); t['pct_wastewater']=comp.get('wastewater/engineered',0)
t['top_habitat']=comp.idxmax(axis=1)
t=t.sort_values('runs',ascending=False); t.to_csv(O+'T4e_FSM_species_habitat_profile.tsv',sep='\t')
print(t.to_string()); print(t.top_habitat.value_counts())
