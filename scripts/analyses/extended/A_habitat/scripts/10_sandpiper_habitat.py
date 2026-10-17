"""Task 4: FSM vs RGM prevalence and relative abundance by habitat class from Sandpiper (GTDB R232 SingleM profiles of SRA runs)."""
import pandas as pd, numpy as np, re, os, glob
O='/path/to/project/analyses/extended/A_habitat/'
C=O+'cache/sp_search/'
def hab(e):
    e=str(e).lower().strip()
    if re.search(r'sea cucumber|sea squirt|sea urchin|zooplankton|phytoplankton|niloticus|nilotica',e): return 'host-associated'
    R=[('wastewater/engineered',r'wastewater|sewage|sludge|bioreactor|digester|activated|bioreactor|wwtp|treatment plant|biogas|landfill|leachate'),
       ('sediment',r'\bsediment\b|tidal flat|\bmud\b|mudflat'),
       ('saline/hypersaline lake, salt marsh',r'hypersaline|salt lake|saline|salt marsh|salt pond|saltern|brine|soda lake'),
       ('lake/reservoir',r'\blake\b|lake water|lacustrine|\bpond\b|reservoir|lentic'),
       ('river/stream',r'\briver|riverine|\bstream\b|\bcreek|\blotic'),
       ('estuary/brackish',r'estuar|brackish|lagoon'),
       ('marine',r'marine|seawater|sea water|ocean|\bsea\b|coral reef|hydrothermal vent|sea ice|plankton'),
       ('groundwater/drinking water',r'groundwater|drinking water|tap water|aquifer|water distribution|spring water|subsurface'),
       ('freshwater (unspecified)',r'^freshwater$|^fresh water$|freshwater metagenome'),
       ('aquatic (unspecified)',r'^aquatic$|^water$|^aqueous'),
       ('soil/peat/wetland',r'soil|rhizosphere|\bpeat|wetland|permafrost|terrestrial|biocrust|\bbog\b|\bfen\b|paddy|grassland|forest'),
       ('host-associated',r'human|homo sapiens|\bgut\b|feces|faeces|fecal|oral|skin|mouse|mus musculus|\bpig\b|bovine|chicken|insect|\bfish|invertebrate|coral|sponge|oyster|mollusc|\broot\b|phyllosphere|leaf|rumen|\bmilk|\bbird|tick|\bbee\b|clinical|respiratory|blood|vagina|nasopharyn|\beye\b|saliva|sputum|lung|urinary|plant|arabidopsis|\brat\b|canine|\bdog|\bcat\b|sheep|goat|horse|equus|bos taurus|gallus|primate|macaca|bat\b|rodent|mammal|\bant\b|annelid|nematode|algae|lichen|moss|tracheal|nervous|gill|shrimp|crab|copepod|zooplankton|daphnia|mosquito|termite|tissue|organ|sus scrofa|host'),
      ]
    for k,rx in R:
        if re.search(rx,e): return k
    return 'other/unspecified'
b=pd.read_csv(C+'d__Bacteria_full.csv.gz',usecols=['sample','taxon_name','study_title','eukaryotic_host_association','low_complexity','release_year'])
b=b.rename(columns={'sample':'run'})
b['environment']=b.taxon_name.astype(str).str.replace(r' metagenome$','',regex=True)
b['habitat_label']=b.environment.map(hab)
b['study_key']=b.study_title.astype(str).map(lambda t: re.split(r' - | – ',t)[0].strip())  # JGI titles carry per-sample suffixes after ' - '
def title_hab(t):
    t=str(t).lower()
    for k,rx in [('wastewater/engineered',r'wastewater|sewage|sludge|bioreactor|digester|treatment plant'),('sediment',r'sediment'),
                 ('saline/hypersaline lake, salt marsh',r'hypersaline|salt lake|soda lake|saltern|salt marsh'),
                 ('lake/reservoir',r'\blakes?\b|reservoir|\bponds?\b|lacustrine|limnolog'),('river/stream',r'\brivers?\b|riverine|\bstreams?\b|\bcreek'),
                 ('estuary/brackish',r'estuar|brackish|lagoon'),('marine',r'marine|ocean|\bseas?\b|seawater|coastal|baltic|pacific|atlantic'),
                 ('groundwater/drinking water',r'groundwater|aquifer|drinking water|tap water'),('soil/peat/wetland',r'\bsoils?\b|rhizosphere|peat|wetland|permafrost')]:
        if re.search(rx,t): return k
    return None
amb=b.habitat_label.isin(['freshwater (unspecified)','aquatic (unspecified)','other/unspecified'])
th=b.loc[amb,'study_key'].drop_duplicates(); thm=dict(zip(th,th.map(title_hab)))
b['habitat']=b.habitat_label
b.loc[amb,'habitat']=b.loc[amb,'study_key'].map(thm).fillna(b.loc[amb,'habitat_label'])
b.loc[(b.habitat=='other/unspecified')&(b.eukaryotic_host_association=='host'),'habitat']='host-associated'
print('runs reassigned by study title:',int((b.habitat!=b.habitat_label).sum()))
envmap=b.groupby(['habitat_label','environment']).size().rename('n_runs').reset_index().sort_values(['habitat_label','n_runs'],ascending=[True,False])
envmap.to_csv(O+'T4_env_to_habitat_map.tsv',sep='\t',index=False)
spg=pd.read_csv(O+'../../../Publication/scripts/analyses/protist/A/species_group_r232.tsv',sep='\t')
spg=spg[spg.group.isin(['FSM','RGM','sisterRGM'])].copy(); spg['clade']=spg.group.replace({'sisterRGM':'RGM'})
rows=[]; status=[]
for _,r in spg.iterrows():
    fn=C+r.species.replace(' ','_')+'.csv'
    if not os.path.exists(fn): status.append((r.species,r.clade,'not fetched',0)); continue
    if open(fn).read(5).startswith('#HTTP'): status.append((r.species,r.clade,'no hits/HTTP error',0)); continue
    x=pd.read_csv(fn)
    if len(x)==0: status.append((r.species,r.clade,'0 runs',0)); continue
    x['species']=r.species; x['clade']=r.clade; rows.append(x); status.append((r.species,r.clade,'ok',len(x)))
pd.DataFrame(status,columns=['species','clade','status','n_runs']).to_csv(O+'T4_species_fetch_status.tsv',sep='\t',index=False)
A=pd.concat(rows)
A.to_csv(O+'cache/sandpiper_species_runs.tsv.gz',sep='\t',index=False,compression='gzip')
# run-level clade abundance (sum over species in clade); runs must be in the bacterial denominator table
keep=['sample','study_title','taxon_name','eukaryotic_host_association','metagenome_base_pairs','low_complexity','latitude','longitude','release_year']
runinfo=A.drop_duplicates('sample')[keep].set_index('sample')
ab=A.groupby(['sample','clade']).relative_abundance.sum().unstack(fill_value=0)
nsp=A.groupby(['sample','clade']).species.nunique().unstack(fill_value=0).add_prefix('nsp_')
R=b.set_index('run').join(ab,how='left').join(nsp,how='left')
for c in ['FSM','RGM']: R[c]=R[c].fillna(0.0)
print('runs in species tables but not in bacterial denominator:',len(set(ab.index)-set(b.run)))
R.loc[(R.FSM>0)|(R.RGM>0)].to_csv(O+'cache/sandpiper_positive_runs.tsv.gz',sep='\t',compression='gzip')
def summarise(D,label):
    out=[]
    for h,g in D.groupby('habitat'):
        row={'run_set':label,'habitat':h,'n_runs':len(g)}
        for c in ['FSM','RGM']:
            v=g[c].values; pos=v[v>0]
            row.update({f'{c}_pos_runs':len(pos),f'{c}_prev_pct':round(100*len(pos)/len(v),2),
                        f'{c}_median_pos':round(np.median(pos),3) if len(pos) else np.nan,
                        f'{c}_p90_pos':round(np.percentile(pos,90),3) if len(pos) else np.nan,
                        f'{c}_max':round(v.max(),2),f'{c}_mean_all':round(v.mean(),4),
                        f'{c}_p90_all':round(np.percentile(v,90),3), f'{c}_p99_all':round(np.percentile(v,99),3),
                        f'{c}_pos_studies':g.loc[g[c]>0,'study_key'].nunique(),
                        f'{c}_runs_ge0.1pct':int((v>=0.1).sum()), f'{c}_runs_ge1pct':int((v>=1).sum())})
            st=g.groupby('study_key')[c].agg(['size','median',lambda q:(q>0).mean()])
            st.columns=['n','med','prev']
            row.update({f'{c}_studies_total':len(st), f'{c}_studies_pos':int((st.prev>0).sum()), f'{c}_studies_pos_pct':round(100*(st.prev>0).mean(),1),
                        f'{c}_studies_prev_ge50pct':int((st.prev>=0.5).sum()),
                        f'{c}_study_median_of_medians_all':round(st.med.median(),3), f'{c}_study_p90_of_medians_all':round(st.med.quantile(0.9),3)})
            if len(pos):
                ps=g.loc[g[c]>0].groupby('study_key')[c].median()
                row[f'{c}_median_of_study_medians']=round(ps.median(),3)
                top=g.loc[g[c]>0,'study_key'].value_counts()
                row[f'{c}_top_study_share_pct']=round(100*top.iloc[0]/len(pos),1)
        out.append(row)
    return pd.DataFrame(out)
R2=R.copy(); R2['habitat']=R2.habitat_label
S=pd.concat([summarise(R,'all runs (organism label + study-title refinement)'),summarise(R[~R.low_complexity.astype(bool)],'excluding low-complexity runs'),summarise(R2,'organism label only')])
order=['lake/reservoir','river/stream','estuary/brackish','saline/hypersaline lake, salt marsh','marine','freshwater (unspecified)','aquatic (unspecified)','groundwater/drinking water','sediment','soil/peat/wetland','wastewater/engineered','host-associated','other/unspecified']
S['o']=S.habitat.map({h:i for i,h in enumerate(order)}); S=S.sort_values(['run_set','o']).drop(columns='o')
S.to_csv(O+'T4_sandpiper_habitat_summary.tsv',sep='\t',index=False)
pd.set_option('display.width',260); pd.set_option('display.max_columns',40)
cols=['habitat','n_runs','FSM_pos_runs','FSM_prev_pct','FSM_median_pos','FSM_p90_pos','FSM_max','FSM_pos_studies','FSM_runs_ge1pct','RGM_pos_runs','RGM_prev_pct','RGM_median_pos','RGM_p90_pos','RGM_max','RGM_pos_studies','RGM_runs_ge1pct']
for rs in S.run_set.unique(): print('##',rs); print(S[S.run_set==rs][cols].to_string(index=False))
cols2=['habitat','FSM_studies_total','FSM_studies_pos','FSM_studies_pos_pct','FSM_studies_prev_ge50pct','FSM_median_of_study_medians','FSM_top_study_share_pct','RGM_studies_pos','RGM_studies_pos_pct','RGM_studies_prev_ge50pct','RGM_median_of_study_medians','RGM_top_study_share_pct']
print(S[S.run_set.str.startswith('all runs')][cols2].to_string(index=False))
# share of all FSM / RGM detections by habitat
for c in ['FSM','RGM']:
    g=R[R[c]>0]; print(c,'positive runs',len(g),'studies (normalised titles)',g.study_key.nunique()); print((g.habitat.value_counts(normalize=True)*100).round(1).to_string())
# FSM:RGM ratio of summed abundance per habitat

# FSM-positive runs labelled only 'freshwater metagenome' / 'aquatic metagenome': split by study-title keyword (descriptive only)
def kw(t):
    t=str(t).lower()
    for k,rx in [('lake/reservoir/pond',r'lake|reservoir|pond|lacustrine'),('river/stream',r'river|stream|creek|watershed|catchment'),('estuary/coast',r'estuar|coast|bay|delta'),('aquaculture',r'aquacult|culture water|fish|shrimp|catfish|pond culture'),('drinking/surface water',r'drinking|surface water|tap')]:
        if re.search(rx,t): return k
    return 'unclear'
g=R[(R.FSM>0)&R.habitat.isin(['freshwater (unspecified)','aquatic (unspecified)'])].copy(); g['kw']=g.study_key.map(kw)
k=g.groupby(['habitat','kw']).agg(runs=('FSM','size'),studies=('study_key','nunique'),median_pos=('FSM','median')).reset_index()
k.to_csv(O+'T4b_unspecified_freshwater_by_title_keyword.tsv',sep='\t',index=False); print(k.to_string(index=False))
# top FSM studies overall
t=R[R.FSM>0].groupby(['study_key']).agg(runs=('FSM','size'),habitat=('habitat',lambda s:s.mode().iloc[0]),median=('FSM','median'),max=('FSM','max')).sort_values('runs',ascending=False)
t.head(30).to_csv(O+'T4c_top_FSM_studies.tsv',sep='\t'); print(t.head(15).to_string())
print('FSM: top 10 studies hold %.1f%% of positive runs'%(100*t.runs.head(10).sum()/t.runs.sum()))
R.reset_index()[['run','habitat','habitat_label','study_key','low_complexity','FSM','RGM']].to_csv(O+'cache/sandpiper_all_runs_clade_abund.tsv.gz',sep='\t',index=False,compression='gzip')
