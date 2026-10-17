import pandas as pd, json, numpy as np
from k2parse import table
from cooccur import run, add_prot
s=pd.read_csv('elbe_abund.tsv',sep='\t'); s=s[s.fraction!='other']
bs=json.load(open('studies/PRJEB54081_bs.json'))
s['lon']=s.sample_accession.map(lambda a: float(bs.get(a,{}).get('geographic location (longitude)','nan')))
k=table(s.run_accession.tolist())
m=s.merge(k,on='run_accession')
num=[c for c in k.columns if c!='run_accession']+['FSM','Mycobacteriaceae','acI_Nanopelagicales','Polynucleobacter','Limnohabitans','Planctomycetota','Midichloriaceae','Parachlamydiaceae','Megaera','Legionella','Rickettsiales','total_cov']
g=m.groupby(['sample_accession','fraction','month','site']).agg({**{c:'sum' for c in num},'lon':'first'}).reset_index()
g['nonFSM_Myco']=g.Mycobacteriaceae-g.FSM
g=add_prot(g); g['stratum']=g.fraction+'_'+g.month
print('samples',len(g), g.fraction.value_counts().to_dict())
g.to_csv('elbe_cooccur_input.tsv',sep='\t',index=False)
r=run(g,'stratum','Elbe (strata: fraction x month; partial: longitude)',cov_col='lon')
r.to_csv('elbe_cooccur_results.tsv',sep='\t',index=False)
pd.set_option('display.width',200,'display.max_rows',300)
print(r[r.bacterium.isin(['FSM'])].to_string(index=False))
piv=r.pivot(index='bacterium',columns='protist_group',values='strat_spearman_rho'); print(piv.round(2).to_string())
print(r[r.BH<0.05].sort_values('strat_spearman_rho').to_string(index=False))
