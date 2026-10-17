import pandas as pd, numpy as np
from k2parse import table
from cooccur import run, add_prot
a=pd.read_csv('pavin_abund.tsv',sep='\t'); a['frac']=np.where(a.lower<1,'S','L')
k=table(a.run_accession.tolist()); m=a.merge(k,on='run_accession')
num=[c for c in k.columns if c!='run_accession']+['FSM','Mycobacteriaceae','acI_Nanopelagicales','Polynucleobacter','Limnohabitans','Planctomycetota','Midichloriaceae','Parachlamydiaceae','Megaera','Legionella','Rickettsiales','total_cov']
g=m.groupby(['sample_alias','frac','depth'])[num].sum().reset_index()
g['nonFSM_Myco']=g.Mycobacteriaceae-g.FSM; g=add_prot(g); g['stratum']=g.frac+'_'+g.depth.astype(str)
print('samples',len(g),g.groupby('stratum').size().to_dict()); print(g.groupby('frac')[['all_protists','Ciliophora','Amoebozoa','Cercozoa','SSU_total','Animalia']].median())
g.to_csv('pavin_cooccur_input.tsv',sep='\t',index=False)
r=run(g,'stratum','Lake Pavin (strata: fraction x depth)')
r.to_csv('pavin_cooccur_results.tsv',sep='\t',index=False)
pd.set_option('display.width',220)
print(r.pivot(index='bacterium',columns='protist_group',values='strat_spearman_rho').round(2).to_string())
print(r[r.bacterium=='FSM'].to_string(index=False))
print(r[r.BH<0.1].to_string(index=False))
