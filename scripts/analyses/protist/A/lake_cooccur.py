import pandas as pd, sys
from k2parse import table
from cooccur import run, add_prot
allres=[]
for s,label in [('PRJEB29238','LakePulse, Canadian lakes (PRJEB29238)'),('PRJEB75477','PELAGICS, European lakes (PRJEB75477)')]:
    a=pd.read_csv(f'{s}_abund.tsv',sep='\t')
    q=pd.read_csv('queue2_meta.tsv',sep='\t'); runs=q[q.study==s].run.tolist()
    k=table(runs); m=a.merge(k,on='run_accession')
    m['nonFSM_Myco']=m.Mycobacteriaceae-m.FSM; m=add_prot(m)
    print(label,'n=',len(m),'FSM present',(m.FSM>0).sum(),'median protist SSU reads',m.all_protists.median(),'ciliates',m.Ciliophora.median())
    m.to_csv(f'{s}_cooccur_input.tsv',sep='\t',index=False)
    allres.append(run(m,None,label))
r=pd.concat(allres); r.to_csv('lake_cooccur_results.tsv',sep='\t',index=False)
pd.set_option('display.width',220,'display.max_rows',400)
for d in r.dataset.unique():
    print(d); print(r[r.dataset==d].pivot(index='bacterium',columns='protist_group',values='strat_spearman_rho').round(2).to_string())
    print(r[(r.dataset==d)&(r.BH<0.05)].sort_values('strat_spearman_rho').to_string(index=False))
