import pandas as pd, numpy as np
from scipy.stats import wilcoxon
m=pd.read_csv('elbe_abund.tsv',sep='\t'); m=m[m.fraction!='other']
taxa=['FSM','RGM','SGM','Mycobacteriaceae','acI_Nanopelagicales','Fonsibacter','Polynucleobacter','Limnohabitans','Flavobacterium','Planctomycetota','Legionella','Rickettsiales','Chlamydiota','Holosporales','Megaera','Caedimonas','Paracaedibacterales','Parachlamydiaceae','Midichloriaceae','Rhabdochlamydiaceae']
# aggregate runs per sample
g=m.groupby(['sample_accession','site','month','rep','fraction'],dropna=False)[taxa+['total_cov']].sum().reset_index()
for k in taxa: g[k+'_rel']=g[k]/g.total_cov
# pair by site, month, rep
key=['site','month','rep']
w=g.pivot_table(index=key,columns='fraction',values=[k+'_rel' for k in taxa])
out=[]
for k in taxa:
    for a,b in [('SPA','FL'),('SINK','FL'),('SINK','SPA')]:
        x=w[(k+'_rel',a)]; y=w[(k+'_rel',b)]
        ok=x.notna()&y.notna()&((x>0)|(y>0))
        if ok.sum()<5: out.append((k,f'{a}/{b}',ok.sum(),np.nan,np.nan,np.nan)); continue
        pc=1e-6
        lr=np.log2((x[ok]+pc)/(y[ok]+pc))
        p=wilcoxon(lr).pvalue
        out.append((k,f'{a}/{b}',int(ok.sum()),round(float(np.median(lr)),2),round(float((lr>0).mean()),2),p))
o=pd.DataFrame(out,columns=['taxon','comparison','n_pairs','median_log2_ratio','frac_pairs_higher_in_first','wilcoxon_p'])
from statsmodels.stats.multitest import multipletests
o['BH']=np.nan; ok=o.wilcoxon_p.notna(); o.loc[ok,'BH']=multipletests(o.loc[ok,'wilcoxon_p'],method='fdr_bh')[1]
o.to_csv('elbe_sizefraction_paired.tsv',sep='\t',index=False)
pd.set_option('display.width',200); print(o.to_string(index=False))
s=g.groupby('fraction')[[k+'_rel' for k in ['FSM','Mycobacteriaceae','acI_Nanopelagicales','Planctomycetota','Flavobacterium']]].median()*100
print(s.round(4)); print(g.groupby('fraction').size())
print('FSM detection by fraction', (g.FSM>0).groupby(g.fraction).agg(['mean','sum','size']))
g.to_csv('elbe_sample_abund.tsv',sep='\t',index=False)
