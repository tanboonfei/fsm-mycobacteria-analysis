import pandas as pd, numpy as np
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests
m=pd.read_csv('pavin_abund.tsv',sep='\t')
taxa=['FSM','RGM','SGM','Mycobacteriaceae','acI_Nanopelagicales','Fonsibacter','Polynucleobacter','Limnohabitans','Flavobacterium','Planctomycetota','Legionella','Rickettsiales','Chlamydiota','Holosporales','Megaera','Caedimonas','Paracaedibacterales','Parachlamydiaceae','Midichloriaceae','Rhabdochlamydiaceae']
m['frac']=np.where(m.lower<1,'S','L')  # S=0.65-10 µm, L=10-50 µm
m['pairkey']=m.sample_alias.str[:2]+m.sample_alias.str[3:]   # drop P/G letter
print(m[['sample_alias','frac','pairkey','depth']].drop_duplicates().head(8))
g=m.groupby(['sample_alias','pairkey','frac','depth'])[taxa+['total_cov']].sum().reset_index()
for k in taxa: g[k+'_rel']=g[k]/g.total_cov
w=g.pivot_table(index=['pairkey','depth'],columns='frac',values=[k+'_rel' for k in taxa])
out=[]
for k in taxa:
    x=w[(k+'_rel','L')]; y=w[(k+'_rel','S')]; ok=x.notna()&y.notna()&((x>0)|(y>0))
    if ok.sum()<5: out.append((k,int(ok.sum()),np.nan,np.nan,np.nan)); continue
    lr=np.log2((x[ok]+1e-6)/(y[ok]+1e-6))
    out.append((k,int(ok.sum()),round(float(np.median(lr)),2),round(float((lr>0).mean()),2),wilcoxon(lr).pvalue))
o=pd.DataFrame(out,columns=['taxon','n_pairs','median_log2_L_over_S','frac_pairs_L_higher','wilcoxon_p'])
ok=o.wilcoxon_p.notna(); o.loc[ok,'BH']=multipletests(o.loc[ok,'wilcoxon_p'],method='fdr_bh')[1]
print(o.to_string(index=False)); o.to_csv('pavin_sizefraction_paired.tsv',sep='\t',index=False)
print((g.FSM>0).groupby(g.frac).mean())
g.to_csv('pavin_sample_abund.tsv',sep='\t',index=False)
