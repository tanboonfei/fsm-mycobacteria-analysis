import pandas as pd, numpy as np, sys
from scipy.stats import spearmanr, rankdata
from statsmodels.stats.multitest import multipletests
from k2parse import table
rng=np.random.default_rng(42)
PROT=['Ciliophora','Cercozoa','Amoebozoa','Discoba','Chrysophyceae','Cryptophyceae','Choanoflagellida','Diatomea','Dinoflagellata','het_protists','all_protists']
BACT=['FSM','nonFSM_Myco','acI_Nanopelagicales','Polynucleobacter','Limnohabitans','Planctomycetota','Midichloriaceae','Parachlamydiaceae','Megaera','Legionella','Rickettsiales']
def add_prot(k):
    k['het_protists']=k[['Ciliophora','Cercozoa','Amoebozoa','Discoba','Choanoflagellida','Bicosoecida','Centrohelida','Telonemia','Picozoa']].sum(axis=1)
    k['all_protists']=k[['Amoebozoa','Rhizaria','Alveolata','Stramenopiles','Discoba','Cryptophyceae','Choanoflagellida','Haptophyta','Centrohelida','Telonemia','Picozoa']].sum(axis=1)
    return k
def within_rank(x,strata):
    out=np.empty(len(x))
    for s in np.unique(strata):
        i=strata==s; r=rankdata(x[i]); out[i]=(r-r.mean())/max(r.std(),1e-9)
    return out
def resid(r,c):
    if c is None: return r
    A=np.c_[np.ones(len(c)),c]; beta=np.linalg.lstsq(A,r,rcond=None)[0]; return r-A@beta
def strat_spearman(x,y,strata,nperm=4999,cov=None):
    rc=within_rank(cov,strata) if cov is not None else None
    rx=resid(within_rank(x,strata),rc); ry=resid(within_rank(y,strata),rc)
    rho=np.corrcoef(rx,ry)[0,1]
    cnt=0; idx=[np.where(strata==s)[0] for s in np.unique(strata)]
    for _ in range(nperm):
        ry2=ry.copy()
        for i in idx: ry2[i]=ry[rng.permutation(i)]
        if abs(np.corrcoef(rx,ry2)[0,1])>=abs(rho): cnt+=1
    return rho,(cnt+1)/(nperm+1)
def run(df,strata_col,label,cov_col=None):
    """df has columns: prot counts, SSU_total, bacterial coverages, total_cov"""
    df=df.copy(); res=[]
    strata=df[strata_col].astype(str).values if strata_col else np.array(['all']*len(df))
    for b in BACT:
        if b not in df or (df[b]>0).sum()<8: continue
        pb=df.loc[df[b]>0,b].min()/2
        xb=np.log((df[b]+pb)/df.total_cov)
        for p in PROT:
            if (df[p]>0).sum()<8: continue
            xp=np.log((df[p]+0.5)/df.SSU_total)
            rho,pv=strat_spearman(xb.values,xp.values,strata,cov=(df[cov_col].values.astype(float) if cov_col else None))
            res.append((label,b,p,len(df),int((df[b]>0).sum()),round(rho,3),pv))
    r=pd.DataFrame(res,columns=['dataset','bacterium','protist_group','n_samples','n_bact_present','strat_spearman_rho','perm_p'])
    r['BH']=multipletests(r.perm_p,method='fdr_bh')[1]
    return r
