import pandas as pd, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
taxa=[('FSM','FSM'),('nonFSM_Myco','Other Mycobacteriaceae'),('acI_Nanopelagicales','acI (free-living)'),('Polynucleobacter','Polynucleobacter (free-living)'),('Planctomycetota','Planctomycetota (particle)'),('Midichloriaceae','Midichloriaceae (protist symbionts)'),('Parachlamydiaceae','Parachlamydiaceae (amoeba symbionts)'),('Megaera','Megaera (ciliate symbiont)')]
def ratios(g,key,a,b):
    g=g.copy(); g['nonFSM_Myco']=g.Mycobacteriaceae-g.FSM
    out={}
    for k,_ in taxa:
        g[k+'_r']=g[k]/g.total_cov
        w=g.pivot_table(index=key,columns='fr',values=k+'_r'); x,y=w[a],w[b]; ok=x.notna()&y.notna()&((x>0)|(y>0))
        out[k]=np.log2((x[ok]+1e-6)/(y[ok]+1e-6)).values
    return out
e=pd.read_csv('elbe_sample_abund.tsv',sep='\t'); e['fr']=e.fraction; e['Mycobacteriaceae']=e['Mycobacteriaceae']
p=pd.read_csv('pavin_sample_abund.tsv',sep='\t'); p['fr']=p.frac
panels=[('Elbe: suspended particles (>5 µm) vs free-living (0.22–5 µm)',ratios(e,['site','month','rep'],'SPA','FL')),
        ('Lake Pavin 9 m: 10–50 µm vs 0.65–10 µm',ratios(p[p.depth==9],'pairkey','L','S')),
        ('Lake Pavin 80 m: 10–50 µm vs 0.65–10 µm',ratios(p[p.depth==80],'pairkey','L','S'))]
fig,axs=plt.subplots(1,3,figsize=(15,5),sharey=True)
for ax,(t,d) in zip(axs,panels):
    data=[d[k] for k,_ in taxa]
    ax.boxplot(data,vert=False,showfliers=False)
    for i,v in enumerate(data): ax.scatter(v,np.full(len(v),i+1)+np.random.uniform(-.15,.15,len(v)),s=8,alpha=.5,color='C3' if i==0 else 'grey')
    ax.axvline(0,color='k',lw=.8); ax.set_title(t,fontsize=10); ax.set_xlabel('log2 ratio (large / small fraction), paired samples')
axs[0].set_yticks(range(1,len(taxa)+1)); axs[0].set_yticklabels([l for _,l in taxa])
plt.tight_layout(); plt.savefig('sizefraction_log2ratios.png',dpi=150); plt.savefig('sizefraction_log2ratios.pdf')
