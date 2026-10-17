"""Reimplementation of the S6c '503 genomes' negative-binomial rows; old labels must reproduce
IRR 0.682 (p 0.0012) and 0.429 (p 1.3e-11). Then the 4 metagenome-derived genomes are dropped."""
import numpy as np, pandas as pd, statsmodels.formula.api as smf, warnings
warnings.filterwarnings('ignore')
G='/path/to/project/analyses/extended/G_culture/'
T='/path/to/project/scripts/inputs/tsv/'
s=pd.read_csv(T+'S1_genome_metadata.tsv',sep='\t'); pg=pd.read_csv(T+'S5_pseudogenes_all.tsv',sep='\t')
asm=pd.read_csv('/path/to/project/analysis_gaps/assembly_stats.tsv',sep='\t')
c=pd.read_csv(G+'corrected_culture_status.tsv',sep='\t')
d=pg.merge(s[['genome_id_or_accession','tip','genome_size_Mb','n_contigs']],on='genome_id_or_accession').merge(c[['tip','new_status']],on='tip')
d=d.merge(asm[['genome_id_or_accession','n_contigs']].rename(columns={'n_contigs':'n_contigs_asm'}),on='genome_id_or_accession',how='left')
d['logMb']=np.log(d.genome_size_Mb)
res=[]
for lab,dd in [('old (503)',d),('corrected (drop 5 metagenome-derived)',d[~((d.kind=='cultured isolate')&(d.new_status=='GTDB MAG'))])]:
    for clades in (['FSM','RGM','SGM'],['FSM','RGM']):
        x=dd[dd.clade3.isin(clades)].copy(); x['clade']=pd.Categorical(x.clade3,categories=['RGM']+[k for k in clades if k!='RGM'])
        for cc in ['n_contigs','n_contigs_asm']:
            x['cpm']=x[cc]/x.genome_size_Mb
            for f in ['pseudogenes_bakta ~ clade','pseudogenes_bakta ~ clade + cpm']:
                y=x.dropna(subset=['cpm'])
                m=smf.negativebinomial(f,data=y,offset=y.logMb).fit(disp=0,maxiter=200)
                b=m.params['clade[T.FSM]']; res.append(dict(labels=lab,clades='+'.join(clades),contigs=cc,model=f.replace('cpm','contigs/Mb'),n=len(y),
                    beta_FSM=round(b,4),IRR=round(np.exp(b),3),p=f"{m.pvalues['clade[T.FSM]']:.2g}"))
r=pd.DataFrame(res); pd.set_option('display.width',250); print(r.to_string())
r.to_csv(G+'S6c_nb_503_old_vs_new.tsv',sep='\t',index=False)
