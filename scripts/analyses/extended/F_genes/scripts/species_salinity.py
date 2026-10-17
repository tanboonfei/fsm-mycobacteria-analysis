"""Species-level test (species collapsed to majority presence among their salinity-annotated MAGs):
species with any medium/very-high-salinity MAG vs species with only low-salinity MAGs."""
from common import *
from scipy.stats import fisher_exact
d=pd.read_csv(F+'FSM_per_genome_calls.tsv',sep='\t',index_col=0,keep_default_na=False)
s=d[d.salinity_habitat!='NA']
sp_sal=s.groupby('species').salinity_habitat.apply(lambda x: 'saline' if (x!='low salinity').any() else 'low')
genes=['betP/betT','mrpA','mrpD','bfr','dps','thiM','thiE','oppB','dppB','sitABCD','mntH','kup','treS','treY','feoB','efeU','znuA','zupT','nthA','ectD','mmpL4/mmpS4','thiXYZ (HMP ABC)']
rows=[]
for g in genes:
    maj=s.groupby('species')[g].mean()>=0.5
    a=int(maj[sp_sal=='saline'].sum()); na=int((sp_sal=='saline').sum()); b=int(maj[sp_sal=='low'].sum()); nb=int((sp_sal=='low').sum())
    allfsm=d[g]; nas=d[d.salinity_habitat=='NA']
    rows.append(dict(gene=g,saline_species_with_gene=f"{a}/{na}",low_species_with_gene=f"{b}/{nb}",species_fisher_p=f"{fisher_exact([[a,na-a],[b,nb-b]])[1]:.2g}",
        saline_MAGs=f"{int(s[s.salinity_habitat!='low salinity'][g].sum())}/{int((s.salinity_habitat!='low salinity').sum())}",low_MAGs=f"{int(s[s.salinity_habitat=='low salinity'][g].sum())}/{int((s.salinity_habitat=='low salinity').sum())}",
        salinityNA_genomes=f"{int(nas[g].sum())}/{len(nas)}",
        saline_species_list=';'.join(f"{k}:{int(s[(s.species==k)&(s.salinity_habitat!='low salinity')][g].sum())}/{int(((s.species==k)&(s.salinity_habitat!='low salinity')).sum())}" for k in sp_sal[sp_sal=='saline'].index)))
df=pd.DataFrame(rows); df.to_csv(F+'TableS8_new_salinity_species_level.tsv',sep='\t',index=False)
print(sp_sal.value_counts().to_dict()); print(df.to_string())
