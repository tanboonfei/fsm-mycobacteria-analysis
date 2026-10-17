"""Cluster representatives (mmseqs 862-genome clusters) carrying each target in FSM and sister RGM,
with the number of genomes each covers; used for TCDB/Swiss-Prot substrate hints."""
from common import *
import collections
d=pickle.load(open(SP+'sets.pkl','rb')); s1=load_meta(); T=targets()
grp={g:'FSM' for g in members(s1,'FSM')}; grp.update({g:'RGM_sister' for g in members(s1,'RGM_sister')})
rows=[]
for t in T:
    keys=['K:'+k for k in t['KOs'].split(',') if k!='-']+['N:'+n for n in t['eggnog_names'].split(',') if n!='-']
    cov=collections.defaultdict(lambda:collections.defaultdict(set))
    for (g,k),v in d['prot'].items():
        if k in keys and g in grp:
            for m,r in v: cov[r][grp[g]].add(g)
    for r,c in cov.items():
        rows.append(dict(point=t['point'],gene=t['gene'],rep=r,FSM_genomes=len(c['FSM']),sister_genomes=len(c['RGM_sister']),eggNOG_desc=d['desc'][r][0][:80],PFAMs=d['desc'][r][1],seed=d['desc'][r][2]))
df=pd.DataFrame(rows)
df=df[(df.FSM_genomes>=3)|(df.sister_genomes>=3)].sort_values(['point','gene','FSM_genomes'],ascending=[True,True,False])
df.to_csv(SP+'reps_targets.tsv',sep='\t',index=False); print(len(df), df.rep.nunique())
