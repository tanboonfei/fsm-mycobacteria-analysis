"""eggNOG call per genome x target: present if any listed KO or eggNOG preferred name
(base name, numeric _N suffix stripped) is propagated to the genome."""
from common import *
import re, pickle
d=pickle.load(open(SP+'sets.pkl','rb'))
s1=load_meta(); T=targets()
rows=[]
for g in s1.gid:
    gn=set(re.sub(r'_\d+$','',x) for x in d['gn'].get(g,()))
    gk=d['gk'].get(g,set())
    r={'gid':g}
    for t in T:
        ks=[k for k in t['KOs'].split(',') if k!='-']; ns=[n for n in t['eggnog_names'].split(',') if n!='-']
        r[t['gene']]=int(any(k in gk for k in ks) or any(n in gn for n in ns))
        r[t['gene']+'|KO']=int(any(k in gk for k in ks)); r[t['gene']+'|name']=int(any(n in gn for n in ns))
    rows.append(r)
df=pd.DataFrame(rows).set_index('gid'); df.to_csv(SP+'calls_eggnog.tsv',sep='\t')
print(df.shape)
