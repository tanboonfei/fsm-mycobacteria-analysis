import collections
from groups import load
G=load()
grp={g:v['group'] for g,v in G.items()}
N=collections.Counter(grp.values())
mem=collections.defaultdict(lambda: collections.defaultdict(set))
for l in open('../../../mmclust862/c862_cluster.tsv'):
    r,m=l.split(); g=m.split('|')[0]; mem[r][grp[g]].add(g)
def spec(target,minfrac,others_max):
    out=[]
    for r,d in mem.items():
        ft=len(d.get(target,()))/N[target]
        no=sum(len(v) for k,v in d.items() if k!=target); nn=sum(N[k] for k in N if k!=target)
        if ft>=minfrac and no/nn<=others_max: out.append((r,ft,no))
    return out
fs=spec('FSM',0.2,0.01); rs=spec('RGM_sister',0.2,0.01)
print('FSM-specific clusters',len(fs),'sister-specific',len(rs))
open('fsm_specific.ids','w').write('\n'.join(r for r,_,_ in fs)+'\n')
open('sister_specific.ids','w').write('\n'.join(r for r,_,_ in rs)+'\n')
