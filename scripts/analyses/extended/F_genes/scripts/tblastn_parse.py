"""Paralogue-aware parsing of tblastn hits (all queries vs FSM and sister-RGM contigs).
Hits on the same contig and strand that overlap by >=50% of the shorter are one locus; each locus
is assigned to the query with the highest bitscore. Credible: >=40% identity over >=50% of the query
(Table S16k criterion); relaxed: >=30% identity over >=40% (Table S16l criterion)."""
from common import *
import numpy as np
import sys, collections
qlen={}
name=None
for l in open(SP+'tbn/q.faa'):
    if l[0]=='>': name=l[1:].strip(); qlen[name]=0
    else: qlen[name]+=len(l.strip())
def parse(fn):
    H=collections.defaultdict(list)
    for l in open(fn):
        q,s,pid,alen,qs,qe,ss,se,ev,bit=l.split('\t')[:10]
        ss,se=int(ss),int(se); strand='+' if se>ss else '-'
        lo,hi=min(ss,se),max(ss,se)
        H[(s,strand)].append(dict(q=q,pid=float(pid),qcov=100*(int(qe)-int(qs)+1)/qlen[q],lo=lo,hi=hi,bit=float(bit),ev=float(ev),contig=s,strand=strand))
    loci=[]
    for k,hs in H.items():
        hs.sort(key=lambda h:-h['bit']); kept=[]
        for h in hs:
            if any(min(h['hi'],x['hi'])-max(h['lo'],x['lo'])>=0.5*min(h['hi']-h['lo'],x['hi']-x['lo']) for x in kept): continue
            kept.append(h)
        loci+=kept
    return loci
res={}
for tag,fn in (('FSM',SP+'tbn/hits_fsm.tsv'),('RGM_sister',SP+'tbn/hits_sis.tsv')):
    L=parse(fn)
    best=collections.defaultdict(dict)  # genome -> query -> best locus
    for h in L:
        g=h['contig'].split('|')[0]
        if h['q'] not in best[g] or h['bit']>best[g][h['q']]['bit']: best[g][h['q']]=h
    res[tag]=best
pickle.dump(res,open(SP+'tbn/loci.pkl','wb'))
s1=load_meta()
rows=[]
for q in qlen:
    r={'query':q,'query_len_aa':qlen[q]}
    for tag,m in (('FSM',members(s1,'FSM')),('RGM_sister',members(s1,'RGM_sister'))):
        b=res[tag]
        cred=[g for g in m if q in b.get(g,{}) and b[g][q]['pid']>=40 and b[g][q]['qcov']>=50]
        rel=[g for g in m if q in b.get(g,{}) and b[g][q]['pid']>=30 and b[g][q]['qcov']>=40]
        r[f'{tag}_credible']=f"{len(cred)}/{len(m)}"; r[f'{tag}_relaxed']=f"{len(rel)}/{len(m)}"
        if tag=='FSM':
            for sub in ('FSM_low','FSM_saline'):
                mm=members(s1,sub); r[f'{sub}_credible']=f"{sum(g in cred for g in mm)}/{len(mm)}"
            hc=members(s1,'FSM',99); r['FSM_comp>=99_credible']=f"{sum(g in cred for g in hc)}/{len(hc)}"
            neg=[g for g in m if g not in cred]
            r['FSM_negative_median_completeness']=round(float(s1.set_index('gid').loc[neg,'comp'].median()),1) if neg else ''
            bh=max((b[g][q] for g in m if q in b.get(g,{})), key=lambda h:h['bit'], default=None)
            r['FSM_best_hit_id/cov']=f"{bh['pid']:.0f}/{bh['qcov']:.0f}" if bh else 'none'
            pids=[b[g][q]['pid'] for g in cred]; r['FSM_credible_median_id']=round(float(np.median(pids)),1) if pids else ''
    rows.append(r)
pd.DataFrame(rows).to_csv(F+'TableS8_new_tblastn_validation.tsv',sep='\t',index=False)
print(pd.DataFrame(rows).to_string())
