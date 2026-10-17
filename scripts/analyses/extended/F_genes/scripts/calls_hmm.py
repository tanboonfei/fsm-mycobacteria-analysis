"""KOfam HMM calls: hmmsearch (HMMER 3.4) of target KOfam profiles vs proteomes of the HMM panel;
a protein is assigned a KO when its score passes the KOfam adaptive threshold (full-sequence score
for score_type 'full', best-domain score for 'domain'), as in KofamScan. Each protein keeps only its
best-scoring passing KO among the target set (to limit paralogue double counting)."""
from common import *
import collections, glob
thr={}
for l in open('/path/to/databases/kofam_2026-10-01/ko_list'):
    f=l.rstrip('\n').split('\t')
    if f[0]=='knum' or f[1]=='-': continue
    thr[f[0]]=(float(f[1]),f[2])
best={}
for fn in sorted(glob.glob(SP+'kofam/part*.tbl')):
    for l in open(fn):
        if l[0]=='#': continue
        f=l.split(); prot,ko=f[0],f[2]; full,dom=float(f[5]),float(f[8])
        if ko not in thr: continue
        t,st=thr[ko]; sc=full if st=='full' else dom
        if sc>=t:
            r=sc/t
            if prot not in best or r>best[prot][1]: best[prot]=(ko,r,sc)
gk=collections.defaultdict(set); gp=collections.defaultdict(list)
for p,(ko,r,sc) in best.items():
    g=p.split('|')[0]; gk[g].add(ko); gp[(g,ko)].append(p)
pickle.dump(dict(gk=dict(gk),gp=dict(gp),best=best),open(SP+'hmm_calls.pkl','wb'))
panel=[l.strip() for l in open(SP+'kofam/hmm_genomes.txt')]
T=targets(); rows=[]
for g in panel:
    r={'gid':g}
    for t in T:
        ks=[k for k in t['KOs'].split(',') if k!='-']
        r[t['gene']]=int(any(k in gk.get(g,()) for k in ks)) if ks else -1
    rows.append(r)
pd.DataFrame(rows).set_index('gid').to_csv(SP+'calls_hmm.tsv',sep='\t')
print(len(best),'proteins with a passing KO;',len(panel),'genomes')
