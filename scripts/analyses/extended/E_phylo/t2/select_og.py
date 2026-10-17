import collections, sys
S1 = {}
for i,l in enumerate(open('../t3/S1_genome_metadata.tsv')):
    f=l.rstrip('\n').split('\t')
    if i==0: h=f; continue
    S1[f[0]]=f[h.index('clade')]
def key(t): return t[3:] if t.startswith(('GB_','RS_')) else t.replace('.fa','')
taxa=open('taxa161.txt').read().split()
k2t={key(t):t for t in taxa}
grp={t:('FSM' if S1[t].startswith('FSM') else 'OUT' if S1[t].startswith('Outgroup') else 'OUTSIDE' if S1[t].startswith('Outside') else S1[t].split()[0]) for t in taxa}
cnt=collections.defaultdict(collections.Counter)
for l in open('../../../../mmclust862/c862_cluster.tsv'):
    rep,mem=l.rstrip('\n').split('\t')
    g=mem.split('|',1)[0]
    if g in k2t: cnt[rep][k2t[g]]+=1
rows=[]
for rep,c in cnt.items():
    n=len(c); sc=all(v==1 for v in c.values())
    nf=sum(1 for t in c if grp[t]=='FSM'); no=sum(1 for t in c if grp[t] in('OUT','OUTSIDE'))
    rows.append((rep,n,int(sc),nf,no,sum(c.values())))
with open('cluster_occupancy161.tsv','w') as o:
    o.write('rep\tn_taxa\tsingle_copy\tn_FSM\tn_outgroup8\tn_seqs\n')
    for r in rows: o.write('\t'.join(map(str,r))+'\n')
import itertools
for thr in (0.9,0.85,0.8,0.7):
    s=[r for r in rows if r[2] and r[1]>=thr*161 and r[3]>=thr*40]
    print(thr,len(s), 'with >=2 of 8 outgroup:', sum(1 for r in s if r[4]>=2))
