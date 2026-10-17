# concatenate trimmed OG alignments over the 161 taxa; write supermatrix + nexus partition + stats
import os
taxa=open('taxa161.txt').read().split(); ogs=open('og_keep.txt').read().split()
sup={t:[] for t in taxa}; parts=[]; pos=1; stats=[]
def rd(fn):
    d={};k=None
    for l in open(fn):
        l=l.strip()
        if l.startswith('>'): k=l[1:].split()[0]; d[k]=[]
        elif k: d[k].append(l)
    return {k:''.join(v) for k,v in d.items()}
for o in ogs:
    fn=f'trim/{o}.trim.faa'
    if not os.path.exists(fn): continue
    a=rd(fn)
    if not a: continue
    L=len(next(iter(a.values())))
    if L<50: stats.append((o,L,len(a),'dropped_short')); continue
    for t in taxa: sup[t].append(a.get(t,'-'*L))
    parts.append((o,pos,pos+L-1)); pos+=L; stats.append((o,L,len(a),'kept'))
with open('concat.faa','w') as f:
    for t in taxa: f.write(f'>{t}\n'+''.join(sup[t])+'\n')
with open('concat.nex','w') as f:
    f.write('#nexus\nbegin sets;\n'+''.join(f'  charset {o} = {a}-{b};\n' for o,a,b in parts)+'end;\n')
with open('og_stats.tsv','w') as f:
    f.write('og\ttrimmed_len\tn_taxa\tstatus\n'+''.join('\t'.join(map(str,s))+'\n' for s in stats))
import statistics
occ={t:1-''.join(sup[t]).count('-')/(pos-1) for t in taxa}
print('genes',len(parts),'sites',pos-1)
S1={}
for i,l in enumerate(open('../t3/S1_genome_metadata.tsv')):
    f=l.rstrip('\n').split('\t')
    if i==0: h=f; continue
    S1[f[0]]=f[h.index('clade')].split()[0]
for g in ('FSM','RGM','SGM','Outside','Outgroup'):
    v=[occ[t] for t in taxa if S1[t]==g]; print(g,len(v),'median occupancy %.3f min %.3f'%(statistics.median(v),min(v)))
