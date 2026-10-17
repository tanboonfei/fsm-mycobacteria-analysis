import collections
def key(t): return t[3:] if t.startswith(('GB_','RS_')) else t.replace('.fa','')
taxa=open('taxa161.txt').read().split(); k2t={key(t):t for t in taxa}
reps=set(open('cand_reps.txt').read().split())
need={}
for l in open('../../../../mmclust862/c862_cluster.tsv'):
    rep,mem=l.rstrip('\n').split('\t')
    if rep in reps:
        g=mem.split('|',1)[0]
        if g in k2t: need[mem]=(rep,k2t[g])
print('need',len(need))
seqs={}
cur=None; buf=[]
def flush():
    if cur in need and cur not in seqs and buf:
        seqs[cur]=''.join(buf).rstrip('*').replace('*','X')
with open('../../../../mmclust862/c862_all_seqs.fasta') as f:
    for l in f:
        if l[0]=='>':
            flush(); cur=l[1:].split()[0]; buf=[]
        else: buf.append(l.strip())
    flush()
print('got',len(seqs))
idx={r:i for i,r in enumerate(sorted(reps))}
out=collections.defaultdict(list)
for m,(rep,t) in need.items():
    if m in seqs: out[rep].append((t,seqs[m]))
with open('og_index.tsv','w') as o:
    for rep in sorted(reps):
        n='og%04d'%idx[rep]; o.write(f'{n}\t{rep}\t{len(out[rep])}\n')
        with open(f'og_fa/{n}.faa','w') as fa:
            for t,s in out[rep]: fa.write(f'>{t}\n{s}\n')
