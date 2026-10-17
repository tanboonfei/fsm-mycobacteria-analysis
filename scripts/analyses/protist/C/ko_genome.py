import collections
ko={}
for l in open('../../../eggnog_ref_out/allfam.emapper.annotations'):
    if l[0]=='#': continue
    f=l.rstrip('\n').split('\t')
    if f[11] not in ('-',''): ko[f[0]]=[k.replace('ko:','') for k in f[11].split(',')]
print('reps with KO',len(ko))
g2ko=collections.defaultdict(set); nrep=0; hit=0
for l in open('../../../mmclust862/c862_cluster.tsv'):
    r,m=l.split()
    nrep+=1
    if r in ko:
        hit+=1; g2ko[m.split('|')[0]].update(ko[r])
print(nrep,hit,len(g2ko))
with open('genome_ko.tsv','w') as o:
    for g,s in g2ko.items(): o.write(g+'\t'+','.join(sorted(s))+'\n')
