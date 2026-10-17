"""Per-genome eggNOG preferred names and KOs for all 862 genomes, propagated from the
862-genome mmseqs cluster representatives (eggnog_ref_out/allfam) to members
(mmclust862/c862_cluster.tsv), exactly as analyses/protist/C (ko_genome.py, hostgenes.py).
Also records protein IDs per (genome, target) for the target table. Read-only on inputs."""
import collections, pickle, re, sys, csv
R='/path/to/project/'
OUT=sys.argv[1]
T=list(csv.DictReader(open(R+'analyses/extended/F_genes/targets.tsv'),delimiter='\t'))
tko=set(k for t in T for k in t['KOs'].split(',') if k!='-')
tnm=set(n for t in T for n in t['eggnog_names'].split(',') if n!='-')
names={}; kos={}; desc={}
for l in open(R+'eggnog_ref_out/allfam.emapper.annotations'):
    if l[0]=='#': continue
    f=l.rstrip('\n').split('\t')
    if f[8] not in ('-',''): names[f[0]]=[x.strip() for x in f[8].split(',')]
    if f[11] not in ('-',''): kos[f[0]]=[k.replace('ko:','') for k in f[11].split(',')]
    desc[f[0]]=(f[7],f[20] if len(f)>20 else '-', f[1])
gn=collections.defaultdict(set); gk=collections.defaultdict(set)
prot=collections.defaultdict(list)   # (genome, key) -> [(member, rep)]
for l in open(R+'mmclust862/c862_cluster.tsv'):
    r,m=l.split(); g=m.split('|')[0]
    for n in names.get(r,()):
        gn[g].add(n); b=re.sub(r'_\d+$','',n)
        if b in tnm or n in tnm: prot[(g,'N:'+b)].append((m,r))
    for k in kos.get(r,()):
        gk[g].add(k)
        if k in tko: prot[(g,'K:'+k)].append((m,r))
pickle.dump(dict(gn=dict(gn),gk=dict(gk),prot=dict(prot),desc={r:desc[r] for (g,k),v in prot.items() for m,r in v}),open(OUT,'wb'))
print(len(gn),len(gk),len(prot))
