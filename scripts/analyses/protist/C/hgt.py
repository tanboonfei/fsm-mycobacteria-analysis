import re,collections
T='/path/to/databases/taxdump/'
name2id={}
for l in open(T+'names.dmp'):
    f=l.split('\t|\t'); 
    if f[3].startswith('scientific name'): name2id[f[1]]=f[0]
parent={};rank={}
for l in open(T+'nodes.dmp'):
    f=l.split('\t|\t'); parent[f[0]]=f[1]; rank[f[0]]=f[2]
id2name={v:k for k,v in name2id.items()}
ids_by_name=collections.defaultdict(list)
for l in open(T+'names.dmp'):
    f=l.split('\t|\t')
    if f[3].startswith('scientific name'): ids_by_name[f[1]].append(f[0])
def lineage(org):
    cands=[c for c in [name2id.get(org),name2id.get(' '.join(org.split()[:2]))] if c]
    cands+= ids_by_name.get(org.split()[0],[])
    best=[]
    for t in cands:
        L=[]
        while t!='1' and t in parent: L.append(id2name.get(t,'')); t=parent[t]
        if 'Bacteria' in L: return L
        if not best: best=L
    return best
fsm=set(l.strip() for l in open('fsm_specific.ids')); sis=set(l.strip() for l in open('sister_specific.ids'))
MYC=re.compile('Mycobact|Mycolicibact|Mycobacteroides|Mycolicibacter|Mycobacteriaceae')
best={}
for l in open('spec_nr.tsv'):
    f=l.rstrip('\n').split('\t'); orgs=re.findall(r'\[([^\]]+)\]',f[8]); org=orgs[-1] if orgs else '?'
    if MYC.search(org) or MYC.search(f[8]): continue
    if f[0] not in best: best[f[0]]=(org,float(f[2]))
SYMB=['Legionellales','Chlamydiota','Chlamydiae','Rickettsiales','Holosporales','Coxiellaceae','Francisellaceae','Amoebophilaceae','Candidatus Amoebophilus','Parachlamydiaceae','Midichloriaceae']
def cls(org):
    L=lineage(org)
    if 'Eukaryota' in L: return 'Eukaryota'
    if any(s in L for s in SYMB) or any(s.split()[-1] in org for s in SYMB[-3:]): return 'protist-symbiont lineage'
    if 'Actinomycetota' in L: return 'Actinomycetota (non-Mycobacteriaceae)'
    if 'Bacteria' in L: return 'other Bacteria'
    if 'Archaea' in L: return 'Archaea'
    if 'Viruses' in L: return 'Viruses'
    return 'unclassified'
res={'FSM':collections.Counter(),'sister':collections.Counter()}; idsum={'FSM':[],'sister':[]}
out=open('spec_bestnonmyc.tsv','w'); out.write('cluster\tset\tbest_non_Mycobacteriaceae_org\tpident\tclass\n')
for q in fsm|sis:
    s='FSM' if q in fsm else 'sister'
    if q in best:
        o,p=best[q]; c=cls(o); idsum[s].append(p)
    else: o,p,c='none','', 'no non-Mycobacteriaceae hit'
    res[s][c]+=1; out.write(f'{q}\t{s}\t{o}\t{p}\t{c}\n')
for s in res:
    n=sum(res[s].values()); print(s,n)
    for c,k in res[s].most_common(): print(f'   {c}: {k} ({100*k/n:.1f}%)')
import statistics
for s in idsum: print(s,'median identity to best non-Mycobacteriaceae hit',statistics.median(idsum[s]))

fam=collections.Counter()
for q in fsm:
    if q in best:
        L=lineage(best[q][0]); 
        f=[x for x in L if x.endswith('aceae') or x.endswith('ales')]
        fam[(f[1] if len(f)>1 else '?')+' / '+(f[0] if f else '?')]+=1
print('FSM best-donor order/family top:',fam.most_common(15))
fam=collections.Counter()
for q in sis:
    if q in best:
        L=lineage(best[q][0]); f=[x for x in L if x.endswith('aceae') or x.endswith('ales')]
        fam[(f[1] if len(f)>1 else '?')+' / '+(f[0] if f else '?')]+=1
print('sister best-donor top:',fam.most_common(10))
