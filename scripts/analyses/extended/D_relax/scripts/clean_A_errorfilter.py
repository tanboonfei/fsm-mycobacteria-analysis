# Re-implementation of HyPhy error-filter.bf (v0.1; the shipped file fails to parse under HyPhy 2.5.59) on a BUSTED-E json.
# Same rule and defaults: error-sink empirical BF >= 100 vs prior odds AND BF(error vs fastest omega class) >= 20;
# terminal branch -> mask that sequence's codon; internal branch -> mask all leaf descendants (complement if > half of leaves,
# as in error-filter.bf), or the whole site if they are >= 40% of sequences. Masking is applied to the ORIGINAL codons.
import json, sys
from Bio import Phylo
from io import StringIO
js, aln, out, rep = sys.argv[1:5]
THR, RATIO, SITE = 100., 20., 0.4
d = json.load(open(js))
rates = d['fits']['Unconstrained model']['Rate Distributions']['Test']
nr = len(rates); p_err = rates['0']['proportion']; p_fast = rates[str(nr-1)]['proportion']
prior = p_err/(1-p_err); prior_ratio = min(1e25, p_err/p_fast)
def readfa(p):
    dd={};n=None
    for l in open(p):
        l=l.rstrip('\n')
        if l.startswith('>'): n=l[1:].split()[0]; dd[n]=[]
        elif n: dd[n].append(l.strip())
    return {k:''.join(v) for k,v in dd.items()}
A = readfa(aln); taxa=list(A); nsite=len(A[taxa[0]])//3
tree = Phylo.read(StringIO(d['input']['trees']['0']+';'), 'newick')
leaves = set(t.name for t in tree.get_terminals())
up = {k.upper(): k for k in taxa}
ba = d['branch attributes']['0']
mask = {t:set() for t in taxa}; whole=set(); events=[]
for cl in tree.find_clades():
    if cl == tree.root or cl.name not in ba: continue
    post = ba[cl.name].get('Posterior prob omega class by site')
    if post is None: continue
    ld = set(t.name for t in cl.get_terminals())
    if not cl.is_terminal() and 2*len(ld) > len(leaves): ld = leaves - ld
    for s in range(nsite):
        pe = post[0][s]; pf = post[nr-1][s]
        if pe is None: continue
        bf = pe/(1-pe)/prior if pe < 1 else 1e25
        bf2 = pe/pf/prior_ratio if pf > 0 else 1e25
        if bf >= THR and bf2 >= RATIO:
            events.append((s+1, cl.name, len(ld)))
            if cl.is_terminal(): mask[up[cl.name.upper()]].add(s)
            elif len(ld)/len(leaves) >= SITE: whole.add(s)
            else:
                for t in ld: mask[up[t.upper()]].add(s)
nm=0
with open(out,'w') as o:
    for t in taxa:
        seq=A[t]; c=[seq[3*i:3*i+3] for i in range(nsite)]
        for s in mask[t]|whole:
            if c[s]!='---': nm+=1
            c[s]='---'
        o.write(f'>{t}\n'+''.join(c)+'\n')
r = dict(error_sink_proportion=p_err, error_sink_omega=rates['0']['omega'], fast_class=rates[str(nr-1)], n_branch_site_events=len(events),
         n_sites_with_event=len(set(e[0] for e in events)), n_whole_sites_masked=len(whole), n_codon_cells_masked=nm,
         per_taxon={t:len(mask[t]) for t in taxa if mask[t]}, events_terminal=sum(1 for e in events if e[2]==1))
json.dump(r, open(rep,'w'), indent=1); print({k:v for k,v in r.items() if k!='per_taxon'}); print(sorted(r['per_taxon'].items(), key=lambda x:-x[1])[:12])
