"""Summarise competitive mapping of Lake Pavin protist SAG reads (first 10M R1 per run) to 862 Mycobacteriaceae genomes.
A read counts if >=95% identity over >=80% of its length; 'unique' adds MAPQ>=20.
Breadth: number of distinct 10-kb windows hit per genome. Reads piled on a few windows are contaminant contigs in
reference assemblies (e.g. a Ralstonia rRNA operon in GCA_001665255.1, a Pseudomonas rRNA operon in GCA_021648845.1,
plant 18S in GCF_002838065.1); a true genome present in a SAG should be hit across >=20 windows."""
import glob,collections,csv,os
fsm=set(x.replace('GB_','').replace('RS_','') for x in open('../fsm_ids.txt').read().split())
clade={}
for l in open('../clade.tsv'):
    t,c=l.rstrip('\n').split('\t'); clade[t.replace('GB_','').replace('RS_','')]=c.split(' (')[0]
runs={r['run_accession']:r for r in csv.DictReader(open('runs.tsv'),delimiter='\t')}
out=open('sag_mapping_summary.tsv','w')
out.write('run\tsample\treads_screened\treads_ge95\treads_ge95_unique\tFSM_ge95\tFSM_ge95_unique\tFSM_ge99_unique\tnonFSM_ge95_unique\ttop_genomes_unique\tgenomes_ge20_windows\tmax_FSM_windows\n')
tot=collections.Counter()
for f in sorted(glob.glob('paf/*.paf')):
    run=os.path.basename(f)[:-4]
    if not os.path.exists(f+'.done'): continue
    c=collections.Counter(); g=collections.Counter(); W=collections.defaultdict(set)
    for l in open(f):
        p=l.split('\t'); ql=int(p[1]); m=int(p[9]); a=int(p[10]); mq=int(p[11])
        if m/a<0.95 or a<0.8*ql: continue
        gen=p[5].split('|')[0]; isf=gen in fsm
        c['ge95']+=1; c['FSM_ge95']+=isf
        if mq>=20:
            c['uniq']+=1; g[gen]+=1; W[gen].add((p[5],int(p[7])//10000))
            if isf: c['FSM_u']+=1; c['FSM_u99']+= m/a>=0.99
            else: c['non_u']+=1
    top=';'.join(f'{k}({clade.get(k,"?")}):{v}reads/{len(W[k])}win' for k,v in g.most_common(5))
    broad=[k for k in W if len(W[k])>=20]; fsmw=max([len(W[k]) for k in W if k in fsm],default=0)
    out.write(f"{run}\t{runs[run]['sample_alias']}\t10000000\t{c['ge95']}\t{c['uniq']}\t{c['FSM_ge95']}\t{c['FSM_u']}\t{c['FSM_u99']}\t{c['non_u']}\t{top}\t{';'.join(broad)}\t{fsmw}\n")
    tot.update(c)
print(dict(tot))
