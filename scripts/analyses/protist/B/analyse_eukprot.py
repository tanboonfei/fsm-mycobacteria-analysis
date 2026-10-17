"""Summarise DIAMOND hits of EukProt v3 proteins to Mycobacteriaceae proteins (862 genomes).
Per protein: best hit; best FSM hit; best non-FSM hit. Protein assigned FSM if best FSM bitscore exceeds best non-FSM
bitscore and FSM pident >= 95. Dataset FSM-positive if >=2 such proteins."""
import collections,csv
fsm=set(x.replace('GB_','').replace('RS_','') for x in open('fsm_ids.txt').read().split())
clade={}
for l in open('clade.tsv'):
    t,c=l.rstrip('\n').split('\t'); clade[t.replace('GB_','').replace('RS_','')]=c.split(' (')[0]
species={}
for r in csv.DictReader(open('../../../Publication/scripts/inputs/tsv/S1_genome_metadata.tsv'),delimiter='\t'):
    species[r['tip'].replace('GB_','').replace('RS_','').replace('.fa','')]=r['species']
P=collections.defaultdict(lambda:{'F':(0,0,''),'N':(0,0,'')})
for l in open('eukprot_vs_myco.tsv'):
    q,s,pid,L,ql,sl,ev,bs=l.rstrip('\n').split('\t'); pid=float(pid); bs=float(bs)
    g=s.split('|')[0]; k='F' if g in fsm else 'N'
    if bs>P[q][k][1]: P[q][k]=(pid,bs,s)
meta={}
for r in csv.DictReader(open('eukprot/EukProt_included_data_sets.txt'),delimiter='\t'):
    meta[r['EukProt_ID']]=r
D=collections.defaultdict(lambda:collections.Counter()); best_g=collections.defaultdict(collections.Counter)
rows=[]
for q,v in P.items():
    ep=q.split('_')[0]
    F,N=v['F'],v['N']
    top=F if F[1]>N[1] else N
    g=top[2].split('|')[0]
    call='FSM' if (F[1]>N[1] and F[0]>=95) else ('nonFSM_myco' if top[0]>=95 else 'myco_90-95')
    D[ep][call]+=1; D[ep]['any90']+=1
    if top[0]>=95: best_g[ep][clade.get(g,'?')]+=1
    rows.append((q,ep,F[0],F[1],F[2],N[0],N[1],N[2],call,clade.get(g,'?'),species.get(g,'')))
with open('eukprot_protein_hits.tsv','w') as o:
    o.write('protein\tEukProt_ID\tbest_FSM_pid\tbest_FSM_bits\tbest_FSM_hit\tbest_nonFSM_pid\tbest_nonFSM_bits\tbest_nonFSM_hit\tcall\ttop_clade\ttop_species\n')
    for r in rows:o.write('\t'.join(map(str,r))+'\n')
with open('eukprot_dataset_summary.tsv','w') as o:
    o.write('EukProt_ID\tname\tsupergroup\ttaxogroup\tdata_type\tsource\tn_prot_myco_ge90\tn_ge95_nonFSM\tn_ge95_FSM_best\tclades_of_ge95_hits\tFSM_positive\n')
    for ep,c in sorted(D.items(),key=lambda x:-x[1]['any90']):
        m=meta.get(ep,{})
        o.write('\t'.join([ep,m.get('Name_to_Use',''),m.get('Supergroup_UniEuk',''),m.get('Taxogroup2_UniEuk',''),m.get('Data_Source_Type',''),m.get('Data_Source_Name','')[:80],
            str(c['any90']),str(c['nonFSM_myco']),str(c['FSM']),';'.join(f'{k}:{n}' for k,n in best_g[ep].most_common()),str(c['FSM']>=2)])+'\n')
print('proteins',len(P),'datasets',len(D),'datasets>=10 prot >=95%',sum(1 for c in D.values() if c['nonFSM_myco']+c['FSM']>=10),
      'FSM-best proteins',sum(c['FSM'] for c in D.values()),'FSM-positive datasets',sum(c['FSM']>=2 for c in D.values()))
