"""Credibility test for protist proteins whose best mycobacterial match is an FSM protein.
For each protist protein P with best FSM hit T: a = identity(P,T); region = P's span aligned to T.
b = best identity of P to any non-FSM mycobacterial protein covering >=80% of that span;
c = best identity of P to a protein from any OTHER FSM genome covering >=80% of the span.
Credible FSM evidence requires a>=95, a-b>=3, and c>=90 (i.e. the FSM gene is consistent across FSM genomes,
not a contaminant contig unique to one MAG)."""
import collections,csv
fsm=set(x.replace('GB_','').replace('RS_','') for x in open('fsm_ids.txt').read().split())
sp={}
for r in csv.DictReader(open('../../../Publication/scripts/inputs/tsv/S1_genome_metadata.tsv'),delimiter='\t'):
    sp[r['tip'].replace('GB_','').replace('RS_','').replace('.fa','')]=r['species']
H=collections.defaultdict(list)
for l in open('cand13_vs_myco_us.tsv'):
    q,s,p,L,qs,qe,ql,b=l.rstrip().split('\t'); H[q].append((s,float(p),int(L),int(qs),int(qe),float(b)))
out=open('protist_candidates_credibility.tsv','w')
out.write('protist_protein\tFSM_target\tFSM_species\ta_id_to_FSM_target\tspan_aa\tb_best_nonFSM_id\tb_species\tc_best_otherFSM_genome_id\tcredible\n')
nc=0
for q,hs in H.items():
    F=[h for h in hs if h[0].split('|')[0] in fsm]
    if not F: continue
    T=max(F,key=lambda h:h[5]); span=T[4]-T[3]+1
    def ov(h): return max(0,min(h[4],T[4])-max(h[3],T[3])+1)>=0.8*span
    N=[h for h in hs if h[0].split('|')[0] not in fsm and ov(h)]
    O=[h for h in F if h[0].split('|')[0]!=T[0].split('|')[0] and ov(h)]
    b=max(N,key=lambda h:h[1]) if N else None; c=max((h[1] for h in O),default=0)
    a=T[1]; bv=b[1] if b else 0
    cred=a>=95 and a-bv>=3 and c>=90; nc+=cred
    out.write(f"{q}\t{T[0]}\t{sp.get(T[0].split('|')[0])}\t{a}\t{span}\t{bv}\t{sp.get(b[0].split('|')[0]) if b else ''}\t{c}\t{cred}\n")
print('credible',nc)
