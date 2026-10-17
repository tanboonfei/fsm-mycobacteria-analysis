"""Read-level classification from minimap2 -c PAF with secondary alignments (-N 30).
Per read: keep alignments covering >=80% of the read at >=95% identity (NM-based: (alnlen-NM)/alnlen via de or
matches/alnlen with -c); take those with the maximal AS. Read = FSM if all best hits are FSM genomes, nonFSM if none,
mixed otherwise. Breadth = distinct 10-kb windows on best-hit genomes (first best hit).
Usage: classify_paf.py <paf> [<paf> ...]  -> prints per file: reads, FSM, nonFSM, mixed, top genomes with windows."""
import sys,collections,re
fsm=set(x.replace('GB_','').replace('RS_','') for x in open('/path/to/project/analyses/protist/B/fsm_ids.txt').read().split())
def run(path):
    best={}
    for l in open(path):
        p=l.rstrip('\n').split('\t')
        ql=int(p[1]); m=int(p[9]); a=int(p[10])
        if a<0.8*ql or m/a<0.95: continue
        AS=int(re.search(r'AS:i:(-?\d+)',l).group(1))
        r=p[0]
        if r not in best or AS>best[r][0]: best[r]=[AS,[(p[5],int(p[7]),m/a)]]
        elif AS==best[r][0]: best[r][1].append((p[5],int(p[7]),m/a))
    c=collections.Counter(); G=collections.Counter(); W=collections.defaultdict(set); FW=collections.defaultdict(set)
    for r,(AS,hits) in best.items():
        gs={h[0].split('|')[0] for h in hits}
        isF=[g in fsm for g in gs]
        cls='FSM' if all(isF) else ('nonFSM' if not any(isF) else 'mixed')
        c[cls]+=1
        if cls=='FSM' and max(h[2] for h in hits)>=0.99: c['FSM_ge99']+=1
        h=hits[0]; g=h[0].split('|')[0]
        if cls!='mixed':
            G[g]+=1; W[g].add((h[0],h[1]//10000))
            if cls=='FSM': FW['FSMclade'].add((h[0],h[1]//10000))
    return c,G,W,FW
if __name__=='__main__':
    for f in sys.argv[1:]:
        c,G,W,FW=run(f)
        print(f,dict(c),'FSM windows',len(FW['FSMclade']),'| top:',';'.join(f'{g}:{n}r/{len(W[g])}w' for g,n in G.most_common(4)))
