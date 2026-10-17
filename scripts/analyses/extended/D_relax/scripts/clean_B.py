# Codon-aware alignment cleaning of the primary 52-orthologue alignment (no row/column re-alignment):
# (1) HmmCleaner-like segment masking: per sequence, per-residue mean BLOSUM62 score against all other residues
#     in the column; 7-residue running mean; residues inside runs (>=5 residues) whose running mean < 0 are masked
#     (codon set to ---). Window/threshold chosen a priori, recorded here.
# (2) trimAl -automated1 on each masked protein alignment; kept protein columns map 1:1 to codon triplets.
import glob, os, subprocess, json, numpy as np
from Bio.Align import substitution_matrices
B='/path/to/project/'; D=B+'analyses/extended/D_relax/clean/'
W=D+'B_work/'; os.makedirs(W,exist_ok=True)
M=substitution_matrices.load('BLOSUM62')
WIN=7; THR=0.0; MINRUN=5
def readfa(p):
    d={};n=None
    for l in open(p):
        l=l.rstrip('\n')
        if l.startswith('>'): n=l[1:].split()[0]; d[n]=[]
        elif n: d[n].append(l.strip())
    return {k:''.join(v) for k,v in d.items()}
prim=readfa(B+'analyses/extended/D_relax/inputs/three_way.codon.fa'); taxa=list(prim)
ogs=sorted(os.path.basename(p)[:-5] for p in glob.glob(B+'analysis_3way/ogx/*.paln'))
out={t:[] for t in taxa}; off=0; stats={'masked_residues':0,'residues':0,'cols_in':0,'cols_kept':0}; per_taxon={t:0 for t in taxa}; per_og={}
for og in ogs:
    pa=readfa(B+f'analysis_3way/ogx/{og}.paln'); L=len(next(iter(pa.values())))
    g2k={k.split('~')[0]:k for k in pa}
    P={t:(pa[g2k[t]] if t in g2k else '-'*L) for t in taxa}
    C={t:[prim[t][3*(off+j):3*(off+j)+3] for j in range(L)] for t in taxa}
    # sanity: protein gaps match codon gaps
    for t in taxa:
        for j in range(L): assert (P[t][j]=='-')==(C[t][j]=='---'), (og,t,j)
    A=np.array([list(P[t]) for t in taxa])
    nmask=0
    for i,t in enumerate(taxa):
        idx=[j for j in range(L) if A[i,j]!='-']
        if not idx: continue
        sc=[]
        for j in idx:
            col=[A[k,j] for k in range(len(taxa)) if k!=i and A[k,j] not in '-X*']
            a=A[i,j]
            sc.append(np.mean([M[a][b] for b in col]) if col and a not in 'X*' else 0.0)
        sc=np.array(sc); n=len(sc)
        rm=np.array([sc[max(0,q-WIN//2):q+WIN//2+1].mean() for q in range(n)])
        low=rm<THR; q=0
        while q<n:
            if low[q]:
                r=q
                while r<n and low[r]: r+=1
                if r-q>=MINRUN:
                    for s in range(q,r): C[t][idx[s]]='---'; A[i,idx[s]]='-'
                    nmask+=r-q
                q=r
            else: q+=1
        stats['residues']+=n; per_taxon[t]+=0
    stats['masked_residues']+=nmask
    for i,t in enumerate(taxa): per_taxon[t]+=sum(1 for j in range(L) if P[t][j]!='-' and A[i,j]=='-')
    # trimAl automated1 on masked protein alignment
    with open(W+og+'.mask.faa','w') as o:
        for i,t in enumerate(taxa): o.write(f'>{t}\n'+''.join(A[i])+'\n')
    r=subprocess.run(['trimal','-in',W+og+'.mask.faa','-out',W+og+'.trim.faa','-automated1','-colnumbering'],capture_output=True,text=True,check=True).stdout
    keep=[int(x) for x in r.split('#ColumnsMap')[1].replace(',',' ').split()]
    stats['cols_in']+=L; stats['cols_kept']+=len(keep); per_og[og]={'cols':L,'kept':len(keep),'masked_residues':nmask}
    for t in taxa: out[t].append(''.join(C[t][j] for j in keep))
    off+=L
assert off*3==len(prim[taxa[0]])
with open(D+'three_way.cleanB.codon.fa','w') as o:
    for t in taxa: o.write(f'>{t}\n'+''.join(out[t])+'\n')
# mask-only alignment (no column removal) for a separate run is not produced; columns kept summary:
stats['per_taxon_masked']=per_taxon; stats['per_og']=per_og
json.dump(stats,open(D+'cleanB_stats.json','w'),indent=1)
print({k:v for k,v in stats.items() if not k.startswith('per')})
