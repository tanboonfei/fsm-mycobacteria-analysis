#!/usr/bin/env python3
"""Direct reduced-Ne test: within-population pi_N/pi_S for FSM vs RGM species (>=4 MAGs).
Higher pi_N/pi_S = less effective purifying selection = smaller Ne. Nei-Gojobori (NG86)."""
import csv, glob, os, subprocess, random, numpy as np
from collections import defaultdict, Counter
random.seed(42)

bases="TCAG"; codons=[a+b+c for a in bases for b in bases for c in bases]
aa="FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
code=dict(zip(codons,aa))
def transl(s):
    return "".join(code.get(s[i:i+3],"X") for i in range(0,len(s)-len(s)%3,3))

# NG86 synonymous/nonsynonymous sites per codon
def sites(cod):
    if cod not in code or code[cod]=="*": return (0.0,0.0)
    syn=0.0
    for i in range(3):
        for b in "ACGT":
            if b==cod[i]: continue
            mut=cod[:i]+b+cod[i+1:]
            if code.get(mut)==code[cod] and code[cod]!="*": syn+=1
    s=syn/3.0; return (s,3-s)
def diff(c1,c2):
    # pathway-averaged syn/nonsyn differences between two codons
    pos=[i for i in range(3) if c1[i]!=c2[i]]
    if not pos: return (0.0,0.0)
    if code.get(c1)=="*" or code.get(c2)=="*" or "X" in (code.get(c1,"X"),code.get(c2,"X")): return (0.0,0.0)
    import itertools
    sd=nd=0.0; paths=0
    for perm in itertools.permutations(pos):
        cur=c1; ok=True; psd=pnd=0
        for i in perm:
            nxt=cur[:i]+c2[i]+cur[i+1:]
            if code.get(nxt)=="*": ok=False; break
            if code.get(nxt)==code.get(cur): psd+=1
            else: pnd+=1
            cur=nxt
        if ok: sd+=psd; nd+=pnd; paths+=1
    if paths==0: return (0.0,0.0)
    return (sd/paths, nd/paths)

def pnps(seqs):
    # seqs: list of equal-length codon strings (gaps already removed columnwise)
    n=len(seqs); L=len(seqs[0])//3
    Sd=Nd=S=N=0.0; npair=0
    for a in range(n):
        for b in range(a+1,n):
            npair+=1
            for k in range(L):
                c1=seqs[a][3*k:3*k+3]; c2=seqs[b][3*k:3*k+3]
                if "-" in c1+c2 or "N" in c1+c2: continue
                s1=sites(c1); s2=sites(c2)
                S+=(s1[0]+s2[0])/2; N+=(s1[1]+s2[1])/2
                d=diff(c1,c2); Sd+=d[0]; Nd+=d[1]
    pS=Sd/S if S else 0; pN=Nd/N if N else 0
    return pN,pS

s1={r["genome_id_or_accession"]:r for r in csv.DictReader(open("../Supplementary_Table_S1.tsv"),delimiter="\t")}
def cl(c): c=str(c); return "FSM" if "FSM" in c else ("RGM" if "RGM" in c else None)
have=set(os.path.basename(d.rstrip("/")) for d in glob.glob("../bakta_out/*/"))
pops=defaultdict(list)
for g,r in s1.items():
    c=cl(r["clade"]); sp=r.get("species","").strip()
    if c in ("FSM","RGM") and g in have and sp and sp!="novel sp.": pops[(c,sp)].append(g)
pops={k:v for k,v in pops.items() if len(v)>=3}
# reproduce the S13a (>=4 MAGs) gene draws exactly: same RNG sequence over the >=4 populations in sorted order
import random as _r
_rng=_r.Random(42); GENESETS={}
_ge4=sorted(k for k,v in pops.items() if len(v)>=4)
_ge3=sorted(k for k,v in pops.items() if len(v)==3)

def load(g):
    # locus_tag -> gene name (single copy), and locus_tag -> nucleotide CDS
    name=Counter(); g2n={}
    for line in open(f"../bakta_out/{g}/{g}.tsv"):
        if line.startswith("#"): continue
        p=line.rstrip("\n").split("\t")
        if len(p)>=7 and p[1]=="cds" and p[6]:
            name[p[6]]+=1; g2n[p[5]]=p[6]
    sc={n for n,c in name.items() if c==1}
    seq={}; lt=None; buf=[]
    for line in open(f"../bakta_out/{g}/{g}.ffn"):
        if line.startswith(">"):
            if lt: seq[lt]="".join(buf)
            lt=line[1:].split()[0]; buf=[]
        else: buf.append(line.strip())
    if lt: seq[lt]="".join(buf)
    # gene name -> CDS (only single-copy named)
    out={}
    for ltag,gn in g2n.items():
        if gn in sc and ltag in seq: out[gn]=seq[ltag]
    return out

NG=80
results=[]
_extra=_r.Random(43)
for (c,sp) in _ge4+_ge3:
    members=pops[(c,sp)]
    gmaps={g:load(g) for g in members}
    shared=set.intersection(*[set(m) for m in gmaps.values()])
    genes=sorted(shared); (_rng if len(members)>=4 else _extra).shuffle(genes); genes=genes[:NG]
    cat={g:"" for g in members}
    used=0
    for gn in genes:
        seqs={g:gmaps[g][gn] for g in members}
        # write prot, mafft, back-translate
        prot={g:transl(seqs[g]) for g in seqs}
        tmp="/path/to/scratch/_pnps.faa"
        with open(tmp,"w") as o:
            for g in members: o.write(f">{g}\n{prot[g]}\n")
        try:
            al=subprocess.run(["mafft","--quiet","--auto",tmp],capture_output=True,text=True,timeout=60).stdout
        except Exception: continue
        ap={}; cur=None
        for line in al.splitlines():
            if line.startswith(">"): cur=line[1:].strip(); ap[cur]=""
            elif cur: ap[cur]+=line.strip()
        if len(ap)!=len(members): continue
        # back-translate
        bt={}
        for g in members:
            s=seqs[g]; p=ap[g]; codonlist=[s[i:i+3] for i in range(0,len(s)-len(s)%3,3)]
            r=[]; idx=0
            for c2 in p:
                if c2=="-": r.append("---")
                else: r.append(codonlist[idx] if idx<len(codonlist) else "---"); idx+=1
            bt[g]="".join(r)
        Lc=min(len(bt[g]) for g in members)//3
        # keep only columns with no gap across members
        keepcods=[]
        for k in range(Lc):
            cods=[bt[g][3*k:3*k+3] for g in members]
            if all("-" not in c2 and "N" not in c2 and len(c2)==3 for c2 in cods): keepcods.append(k)
        if not keepcods: continue
        for g in members: cat[g]+="".join(bt[g][3*k:3*k+3] for k in keepcods)
        used+=1
    if not cat[members[0]]: continue
    pN,pS=pnps([cat[g] for g in members])
    ratio=pN/pS if pS>0 else float('nan')
    results.append((c,sp,len(members),used,round(pS,4),round(pN,4),round(ratio,3)))
    print(f"{c} {sp[:28]:28s} n={len(members):2d} genes={used:3d} piS={pS:.4f} piN={pN:.4f} piN/piS={ratio:.3f}",flush=True)

print("\n=== SUMMARY: within-population piN/piS ===")
import statistics as st
for c in ("FSM","RGM"):
    rs=[r[6] for r in results if r[0]==c and r[6]==r[6]]
    print(f"{c}: {len(rs)} populations, median piN/piS = {st.median(rs):.3f}, values={sorted(round(x,3) for x in rs)}")
from scipy import stats
F=[r[6] for r in results if r[0]=="FSM" and r[6]==r[6]]; R=[r[6] for r in results if r[0]=="RGM" and r[6]==r[6]]
if F and R:
    u=stats.mannwhitneyu(F,R,alternative="greater")
    print(f"Mann-Whitney FSM>RGM piN/piS: p={u.pvalue:.3f}")
with open("../analyses/pnps_consistent.tsv","w") as o:
    o.write("clade\tspecies\tn_genomes\tn_genes\tpiS\tpiN\tpiN_piS\n")
    for r in results: o.write("\t".join(map(str,r))+"\n")
print("-> wrote pnps_popgen.tsv")
