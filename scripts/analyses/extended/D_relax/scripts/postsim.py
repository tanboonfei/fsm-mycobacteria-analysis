# Restore taxon names and impose the real alignment's missing-data pattern on a simulated codon alignment.
# usage: postsim.py <sim.fa> <real.fa> <out.fa>
import sys
def readfa(p):
    d={};n=None
    for l in open(p):
        l=l.rstrip('\n')
        if l.startswith('>'): n=l[1:].split()[0]; d[n]=[]
        elif n: d[n].append(l.strip())
    return {k:''.join(v) for k,v in d.items()}
sim=readfa(sys.argv[1]); real=readfa(sys.argv[2])
up={k.upper():k for k in real}
with open(sys.argv[3],'w') as o:
    for k,s in sim.items():
        name=up[k.upper()]; r=real[name]; assert len(r)==len(s),(name,len(r),len(s))
        out=[]
        for i in range(0,len(s),3):
            rc=r[i:i+3].upper()
            out.append(s[i:i+3] if all(c in 'ACGT' for c in rc) else '---')
        o.write(f'>{name}\n'+''.join(out)+'\n')
