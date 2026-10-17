# Build one extended codon alignment per control clade: primary 37-taxon alignment (unchanged rows)
# + the clade's genomes added by MAFFT --add to each of the 52 primary protein alignments.
import json, os, glob, subprocess, collections, sys
from Bio.Seq import Seq
B='/path/to/project/'
D=B+'analyses/extended/D_relax/'
W=D+'controls/build/'; os.makedirs(W,exist_ok=True); os.chdir(W)
cl=json.load(open(D+'controls/control_clades.json'))
groups=json.load(open(B+'analysis_3way/groups.json')); refpanel=set(groups['RGM']+groups['SGM'])
def readfa(p):
    d={};n=None
    for l in open(p):
        l=l.rstrip('\n')
        if l.startswith('>'): n=l[1:].split()[0]; d[n]=[]
        elif n: d[n].append(l.strip())
    return {k:''.join(v) for k,v in d.items()}
newg=[g for c in cl.values() for g in c['tips']]
nt={}
with open('ctrlprot.faa','w') as o:
    for g in newg:
        for k,s in readfa(B+'genes/'+g+'.ffn').items():
            s=s.upper()
            if len(s)%3: continue
            p=str(Seq(s).translate(table=11)).rstrip('*')
            if '*' in p: continue
            nid=g+'~'+k; nt[nid]=s; o.write('>'+nid+'\n'+p+'\n')
subprocess.run('diamond makedb --in ctrlprot.faa -d ctrlprot --quiet',shell=True,check=True)
ogs=sorted(os.path.basename(p)[:-5] for p in glob.glob(B+'analysis_3way/ogx/*.paln'))
with open('queries.faa','w') as o:
    for og in ogs:
        for k,s in readfa(B+f'analysis_3way/ogx/{og}.paln').items():
            if k.split('~')[0] in refpanel: o.write(f'>{og}|{k}\n'+s.replace('-','')+'\n')
subprocess.run('diamond blastp -q queries.faa -d ctrlprot -o hits.tsv --more-sensitive -k 200 -e 1e-10 -f 6 qseqid sseqid pident length qlen slen bitscore --threads 2 --quiet',shell=True,check=True)
best=collections.defaultdict(dict)
for l in open('hits.tsv'):
    q,s,pid,ln,ql,sl,bs=l.split('\t'); og=q.split('|')[0]; g=s.split('~')[0]
    if float(pid)<50 or int(ln)<0.7*int(ql): continue
    d=best[og].setdefault(g,{}); d[s]=max(d.get(s,0),float(bs))
chosen={}; stats=collections.Counter()
for og in ogs:
    chosen[og]={}
    for g in newg:
        d=best[og].get(g,{})
        if not d: stats['missing']+=1; continue
        r=sorted(d.items(),key=lambda x:-x[1])
        if len(r)>1 and r[1][1]>=0.9*r[0][1]: stats['ambiguous']+=1; continue
        chosen[og][g]=r[0][0]; stats['ok']+=1
print(stats)
json.dump(chosen,open('chosen.json','w'),indent=0)
prim=readfa(D+'inputs/three_way.codon.fa'); taxa37=list(prim)
occ={}
for cname,c in cl.items():
    concat=collections.defaultdict(list)
    for og in ogs:
        pa=readfa(B+f'analysis_3way/ogx/{og}.paln'); ntd=readfa(B+f'analysis_3way/ogx/{og}.nt')
        L0=len(next(iter(pa.values())))
        add={g:chosen[og][g] for g in c['tips'] if g in chosen[og]}
        if add:
            with open('tmp.add.faa','w') as o:
                for g,tid in add.items(): o.write(f'>NEW_{g}\n'+str(Seq(nt[tid]).translate(table=11)).rstrip('*')+'\n')
            with open('tmp.ex.faa','w') as o:
                for k,s in pa.items(): o.write(f'>{k}\n{s}\n')
            out=subprocess.run('mafft --add tmp.add.faa --quiet --thread 1 tmp.ex.faa',shell=True,capture_output=True,text=True,check=True).stdout
            al={};n=None
            for l in out.splitlines():
                if l.startswith('>'): n=l[1:].strip(); al[n]=''
                else: al[n]+=l.strip()
        else:
            al=dict(pa)
        L=len(next(iter(al.values())))
        # backtranslate every row on the full (possibly widened) alignment
        cod={}
        for k,a in al.items():
            if k.startswith('NEW_'): g=k[4:]; s=nt[add[g]]
            else: g=k.split('~')[0]; s=ntd[k].upper()
            i=0; o=[]
            for ch in a:
                if ch=='-': o.append('---')
                else: o.append(s[3*i:3*i+3]); i+=1
            cod[g]=o
        # drop columns gapped in every original sequence (insertions relative to the fixed primary alignment)
        orig=[g for g in cod if g in prim]
        keep=[j for j in range(L) if any(al_k[j]!='-' for k,al_k in al.items() if not k.startswith('NEW_'))]
        assert len(keep)==L0,(og,len(keep),L0)
        for g in taxa37+c['tips']:
            if g in cod: concat[g].append(''.join(cod[g][j] for j in keep))
            else: concat[g].append('---'*L0)
    seqs={g:''.join(v) for g,v in concat.items()}
    for g in taxa37: assert seqs[g]==prim[g], (cname,g)
    stops={'TAA','TAG','TGA'}
    with open(D+f'controls/{cname}.codon.fa','w') as o:
        for g in taxa37+c['tips']:
            s=seqs[g]
            if g not in prim: s=''.join('---' if s[i:i+3] in stops else s[i:i+3] for i in range(0,len(s),3))
            o.write(f'>{g}\n{s}\n')
    occ[cname]={g:sum(1 for og in ogs if g in chosen[og]) for g in c['tips']}
    print(cname, 'taxa', len(seqs), 'codons', len(seqs[taxa37[0]])//3, 'OG occupancy', occ[cname], flush=True)
json.dump(occ,open(D+'controls/occupancy.json','w'),indent=1)
