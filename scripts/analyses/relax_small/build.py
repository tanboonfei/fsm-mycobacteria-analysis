import json, os, glob, subprocess, collections
from Bio.Seq import Seq
B='/path/to/project/'
new=json.load(open('../small14.json')); origFSM=set(json.load(open(B+'analysis_3way/groups.json'))['FSM'])
def readfa(p):
    d={};n=None
    for l in open(p):
        l=l.rstrip('\n')
        if l.startswith('>'): n=l[1:].split()[0]; d[n]=[]
        elif n: d[n].append(l.strip())
    return {k:''.join(v) for k,v in d.items()}
# translate new genomes
nt={}; 
with open('newprot.faa','w') as o:
    for g in new:
        for k,s in readfa(B+'genes/'+g+'.ffn').items():
            s=s.upper(); 
            if len(s)%3: continue
            p=str(Seq(s).translate(table=11)).rstrip('*')
            if '*' in p: continue
            nid=g+'~'+k; nt[nid]=s; o.write('>'+nid+'\n'+p+'\n')
subprocess.run('diamond makedb --in newprot.faa -d newprot --quiet',shell=True,check=True)
ogs=sorted(set(os.path.basename(p).split('.')[0] for p in glob.glob(B+'analysis_3way/ogx/*.paln')))
with open('queries.faa','w') as o:
    for og in ogs:
        for k,s in readfa(B+f'analysis_3way/ogx/{og}.paln').items():
            if k.split('~')[0] in origFSM: o.write(f'>{og}|{k}\n'+s.replace('-','')+'\n')
subprocess.run('diamond blastp -q queries.faa -d newprot -o hits.tsv --more-sensitive -k 50 -e 1e-10 -f 6 qseqid sseqid pident length qlen slen bitscore --threads 4 --quiet',shell=True,check=True)
best=collections.defaultdict(dict)  # og -> genome -> {target: bitscore}
for l in open('hits.tsv'):
    q,s,pid,ln,ql,sl,bs=l.split('\t'); og=q.split('|')[0]; g=s.split('~')[0]
    if float(pid)<50 or int(ln)<0.7*int(ql): continue
    d=best[og].setdefault(g,{}); d[s]=max(d.get(s,0),float(bs))
chosen={}; stats=collections.Counter()
for og in ogs:
    chosen[og]={}
    for g in new:
        d=best[og].get(g,{})
        if not d: stats['missing']+=1; continue
        r=sorted(d.items(),key=lambda x:-x[1])
        if len(r)>1 and r[1][1]>=0.9*r[0][1]: stats['ambiguous']+=1; continue
        chosen[og][g]=r[0][0]; stats['ok']+=1
print(stats)
keep=[og for og in ogs if len(chosen[og])>=12]
print('OGs with >=12 of 14 new genomes:', len(keep), 'of', len(ogs))
json.dump({og:chosen[og] for og in keep}, open('chosen.json','w'), indent=0)
# codon alignments
others=None; concat=collections.defaultdict(str); taxa=None
for og in keep:
    prot=readfa(B+f'analysis_3way/ogx/{og}.paln'); ntd=readfa(B+f'analysis_3way/ogx/{og}.nt')
    seqs={}; nts={}
    for k,s in prot.items():
        t=k.split('~')[0]
        if t in origFSM: continue
        seqs[t]=s.replace('-',''); nts[t]=ntd[k].replace('-','').upper()
    for g,tid in chosen[og].items():
        seqs[g]=str(Seq(nt[tid]).translate(table=11)).rstrip('*'); nts[g]=nt[tid]
    with open(f'{og}.in.faa','w') as o:
        for t,s in seqs.items(): o.write(f'>{t}\n{s}\n')
    aln=subprocess.run(f'mafft --auto --quiet --thread 2 {og}.in.faa',shell=True,capture_output=True,text=True,check=True).stdout
    pa={}; n=None
    for l in aln.splitlines():
        if l.startswith('>'): n=l[1:].strip(); pa[n]=''
        else: pa[n]+=l.strip()
    L=len(next(iter(pa.values())))
    cod={}
    for t,a in pa.items():
        s=nts[t]; i=0; out=[]
        for ch in a:
            if ch=='-': out.append('---')
            else: out.append(s[3*i:3*i+3]); i+=1
        cod[t]=''.join(out)
    if taxa is None: taxa=sorted(set(cod)|set(new))
    for t in set(taxa)|set(cod):
        concat[t]+=cod.get(t,'---'*L)
    lens=set(len(v) for v in concat.values()); assert len(lens)==1, (og,lens)
# remove columns that are stop codons anywhere -> replace stop codons with gaps
stops={'TAA','TAG','TGA'}
with open('small_panel.codon.fa','w') as o:
    for t,s in concat.items():
        s=''.join('---' if s[i:i+3] in stops else s[i:i+3] for i in range(0,len(s),3))
        o.write(f'>{t}\n{s}\n')
print('taxa', len(concat), 'codons', len(next(iter(concat.values())))//3)
