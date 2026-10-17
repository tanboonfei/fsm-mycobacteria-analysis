"""Classify SILVA/nt 16S hits to FSM and fetch INSDC source metadata.
Usage: python3 analyse_16S.py <hits.tsv> <db: silva|nt> <outprefix>
A hit is retained at >=98.7% identity over >=300 bp and >=90% of the shorter of query/subject.
FSM-specific: the hit sequence's best match in the local 16S panel (616 reference + 72 MAG 16S) is an FSM 16S,
and exceeds its best non-FSM match by >=0.5 points."""
import sys,subprocess,json,collections,urllib.request,urllib.parse,time,re,os
hits,db,pre=sys.argv[1:4]
grp=json.load(open('panel_groups.json'))
fsm=set(open('fsm_ids.txt').read().split())
def isfsm(sid):
    return grp.get(sid,'').startswith('FSM') or sid.split('_16S')[0] in fsm
qfsm={l[1:].split()[0] for l in open('q_fsm16S.fa') if l[0]=='>'}
keep=collections.defaultdict(list)
for l in open(hits):
    p=l.rstrip('\n').split('\t')
    q,s,pid,L,ql,sl=p[0],p[1],float(p[2]),int(p[3]),int(p[4]),int(p[5])
    if pid<98.7 or L<300 or L<0.9*min(ql,sl): continue
    keep[s].append((q,pid,L))
print(len(keep),'subjects retained')
# fetch subject sequences
ids=list(keep)
fa=pre+'_subjects.fa'
if os.path.exists(fa) and os.path.getsize(fa)>0: pass
elif db=='silva':
    open(pre+'_ids.txt','w').write('\n'.join(ids)+'\n')
    subprocess.run(['seqkit','grep','-f',pre+'_ids.txt','silva/silva144.fa','-o',fa,'-j','8'],check=True)
else:
    with open(fa,'w') as o:
        for i in range(0,len(ids),200):
            u='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id='+','.join(ids[i:i+200])
            o.write(urllib.request.urlopen(u,timeout=120).read().decode()); time.sleep(0.4)
# panel classification of each subject
out=subprocess.run(['blastn','-query',fa,'-db','panel16S','-max_target_seqs','2000','-outfmt','6 qseqid sseqid pident length qlen','-num_threads','8'],capture_output=True,text=True).stdout
bestF=collections.defaultdict(float); bestN=collections.defaultdict(lambda:(0.0,''))
for l in out.splitlines():
    q,s,pid,L,ql=l.split('\t'); pid=float(pid)
    if int(L)<min(300,0.9*int(ql)): continue
    q=q.split()[0]
    if isfsm(s): bestF[q]=max(bestF[q],pid)
    elif pid>bestN[q][0]: bestN[q]=(pid,s+':'+grp.get(s,'?'))
def key(s): return s.split()[0]
# metadata
def acc_of(s):
    a=s.split('.')[0] if db=='silva' else s.split('.')[0]
    return a
accs=sorted({acc_of(s) for s in ids})
mcache=pre+'_meta.json'
meta=json.load(open(mcache)) if os.path.exists(mcache) else {}
todo=[a for a in accs if a not in meta]
for i in range(0,len(todo),100):
    chunk=todo[i:i+100]
    data=urllib.parse.urlencode({'db':'nuccore','id':','.join(chunk),'rettype':'gb','retmode':'text'}).encode()
    for t in range(3):
        try: txt=urllib.request.urlopen('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi',data,timeout=180).read().decode(); break
        except Exception as e: time.sleep(5); txt=''
    for rec in txt.split('\n//'):
        m=re.search(r'^ACCESSION\s+(\S+)',rec,re.M)
        if not m: continue
        a=m.group(1)
        def q(f):
            mm=re.search(r'/'+f+r'="([^"]*)"',rec,re.S); return re.sub(r'\s+',' ',mm.group(1)) if mm else ''
        de=re.search(r'^DEFINITION\s+(.*?)\nACCESSION',rec,re.M|re.S)
        ti=re.search(r'^\s+TITLE\s+(.*?)\n\s+JOURNAL',rec,re.M|re.S)
        meta[a]=dict(definition=re.sub(r'\s+',' ',de.group(1)) if de else '',organism=q('organism'),isolation_source=q('isolation_source'),host=q('host'),country=q('country') or q('geo_loc_name'),title=re.sub(r'\s+',' ',ti.group(1)) if ti else '',env=('environmental_sample' in rec))
    time.sleep(0.4)
EUK=re.compile(r'amoeb|acanthamoeba|vermamoeba|hartmannella|dictyostel|naegleri|ciliat|tetrahymena|paramecium|euplotes|protist|protozo|flagellat|endosymb|symbiont|intracellular|co-cultur|coculture|alga|diatom|chlamydomonas|eukaryot|host cell|sponge|nematode|rotifer|daphnia|insect|fish|plant|root|rhizosph|leaf|animal|gut|feces|faec|human|patient|clinical|sputum',re.I)
json.dump(meta,open(mcache,'w'))
rows=[]
for s in ids:
    k=key(s); a=acc_of(s); md=meta.get(a,{})
    qbest=max(keep[s],key=lambda x:x[1])
    spec = bestF[k]>0 and bestF[k]>=bestN[k][0]+0.5
    ctx=' | '.join(md.get(f,'') for f in ('isolation_source','host','definition','title'))
    rows.append(dict(subject=s,accession=a,best_query=qbest[0],query_is_FSM=qbest[0] in qfsm,best_pid=qbest[1],aln=qbest[2],
        panel_best_FSM=bestF[k],panel_best_nonFSM=bestN[k][0],panel_best_nonFSM_hit=bestN[k][1],FSM_specific=spec,
        eukaryote_or_host_flag=bool(EUK.search(ctx)),flag_terms=';'.join(sorted({m.lower() for m in EUK.findall(ctx)})),**md))
cols=list(rows[0].keys()) if rows else []
with open(pre+'_hits_annotated.tsv','w') as o:
    o.write('\t'.join(cols)+'\n')
    for r in rows: o.write('\t'.join(str(r.get(c,'')) for c in cols)+'\n')
f=[r for r in rows if r['query_is_FSM'] and r['FSM_specific']]
print('FSM-query hits',sum(r['query_is_FSM'] for r in rows),'FSM-specific',len(f),'>=99%',sum(r['best_pid']>=99 for r in f),
      'flagged euk/host',sum(r['eukaryote_or_host_flag'] for r in f))
