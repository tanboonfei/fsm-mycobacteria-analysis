"""Submit one query per request to the NCBI BLAST URL API and poll (>=60 s between polls per RID).
Usage: urlblast.py <fasta> <program> <database> <entrez_query or ''> <out.tsv> [megablast on/off]"""
import sys,time,re,urllib.request,urllib.parse,json,os
def get(u,data=None,timeout=300):
    for i in range(8):
        try: return urllib.request.urlopen(u,data,timeout=timeout).read().decode()
        except Exception as e: print('retry',i,e,flush=True); time.sleep(30)
    raise SystemExit('network failure')
fa,prog,db,ent,out=sys.argv[1:6]; mb=sys.argv[6] if len(sys.argv)>6 else 'on'
recs=['>'+r for r in open(fa).read().split('>')[1:]]
def put(q):
    p={'CMD':'Put','PROGRAM':prog,'DATABASE':db,'QUERY':q,'HITLIST_SIZE':'100','EXPECT':'1e-10','TOOL':'protist_screen'}
    if prog=='blastn' and mb=='on': p['MEGABLAST']='on'
    if ent: p['ENTREZ_QUERY']=ent
    r=get('https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi',urllib.parse.urlencode(p).encode())
    m=re.search(r'RID = (\S+)',r)
    if not m:
        print('no RID, retrying',flush=True); time.sleep(60); return put(q)
    return m.group(1)
ridf=out+'.rids.json'
done=set(l.split('\t')[0] for l in open(out)) if os.path.exists(out) else set()
pending=json.load(open(ridf)) if os.path.exists(ridf) else {}
for r in recs:
    q=r.split()[0][1:]
    if q in pending or q in done: continue
    pending[q]=put(r); json.dump(pending,open(ridf,'w')); time.sleep(11)
print('submitted',len(pending),flush=True)
o=open(out,'a')
while pending:
    time.sleep(60)
    for qid,rid in list(pending.items()):
        s=get(f'https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID={rid}')
        st=re.search(r'Status=(\w+)',s); st=st.group(1) if st else '?'
        if st=='READY':
            t=get(f'https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_TYPE=Text&ALIGNMENT_VIEW=Tabular&ALIGNMENTS=500&DESCRIPTIONS=500&RID={rid}')
            rows=[l for l in t.splitlines() if l and not l.startswith('#') and l.count('\t')>=9]
            for l in rows: o.write(qid+'\t'+l.split('\t',1)[1]+'\n')
            o.write(qid+'\t#READY\t'+str(len(rows))+'\n'); o.flush(); print(qid,'READY',len(rows),flush=True); del pending[qid]; json.dump(pending,open(ridf,'w'))
        elif st in('FAILED','UNKNOWN'):
            print(qid,st,flush=True); del pending[qid]
        time.sleep(3)
print('done',flush=True)
