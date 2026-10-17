"""Fetch BioProject titles/descriptions for every GTDB-reference genome's BioProject(s). No email sent."""
import json, os, time, requests, pandas as pd
G='/path/to/project/analyses/extended/G_culture/'
r=json.load(open(G+'cache/assembly_esummary.json'))
bps=sorted({p['bioprojectaccn'] for x in r.values() for p in x.get('gb_bioprojects',[])})
C=G+'cache/bioproject_esummary.json'; res=json.load(open(C)) if os.path.exists(C) else {}
E='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'; H={'User-Agent':'Mozilla/5.0 (research script)'}
def post(f,**kw):
    for t in range(5):
        try:
            q=requests.post(E+f,data=dict(tool='culture_check',**kw),headers=H,timeout=120); q.raise_for_status(); time.sleep(0.4); return q.json()
        except Exception as e: print('retry',e); time.sleep(3)
todo=[b for b in bps if b not in res]
for i in range(0,len(todo),100):
    b=todo[i:i+100]
    ids=post('esearch.fcgi',db='bioproject',term=' OR '.join(f'{x}[Project Accession]' for x in b),retmax=1000,retmode='json')['esearchresult']['idlist']
    sm=post('esummary.fcgi',db='bioproject',id=','.join(ids),retmode='json')['result']
    for u in sm['uids']:
        x=sm[u]; res[x['project_acc']]=dict(title=x.get('project_title',''),description=x.get('project_description',''),
            data_type=x.get('project_data_type',''),target_material=x.get('project_target_material',''),target_capture=x.get('project_target_capture',''),
            sample_scope=x.get('project_target_scope',''),method=x.get('project_methodtype',''))
    print(len(res),'/',len(bps))
json.dump(res,open(C,'w')); print('missing',[b for b in bps if b not in res])
