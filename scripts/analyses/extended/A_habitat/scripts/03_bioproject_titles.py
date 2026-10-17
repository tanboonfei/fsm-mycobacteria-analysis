"""Fetch BioProject titles/submitters for PRJ* study ids (to collapse linked projects). Cached."""
import pandas as pd, requests, json, time, os
O='/path/to/project/analyses/extended/A_habitat/'
C=O+'cache/ncbi/bioproject.json'
d=pd.read_csv(O+'mag_master.tsv',sep='\t',dtype=str)
ids=sorted({s for s in d.study if s.startswith('PRJ')})
res=json.load(open(C)) if os.path.exists(C) else {}
todo=[i for i in ids if i not in res]
E='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
for i in range(0,len(todo),100):
    b=todo[i:i+100]
    r=requests.post(E+'esearch.fcgi',data=dict(db='bioproject',term=' OR '.join(f'{x}[PRJA]' for x in b),retmax=500,retmode='json',tool='extended_habitat')).json()
    time.sleep(0.5)
    s=requests.post(E+'esummary.fcgi',data=dict(db='bioproject',id=','.join(r['esearchresult']['idlist']),retmode='json',tool='extended_habitat')).json()['result']
    time.sleep(0.5)
    for u in s['uids']:
        x=s[u]; res[x['project_acc']]=dict(title=x.get('project_title'),submitter=x.get('submitter_organization') or ';'.join(x.get('organization_list',[]) if isinstance(x.get('organization_list'),list) else []),date=x.get('registration_date'))
json.dump(res,open(C,'w'),indent=0)
print(len(res),len(ids),[i for i in ids if i not in res])
