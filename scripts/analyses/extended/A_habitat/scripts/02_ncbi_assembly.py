"""Fetch NCBI assembly summaries (BioProject, BioSample, 'derived from metagenome' flag) for all GTDB-reference genomes labelled MAG in S1. Cached; no email sent."""
import pandas as pd, requests, json, time, os
O='/path/to/project/analyses/extended/A_habitat/'
C=O+'cache/ncbi/assembly_batch.json'
d=pd.read_csv(O+'mag_master.tsv',sep='\t',dtype=str)
acc=sorted(d[d.origin=='GTDB reference'].t.str.replace(r'^(GB|RS)_','',regex=True))
res=json.load(open(C)) if os.path.exists(C) else {}
todo=[a for a in acc if a not in res]
E='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
for i in range(0,len(todo),80):
    b=todo[i:i+80]
    term=' OR '.join(f'{a}[Assembly Accession]' for a in b)
    r=requests.post(E+'esearch.fcgi',data=dict(db='assembly',term=term,retmax=500,retmode='json',tool='extended_habitat')).json()
    ids=r['esearchresult']['idlist']; time.sleep(0.5)
    s=requests.post(E+'esummary.fcgi',data=dict(db='assembly',id=','.join(ids),retmode='json',tool='extended_habitat')).json()['result']
    time.sleep(0.5)
    for u in s['uids']:
        x=s[u]
        for a in [x.get('assemblyaccession'),x.get('synonym',{}).get('genbank'),x.get('synonym',{}).get('refseq')]:
            if a in b: res[a]=dict(biosample=x.get('biosampleaccn'),bioproject=';'.join(g.get('bioprojectaccn') for g in x.get('gb_bioprojects',[])),
                     submitter=x.get('submitterorganization'),excl=';'.join(x.get('exclfromrefseq',[])),isolate=(x.get('biosource') or {}).get('isolate'),
                     strain=(x.get('biosource') or {}).get('infraspecieslist',[{}])[0].get('sub_value') if (x.get('biosource') or {}).get('infraspecieslist') else None)
    print(len(res),'/',len(acc))
json.dump(res,open(C,'w'))
missing=[a for a in acc if a not in res]; print('missing',missing)
