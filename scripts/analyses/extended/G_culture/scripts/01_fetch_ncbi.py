"""Fetch full NCBI assembly esummary and BioSample XML for every GTDB-reference genome in S1.
No email parameter sent; generic tool name. Cached in ../cache/."""
import pandas as pd, requests, json, time, os
G='/path/to/project/analyses/extended/G_culture/'
S1='/path/to/project/scripts/inputs/tsv/S1_genome_metadata.tsv'
s=pd.read_csv(S1,sep='\t')
acc=sorted(s[s.origin=='GTDB reference'].genome_id_or_accession)
E='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
H={'User-Agent':'Mozilla/5.0 (research script)'}
CA=G+'cache/assembly_esummary.json'
res=json.load(open(CA)) if os.path.exists(CA) else {}
todo=[a for a in acc if a not in res]
def post(f,**kw):
    for t in range(5):
        try:
            r=requests.post(E+f,data=dict(tool='culture_check',**kw),headers=H,timeout=120); r.raise_for_status(); time.sleep(0.4); return r
        except Exception as e: print('retry',e); time.sleep(3)
    raise SystemExit('fail')
for i in range(0,len(todo),100):
    b=todo[i:i+100]
    term=' OR '.join(f'{a}[Assembly Accession]' for a in b)
    ids=post('esearch.fcgi',db='assembly',term=term,retmax=1000,retmode='json').json()['esearchresult']['idlist']
    sm=post('esummary.fcgi',db='assembly',id=','.join(ids),retmode='json').json()['result']
    for u in sm['uids']:
        x=sm[u]
        for a in {x.get('assemblyaccession'),x.get('synonym',{}).get('genbank'),x.get('synonym',{}).get('refseq')}:
            if a in b:
                x.pop('meta',None) if False else None
                res[a]=x
    print(len(res),'/',len(acc),flush=True)
json.dump(res,open(CA,'w'))
miss=[a for a in acc if a not in res]; print('missing assemblies',miss)
# versioned mismatch: try without version
CB=G+'cache/biosample_xml.json'
bs=json.load(open(CB)) if os.path.exists(CB) else {}
bsl=sorted({res[a]['biosampleaccn'] for a in res if res[a].get('biosampleaccn')} - set(bs))
import xml.etree.ElementTree as ET
for i in range(0,len(bsl),100):
    b=bsl[i:i+100]
    ids=post('esearch.fcgi',db='biosample',term=' OR '.join(f'{x}[accn]' for x in b),retmax=1000,retmode='json').json()['esearchresult']['idlist']
    xml=post('efetch.fcgi',db='biosample',id=','.join(ids)).text
    root=ET.fromstring(xml)
    for e in root.findall('BioSample'):
        a=e.get('accession'); bs[a]=ET.tostring(e,encoding='unicode')
    print('biosamples',len(bs),flush=True)
json.dump(bs,open(CB,'w'))
print('missing biosamples',[x for x in bsl if x not in bs])
