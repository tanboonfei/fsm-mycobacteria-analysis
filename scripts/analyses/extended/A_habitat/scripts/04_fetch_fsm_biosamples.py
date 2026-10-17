"""Fetch full NCBI BioSample records (title + all attributes) for FSM source samples (and MAG biosamples where no source sample). Cached; no email."""
import pandas as pd, requests, json, time, os, re, xml.etree.ElementTree as ET
O='/path/to/project/analyses/extended/A_habitat/'
C=O+'cache/ncbi/fsm_biosamples.json'
m=pd.read_csv(O+'../../../Publication/scripts/analyses/protist/A/mag_sample_metadata.tsv',sep='\t',dtype=str)
d=pd.read_csv(O+'mag_master.tsv',sep='\t',dtype=str)
asm=json.load(open(O+'cache/ncbi/assembly_batch.json'))
fs=d[d.clade3=='FSM']
ids=set(s for s in m[m.t.isin(fs.t)]['sample'].dropna() if re.match(r'SAM[NED]',s))
ids|={asm[a]['biosample'] for a in fs.acc if a in asm}
extra=['ERS4065084','ERS6344551','ERS4065080','ERS6344560','ERS11712895','ERS4649301','ERS4649307']
res=json.load(open(C)) if os.path.exists(C) else {}
todo=sorted(i for i in ids|set(extra) if i not in res)
E='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
for i in range(0,len(todo),100):
    b=todo[i:i+100]
    r=requests.post(E+'esearch.fcgi',data=dict(db='biosample',term=' OR '.join(x+'[accn]' for x in b),retmax=500,retmode='json',tool='extended_habitat')).json(); time.sleep(0.5)
    x=requests.post(E+'efetch.fcgi',data=dict(db='biosample',id=','.join(r['esearchresult']['idlist']),tool='extended_habitat')).text; time.sleep(0.5)
    root=ET.fromstring(x)
    for bs in root.findall('BioSample'):
        acc=bs.get('accession'); ids_={i.text for i in bs.findall('.//Id')}
        rec={'_title':(bs.findtext('.//Title') or ''),'_comment':' '.join((p.text or '') for p in bs.findall('.//Description/Comment/Paragraph'))}
        for a in bs.findall('.//Attribute'): rec[a.get('harmonized_name') or a.get('attribute_name')]=a.text
        rec['_ids']=sorted(i for i in ids_ if i)
        res[acc]=rec
        for e in extra:
            if e in ids_: res[e]=rec
json.dump(res,open(C,'w'),indent=0)
print(len(res),'missing',[i for i in todo if i not in res])
