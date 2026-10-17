import pandas as pd, requests, time, json, xml.etree.ElementTree as ET, re
E='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
m=pd.read_csv('mag_samples.tsv',sep='\t')
m=m[m.culture_status=='MAG / uncultured']
# assemblies lacking spire sample
acc=m[m.derived_from_sample.isna()].genome_id_or_accession.tolist()
asm2bs={}
for i in range(0,len(acc),100):
    chunk=acc[i:i+100]
    r=requests.get(E+'esearch.fcgi',params={'db':'assembly','term':' OR '.join(a+'[asac]' for a in chunk),'retmax':500,'retmode':'json'}).json()
    ids=r['esearchresult']['idlist']; time.sleep(0.4)
    s=requests.post(E+'esummary.fcgi',data={'db':'assembly','id':','.join(ids),'retmode':'json'}).json()['result']
    for k in s.get('uids',[]):
        d=s[k]; asm2bs[d['assemblyaccession']]=d.get('biosampleaccn')
    time.sleep(0.4)
print(len(acc),len(asm2bs))
json.dump(asm2bs,open('asm2biosample.json','w'))
