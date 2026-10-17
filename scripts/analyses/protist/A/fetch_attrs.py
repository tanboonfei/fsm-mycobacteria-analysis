import pandas as pd, requests, time, json, xml.etree.ElementTree as ET
E='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
def bs_attrs(accs):
    out={}
    accs=list(dict.fromkeys(a for a in accs if isinstance(a,str)))
    for i in range(0,len(accs),100):
        ch=accs[i:i+100]
        r=requests.get(E+'esearch.fcgi',params={'db':'biosample','term':' OR '.join(a+'[accn]' for a in ch),'retmax':500,'retmode':'json'}).json()
        ids=r['esearchresult']['idlist']; time.sleep(0.4)
        x=requests.post(E+'efetch.fcgi',data={'db':'biosample','id':','.join(ids)}).text; time.sleep(0.4)
        root=ET.fromstring(x)
        for b in root.findall('BioSample'):
            a=b.get('accession'); d={'_title':(b.findtext('Description/Title') or ''),'_organism':b.find('Description/Organism').get('taxonomy_name') if b.find('Description/Organism') is not None else ''}
            for ids_ in b.findall('Ids/Id'):
                if ids_.get('db')=='SRA': d['_sra']=ids_.text
            for at in b.findall('Attributes/Attribute'):
                d[at.get('harmonized_name') or at.get('attribute_name')]=at.text
            l=b.find('Links')
            if l is not None:
                d['_bioproject']=';'.join(x.get('label') or x.text or '' for x in l.findall('Link') if x.get('target')=='bioproject')
            out[a]=d
    return out
m=pd.read_csv('mag_samples.tsv',sep='\t'); m=m[m.culture_status=='MAG / uncultured']
a2b=json.load(open('asm2biosample.json'))
m['mag_biosample']=m.genome_id_or_accession.map(a2b)
A=bs_attrs(m.mag_biosample.tolist())
# derived-from source sample for NCBI MAGs
import re
def src(d):
    for k,v in d.items():
        if v and re.search(r'deriv|source_mat|sample_derived',k,re.I):
            s=re.findall(r'SAM[NED][A-Z]?\d+',v)
            if s: return s
    return []
m['source_sample']=[ (';'.join(src(A.get(b,{}))) or None) if isinstance(b,str) else s for b,s in zip(m.mag_biosample,m.derived_from_sample)]
m.loc[m.derived_from_sample.notna(),'source_sample']=m.derived_from_sample
srcs=[x for s in m.source_sample.dropna() for x in s.split(';')]
B=bs_attrs(srcs)
json.dump({'mag_bs':A,'src_bs':B},open('biosample_attrs.json','w'))
m.to_csv('mag_samples2.tsv',sep='\t',index=False)
print(m.source_sample.notna().groupby(m.group).mean(), len(B))
