import requests, time, xml.etree.ElementTree as ET
E='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
def bs_attrs(accs):
    out={}
    accs=list(dict.fromkeys(a for a in accs if isinstance(a,str)))
    for i in range(0,len(accs),100):
        ch=accs[i:i+100]
        for attempt in range(3):
            try:
                r=requests.post(E+'esearch.fcgi',data={'db':'biosample','term':' OR '.join(a+'[accn]' for a in ch),'retmax':500,'retmode':'json'}).json()
                ids=r['esearchresult']['idlist']; time.sleep(0.4)
                x=requests.post(E+'efetch.fcgi',data={'db':'biosample','id':','.join(ids)}).text; time.sleep(0.4)
                root=ET.fromstring(x); break
            except Exception as e: time.sleep(3)
        for b in root.findall('BioSample'):
            a=b.get('accession'); d={'_title':(b.findtext('Description/Title') or '')}
            o=b.find('Description/Organism'); d['_organism']=o.get('taxonomy_name') if o is not None else ''
            for at in b.findall('Attributes/Attribute'):
                d[at.get('harmonized_name') or at.get('attribute_name')]=at.text
            out[a]=d
    return out
