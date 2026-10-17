import sys, urllib.request, urllib.parse, json, time, re
for doi in sys.argv[1:]:
    u='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmode=json&term='+urllib.parse.quote(doi+'[doi]')
    ids=json.load(urllib.request.urlopen(u,timeout=30))['esearchresult']['idlist']; time.sleep(0.4)
    if not ids: print('##',doi,'NO PMID'); continue
    t=urllib.request.urlopen('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&rettype=abstract&retmode=text&id='+ids[0],timeout=30).read().decode(); time.sleep(0.4)
    open(re.sub(r'[^A-Za-z0-9]','_',doi)+'.txt','w').write(t)
    print('##',doi,'PMID',ids[0]); print(t[:3000])
