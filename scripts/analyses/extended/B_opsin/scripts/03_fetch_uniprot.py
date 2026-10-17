"""Paginated UniProt REST download (no user identity sent). Usage: query outprefix"""
import sys, re, urllib.request, time
q, out = sys.argv[1], sys.argv[2]
for fmt, ext, extra in [('tsv','tsv','&fields=accession,id,protein_name,organism_name,organism_id,lineage,length,reviewed,xref_pfam'),('fasta','faa','')]:
    url=f"https://rest.uniprot.org/uniprotkb/search?query={q}&format={fmt}&size=500{extra}"
    first=True
    with open(f"{out}.{ext}","w") as fh:
        while url:
            for t in range(5):
                try:
                    r=urllib.request.urlopen(url,timeout=120); body=r.read().decode(); break
                except Exception as e: time.sleep(5)
            if fmt=='tsv' and not first: body=body.split('\n',1)[1]
            fh.write(body); first=False
            m=re.search(r'<([^>]+)>; rel="next"', r.headers.get('Link','') or ''); url=m.group(1) if m else None
            print(r.headers.get('X-UniProt-Release'), file=sys.stderr)
