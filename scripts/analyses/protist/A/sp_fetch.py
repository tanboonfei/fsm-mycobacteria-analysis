import requests, json, sys, os, concurrent.futures as cf
B='https://sandpiper.qut.edu.au/api'
def flat(node,path,out):
    nm=node.get('name'); p=path+[nm] if nm else path
    if 'size' in node and node.get('size') is not None: out.append((';'.join(p),node['size']))
    for c in node.get('children',[]) or []: flat(c,p,out)
def get(run):
    fn=f'sp_profiles/{run}.json'
    if os.path.exists(fn): return run,'cached'
    try:
        r=requests.get(f'{B}/condensed/{run}?taxonomy_type=gtdb',timeout=120)
        if r.status_code!=200: json.dump({'error':r.status_code},open(fn,'w')); return run,r.status_code
        j=r.json(); out=[]; flat(j['condensed'],[],out)
        json.dump({'profile':out},open(fn,'w')); return run,len(out)
    except Exception as e: return run,str(e)
runs=[l.strip() for l in open(sys.argv[1]) if l.strip()]
with cf.ThreadPoolExecutor(8) as ex:
    for r,s in ex.map(get,runs): pass
print('done',len(runs))
