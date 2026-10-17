"""Task 4: fetch Sandpiper (sandpiper.qut.edu.au, GTDB R232) per-species run tables for FSM and RGM r232 species,
# NOTE: the denominator actually used by script 10 is the full d__Bacteria table, fetched once with:
#   curl -s -A "academic-research-script/1.0" "https://sandpiper.qut.edu.au/api/taxonomy_search_csv/d__Bacteria?taxonomy_type=gtdb" | gzip > cache/sp_search/d__Bacteria_full.csv.gz
# The minimal table fetched below was used only for exploration and was deleted afterwards.
plus the d__Bacteria minimal table (all runs with any bacterial signal) as the denominator.
Sequential requests, >=1 s apart, cached to disk; no personal information in requests."""
import pandas as pd, requests, time, os, sys, urllib.parse
O='/path/to/project/analyses/extended/A_habitat/'
C=O+'cache/sp_search/'; os.makedirs(C,exist_ok=True)
B='https://sandpiper.qut.edu.au/api/'
H={'User-Agent':'academic-research-script/1.0 (low request rate)'}
s=pd.read_csv(O+'../../../Publication/scripts/analyses/protist/A/species_group_r232.tsv',sep='\t')
s=s[s.group.isin(['FSM','RGM','sisterRGM'])]
def get(url,fn,tries=4):
    if os.path.exists(fn) and os.path.getsize(fn)>0: return 'cached'
    for i in range(tries):
        try:
            r=requests.get(url,headers=H,timeout=600)
            if r.status_code==200:
                open(fn,'wb').write(r.content); return f'ok {len(r.content)}'
            if r.status_code in (404,400): open(fn,'w').write(f'#HTTP {r.status_code}\n'); return f'http {r.status_code}'
        except Exception as e: err=str(e)
        time.sleep(10*(i+1))
    return 'failed'
jobs=[('d__Bacteria_minimal',B+'taxonomy_search_csv_minimal/d__Bacteria?taxonomy_type=gtdb')]
for sp in s.species: jobs.append((sp.replace(' ','_'),B+'taxonomy_search_csv/'+urllib.parse.quote('s__'+sp)+'?taxonomy_type=gtdb'))
for name,url in jobs:
    st=get(url,C+name+'.csv'); print(name,st,flush=True)
    if st!='cached': time.sleep(1.5)
print('done')
