import pandas as pd, re, json
t=pd.read_csv('mag_sample_metadata.tsv',sep='\t')
d=json.load(open('biosample_attrs.json'))
def cat(r):
    s=' '.join(str(x) for x in [r.organism,r.iso,r.env,r.medium,r.local,r.title,r.gold] if isinstance(x,str)).lower()
    if re.search(r'sludge|wastewater|digester|bioreactor|sewage',s): return 'wastewater/engineered'
    if re.search(r'aquarium|tank|aquaculture',s): return 'aquarium/aquaculture'
    if re.search(r'sponge|coral|pancrustacean|gut|feces|faeces|human|sputum|root|rhizo|plant|leaf|host|daphnia|fish|insect',s): return 'host-associated'
    if re.search(r'sediment',s): return 'sediment'
    if re.search(r'soil|wetland|peat|permafrost',s): return 'soil/wetland'
    if re.search(r'lake|reservoir|pond|freshwater|river|estuar|bay|brackish|aquatic|water|sea|marine|ocean',s): return 'water column'
    return 'other/unknown'
t['category']=t.apply(cat,axis=1)
def flags(s):
    a=d['src_bs'].get(s) or d['mag_bs'].get(s) or {}
    txt=' '.join(f'{k}={v}' for k,v in a.items() if v).lower()
    return pd.Series({'sizefrac_info':bool(re.search(r'size.fraction|prefilt|µm|micron|filter',txt)),
                      'protist_or_euk_host':bool(re.search(r'protist|amoeb|ciliat|flagellat|acanthamoeba|paramecium|tetrahymena|dictyostelium|algal culture|co-?culture',txt)),
                      'host_field':(a.get('host') or '').strip()})
t=pd.concat([t,t['sample'].apply(flags)],axis=1)
# per MAG: majority category
per=t.groupby(['t','group']).agg(category=('category',lambda x:x.mode().iloc[0]),n_source_samples=('sample','nunique'),
     sizefrac_info=('sizefrac_info','any'),protist_or_euk_host=('protist_or_euk_host','any'),
     host=('host_field',lambda x:';'.join(sorted({h for h in x if h and h.lower() not in('not applicable','missing','none','na')}))),
     example_isolation_source=('iso',lambda x:next((v for v in x if isinstance(v,str)),'')),bioproject=('bioproject',lambda x:';'.join(sorted({str(v) for v in x if isinstance(v,str)}))[:200])).reset_index()
per.to_csv('mag_source_metadata_summary.tsv',sep='\t',index=False)
grp=per.group.replace({'sisterRGM':'RGM'})
tab=pd.crosstab(per.category,grp); print(tab); print((tab/tab.sum()*100).round(1))
print(per.groupby(grp)[['sizefrac_info','protist_or_euk_host']].sum())
print(per[(per.group=='FSM')&(per.category=='host-associated')][['t','host','example_isolation_source']])
print(per[per.protist_or_euk_host][['t','group','example_isolation_source']])
tab.to_csv('mag_source_category_by_clade.tsv',sep='\t')
