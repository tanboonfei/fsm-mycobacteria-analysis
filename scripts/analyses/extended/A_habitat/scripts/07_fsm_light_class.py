"""Task 3: classify each FSM genome's source environment by light/oxygen regime from sample metadata (rules documented in RESULTS.md)."""
import pandas as pd, json, re
O='/path/to/project/analyses/extended/A_habitat/'
m=pd.read_csv(O+'../../../Publication/scripts/analyses/protist/A/mag_sample_metadata.tsv',sep='\t',dtype=str)
d=pd.read_csv(O+'mag_master.tsv',sep='\t',dtype=str)
bs=json.load(open(O+'cache/ncbi/fsm_biosamples.json'))
asm=json.load(open(O+'cache/ncbi/assembly_batch.json'))
f=d[d.clade3=='FSM'].copy()
# explicit source-sample corrections (documented)
src_override={'GB_GCA_036782595.1':['SAMN04421877'],   # MAG re-binned from SRR3139094 = SAMN04421877 (Western Scheldt water, PRJNA309119)
              'GB_GCA_018401645.1':['SAMN19005952'],   # no source sample in table; MAG BioSample carries depth/layer
              'GB_GCA_026393245.1':['SAMN30025985']}
mag_bs_extra={'GB_GCA_943331415.1':'ERS11712895','GB_GCA_945882475.1':'ERS4065084','GB_GCA_945897705.1':'ERS6344551',
              'GB_GCA_964202555.1':'ERS4065080','GB_GCA_964203005.1':'ERS6344560','GB_GCA_903862565.1':'ERS4649301',
              'GB_GCA_903884065.1':'ERS4649301','GB_GCA_903877095.1':'ERS4649307'}  # 'sample derived from' of the MAG BioSample
def num(x):
    if x is None: return None
    x=str(x).replace(',','.'); mm=re.search(r'(-?\d+(?:\.\d+)?)',x)
    return float(mm.group(1)) if mm else None
rows=[]
for _,r in f.iterrows():
    samples=src_override.get(r.t) or [s for s in m[m.t==r.t]['sample'].dropna().unique()]
    if r.t in mag_bs_extra: samples=[mag_bs_extra[r.t]]
    recs=[bs[s] for s in samples if s in bs]
    a=re.sub(r'^(GB|RS)_','',r.t)
    if a in asm and asm[a]['biosample'] in bs and not recs: recs=[bs[asm[a]['biosample']]]
    txt=' || '.join(' ; '.join(f'{k}={v}' for k,v in x.items() if k!='_ids' and v) for x in recs)
    if not txt: txt=' '.join(str(r[c]) for c in ['title','iso','env','medium','local','depth','geo'] if isinstance(r[c],str))
    low=txt.lower().replace('to get rid of sediment','')  # filtration note, not a sediment sample
    depths=[]; oxy=[]
    for x in recs:
        for k in ['depth','geographic location (depth)']:
            if x.get(k) and num(x[k]) is not None and 'soil' not in str(x[k]).lower(): depths.append(num(x[k]))
        for mm in re.finditer(r'(\d+(?:\.\d+)?) ?m(?:eter)?(?: depth|\b)',x.get('_title','')+' '+x.get('_comment','')+' '+x.get('isolation_source','')):
            depths.append(float(mm.group(1)))
        mm=re.search(r'_(\d+)m$',x.get('_title',''))
        if mm: depths.append(float(mm.group(1)))
        if x.get('oxygen'): oxy.append(('mg/L',num(x['oxygen'])))
        if x.get('diss_oxygen'):
            u='mg/L' if 'mg' in x['diss_oxygen'] else 'umol'; oxy.append((u,num(x['diss_oxygen'])))
    mm=re.search(r'F1-(\d+)',txt)            # Lake Fuxian station F1, depth-coded sample names (F1-140, F1-20)
    if mm: depths.append(float(mm.group(1)))
    depth=max(depths) if depths else None
    anoxic_o2=any((u=='mg/L' and v is not None and v<1) or (u=='umol' and v is not None and v<31) for u,v in oxy)  # <1 mg/L ~ <31 umol/L
    rule=''
    if re.search(r'sludge|wastewater|digester|sewage|bioreactor',low): cat,rule='wastewater/engineered','keyword'
    elif re.search(r'sediment|pore water|mbsf',low) and r.t!='GB_GCA_036782595.1': cat,rule='sediment','keyword (sediment / pore water / mbsf)'
    elif re.search(r'soil|wetland',low): cat,rule='other: soil/wetland','keyword'
    elif re.search(r'sponge|haliclona',low): cat,rule='other: host-associated','keyword'
    elif 'pitcherplants' in r.study: cat,rule='other: host-associated','SPIRE study (pitcher-plant fluid)'
    elif 'tapwater' in r.study: cat,rule='other: drinking/tap water','SPIRE study name (tap water)'
    elif re.search(r'coassembly',low): cat,rule='water: depth unresolved','co-assembly across a depth profile'
    elif re.search(r'anoxic zone|hypolimnion|monimolimnion|anoxic water',low): cat,rule='dark/anoxic water','label (anoxic zone/hypolimnion/monimolimnion)'
    elif anoxic_o2: cat,rule='dark/anoxic water','dissolved O2 < 1 mg/L (31 umol/L)'
    elif depth is not None and depth>20: cat,rule='dark/anoxic water','depth > 20 m (below typical euphotic depth)'
    elif depth is not None and depth>10: cat,rule='transition (10-20 m)','10 m < depth <= 20 m'
    elif depth is not None or re.search(r'epilimnion|euphotic|surface layer|bloom|mixolimnion',low): cat,rule='photic water column','depth <= 10 m or epilimnion/surface/bloom label'
    elif re.search(r'lake|river|estuar|bay|sea water|freshwater|reservoir|water',low): cat,rule='photic water column (depth not reported)','water sample without depth: assumed near-surface'
    else: cat,rule='other: unknown','no usable metadata'
    rows.append(dict(t=r.t,study=r.study,species_key=r.species_key,origin=r.origin,rhodopsin=r.rhodopsin,salinity=r.salinity_habitat,
                     light_class=cat,rule=rule,depth_m=depth,o2=';'.join(f'{v}{u}' for u,v in oxy),samples_used=';'.join(samples[:5]),evidence=txt[:400]))
x=pd.DataFrame(rows)
x['light_class_main']=x.light_class.replace({'photic water column (depth not reported)':'photic water column'})
x.to_csv(O+'T3a_fsm_light_classification.tsv',sep='\t',index=False)
print(pd.crosstab(x.light_class,x.rhodopsin,margins=True))
print(x.groupby(['light_class','study']).agg(n=('t','size'),opsin=('rhodopsin',lambda s:(s=='present').sum()),depth=('depth_m',lambda s:','.join(sorted({str(v) for v in s.dropna()})))).to_string())
