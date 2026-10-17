"""Build per-MAG master table (FSM/RGM/SGM MAGs) with source study assignment."""
import pandas as pd, json, re, os
R='/path/to/project/'
A=R+'Publication/scripts/analyses/protist/A/'
O=R+'analyses/extended/A_habitat/'
s1=pd.read_excel(R+'Publication/Additional_file_2.xlsx','S1')
s1=s1[s1.clade.isin(['FSM (small-genome)','RGM (large-genome)','SGM'])].copy()
s1['t']=s1.tip.str.replace(r'\.fa$','',regex=True)
s1['clade3']=s1.clade.str.split(' ').str[0]
m=pd.read_csv(A+'mag_sample_metadata.tsv',sep='\t',dtype=str)
m2=pd.read_csv(A+'mag_samples2.tsv',sep='\t',dtype=str).set_index('t')
mo=pd.read_csv(R+'microntology.tsv.gz',sep='\t',dtype=str)
samp2study=dict(zip(mo.ena_ers_sample_id,mo.study_name)); samp2study.update({k:v for k,v in zip(mo.sample_name,mo.study_name) if isinstance(k,str) and k!='NA'})
asm=json.load(open(O+'cache/ncbi/assembly_batch.json')) if os.path.exists(O+'cache/ncbi/assembly_batch.json') else {}
bp2study=mo.groupby('ena_project_id').study_name.agg(lambda s:s.mode().iloc[0]).to_dict()
bs=json.load(open(R+'analyses/protist/A/biosample_attrs.json'))
sp=pd.read_csv(R+'mycobacterium_spire_genomes.tsv',sep='\t',dtype=str)
spire_src=dict(zip(sp.genome_id,sp.derived_from_sample))
# manual collapse of linked BioProjects not in SPIRE (same title stem / submitter); documented in RESULTS.md
manual={'PRJNA648801':'Wang_2018_wastewater',   # same HKU lab, same Shatin activated-sludge plant
        'PRJNA375617':'Maresca_2018_Chesapeake'}  # paired Delaware Bay project
genome_override={'GB_GCA_036782595.1':('PRJNA309119_Eurytemora_copepod_water','manual: reanalysis of SRX1557323 (SRP069133 = PRJNA309119)')}
BPT=json.load(open(O+'cache/ncbi/bioproject.json')) if os.path.exists(O+'cache/ncbi/bioproject.json') else {}
def study_of(t,rows):
    if t in genome_override: return genome_override[t]
    samples=[s for s in rows['sample'].dropna()]
    if t in spire_src and isinstance(spire_src[t],str): samples=[spire_src[t]]+samples
    st=[samp2study[s] for s in samples if s in samp2study]
    if st: return pd.Series(st).mode().iloc[0],'SPIRE study (sample)'
    bps=[b for b in rows['bioproject'].dropna() for b in b.split(';')]
    st=[bp2study[b] for b in bps if b in bp2study]
    if st: return pd.Series(st).mode().iloc[0],'SPIRE study (BioProject)'
    if bps:
        b=pd.Series(bps).mode().iloc[0]; return manual.get(b,b),'source BioProject'
    a=re.sub(r'^(GB|RS)_','',t)
    if a in asm and asm[a]['bioproject']:
        b=asm[a]['bioproject'].split(';')[0]; st=bp2study.get(b)
        return (st or manual.get(b,b)),'NCBI assembly BioProject'
    mb=m2.mag_biosample.get(t) if t in m2.index else None
    if isinstance(mb,str) and mb in bs['mag_bs'] and bs['mag_bs'][mb].get('_bioproject'):
        b=bs['mag_bs'][mb]['_bioproject']
        st=bp2study.get(b)
        return (st or manual.get(b,b)),'MAG BioSample BioProject'
    return 'unknown_'+t,'unresolved (own study)'
out=[]
for t,rows in m.groupby('t'):
    st,how=study_of(t,rows)
    txt=lambda c:' | '.join(sorted({str(x) for x in rows[c].dropna() if str(x) not in('None','nan','not applicable','missing')}))[:300]
    out.append(dict(t=t,study=st,study_source=how,n_source_samples=rows['sample'].nunique(),
        bioprojects=';'.join(sorted({b for b in rows['bioproject'].dropna()})),
        title=txt('title'),organism=txt('organism'),iso=txt('iso'),env=txt('env'),medium=txt('medium'),local=txt('local'),
        gold=txt('gold'),depth=txt('depth'),geo=txt('geo'),host=txt('host'),size=txt('size'),
        samples=';'.join(sorted(rows['sample'].dropna().unique()))[:2000]))
d=pd.DataFrame(out)
# collapse linked BioProjects: identical normalised title + submitter
def norm(t): return re.sub(r'[0-9]+','#',re.split(r' - |: ',t)[0].lower()).strip()
key={}
for st in d.study.unique():
    if st.startswith('PRJ') and st in BPT and BPT[st].get('title'):
        k=(norm(BPT[st]['title']),str(BPT[st].get('submitter')))
        key.setdefault(k,[]).append(st)
coll={}
for k,v in key.items():
    if len(v)>1 and not k[0].startswith('mycobacterium') and 'genome sequencing' not in k[0]:
        for x in v: coll[x]=sorted(v)[0]+'_linked'
d['study_raw']=d.study; d['study']=d.study.map(lambda x:coll.get(x,x))
print('collapsed by title/submitter:',{k:v for k,v in key.items() if len(v)>1})
d=d.merge(s1[['t','origin','clade3','species','biome','salinity_habitat','rhodopsin']],on='t',how='inner')
gs=pd.read_csv(A+'genome_species_r232.tsv',sep='\t').set_index('t').sp
d['species_r232']=d.t.map(gs)
d['acc']=d.t.str.replace(r'^(GB|RS)_','',regex=True)
d['ncbi_metagenome_derived']=d.acc.map(lambda a:('derived from metagenome' in asm[a]['excl']) if a in asm else None)
d['ncbi_bioproject']=d.acc.map(lambda a:asm[a]['bioproject'] if a in asm else None)
d['ncbi_submitter']=d.acc.map(lambda a:asm[a]['submitter'] if a in asm else None)
d['source_set']=d.origin.map({'SPIRE MAG':'SPIRE MAG'}).fillna(d.ncbi_metagenome_derived.map({True:'GTDB MAG (NCBI metagenome-derived)',False:'GTDB genome not flagged metagenome-derived'}))
# species unit: GTDB species name from S1; for 'novel sp.' FSM genomes the skani 95% ANI clusters
# (analyses/unassigned_ani/clusters.tsv); for other novel SPIRE MAGs the SPIRE 95%-ANI species cluster.
ani=pd.read_csv(R+'Publication/scripts/analyses/unassigned_ani/clusters.tsv',sep='\t'); ani['t']=ani.tip.str.replace(r'\.fa$','',regex=True)
ani=dict(zip(ani.t,'FSM_'+ani.cluster)); spc=dict(zip(sp.genome_id,sp.spire_cluster))
def spkey(r):
    if isinstance(r.species,str) and r.species!='novel sp.': return r.species
    if r.t in ani: return ani[r.t]
    if r.t in spc and isinstance(spc[r.t],str): return 'SPIRE_'+spc[r.t]
    return 'single_'+r.t
d['species_key']=d.apply(spkey,axis=1)
d.to_csv(O+'mag_master.tsv',sep='\t',index=False)
print(d.groupby('clade3').size()); print(d.study_source.value_counts())
print(d[d.study_source.str.startswith('unres')][['t','clade3']].groupby('clade3').size())
