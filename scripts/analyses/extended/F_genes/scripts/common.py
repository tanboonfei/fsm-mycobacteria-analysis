import pandas as pd, csv, re, pickle
R='/path/to/project/'
F=R+'analyses/extended/F_genes/'
SP='/path/to/scratch/'
def load_meta():
    s1=pd.read_excel(R+'Publication/Additional_file_2.xlsx',sheet_name='S1',nrows=862)
    sis=set(open(R+'analyses/sister_tips.txt').read().split())
    cl={}
    for l in open(R+'Publication/scripts/analyses/unassigned_ani/clusters.tsv'):
        if l.startswith('tip'): continue
        t,c=l.split(); cl[t.replace('.fa','')]=c
    def grp(r):
        c=r['clade']
        if c.startswith('FSM'): return 'FSM'
        if c.startswith('RGM'): return 'RGM_sister' if r['tip'] in sis else 'RGM_other'
        if c=='SGM': return 'SGM'
        return 'other'
    s1['group']=s1.apply(grp,axis=1)
    s1['gid']=s1['genome_id_or_accession']
    s1['sal']=s1['salinity_habitat'].fillna('NA')
    s1['species2']=[ (cl.get(g,'novel sp.') if sp=='novel sp.' else sp.replace('Mycobacterium ','')) for g,sp in zip(s1.gid,s1.species)]
    s1['comp']=s1['completeness_used_pct'].astype(float)
    return s1
def targets():
    return list(csv.DictReader(open(F+'targets.tsv'),delimiter='\t'))
GROUPS=['FSM','FSM_low','FSM_medium','FSM_veryhigh','FSM_salNA','RGM_sister','RGM_other','SGM']
def members(s1,g,mincomp=0):
    s=s1[s1.comp>=mincomp]
    if g=='FSM_low': s=s[(s.group=='FSM')&(s.sal=='low salinity')]
    elif g=='FSM_medium': s=s[(s.group=='FSM')&(s.sal=='medium salinity')]
    elif g=='FSM_veryhigh': s=s[(s.group=='FSM')&(s.sal=='very high salinity')]
    elif g=='FSM_salNA': s=s[(s.group=='FSM')&(s.sal=='NA')]
    elif g=='FSM_saline': s=s[(s.group=='FSM')&(s.sal.isin(['medium salinity','very high salinity']))]
    else: s=s[s.group==g]
    return list(s.gid)
