import json, pandas as pd, os
spg=pd.read_csv('species_group_r232.tsv',sep='\t'); G=dict(zip(spg.species,spg.group))
def run_abund(run):
    fn=f'sp_profiles/{run}.json'
    if not os.path.exists(fn): return None
    j=json.load(open(fn))
    if 'profile' not in j: return None
    tot=sum(s for p,s in j['profile']); d={'total_cov':tot}
    for k in ['FSM','RGM','sisterRGM','SGM']: d[k]=0.0
    myco=0; ctrl={'acI_Nanopelagicales':0,'Fonsibacter':0,'Polynucleobacter':0,'Limnohabitans':0,'Flavobacterium':0,'Planctomycetota':0,'Legionella':0,'Rickettsiales':0,'Chlamydiota':0,'Holosporales':0,'Megaera':0,'Caedimonas':0,'Paracaedibacterales':0,'Parachlamydiaceae':0,'Midichloriaceae':0,'Rhabdochlamydiaceae':0}
    for p,s in j['profile']:
        names=p.split(';'); leaf=names[-1]
        if 'Mycobacteriaceae' in p: myco+=s
        if leaf in G and G[leaf] in d: d[G[leaf]]+=s
        if 'Nanopelagicales' in p: ctrl['acI_Nanopelagicales']+=s
        if 'Fonsibacter' in p: ctrl['Fonsibacter']+=s
        if 'Polynucleobacter' in p: ctrl['Polynucleobacter']+=s
        if 'Limnohabitans' in p: ctrl['Limnohabitans']+=s
        if ';Flavobacterium' in p: ctrl['Flavobacterium']+=s
        if 'Planctomycetota' in p: ctrl['Planctomycetota']+=s
        if 'Legionellales' in p or 'Legionella' in p: ctrl['Legionella']+=s
        if 'Rickettsiales' in p: ctrl['Rickettsiales']+=s
        if 'Chlamydiota' in p: ctrl['Chlamydiota']+=s
        if 'Holosporales' in p: ctrl['Holosporales']+=s
        for k in ['Megaera','Caedimonas','Paracaedibacterales','Parachlamydiaceae','Midichloriaceae','Rhabdochlamydiaceae']:
            if k in p: ctrl[k]+=s
    d['Mycobacteriaceae']=myco; d.update(ctrl)
    return d
