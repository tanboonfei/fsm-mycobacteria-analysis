import csv
S1='../../../Publication/scripts/inputs/tsv/S1_genome_metadata.tsv'
sister=set(l.strip() for l in open('../../sister_tips.txt'))
def load():
    G={}
    for r in csv.DictReader(open(S1),delimiter='\t'):
        gid=r['genome_id_or_accession']; c=r['clade']; p=r['pathogenicity']
        if c.startswith('FSM'): g='FSM'
        elif c.startswith('RGM'): g='RGM_sister' if r['tip'] in sister else 'RGM_other'
        elif c=='SGM':
            if r['species'] in ('Mycobacterium leprae','Mycobacterium lepromatosis'): g='Leprosy'
            elif p=='obligate_pathogen': g='SGM_obligate'
            elif p=='opportunistic_pathogen': g='NTM_opportunistic'
            else: g='SGM_env'
        elif c.startswith('Outside'): g='Outside'
        else: g='Outgroup'
        if g!='Leprosy' and p=='opportunistic_pathogen' and c.startswith('RGM'): g='NTM_opportunistic'
        try: comp=float(r['completeness_used_pct'])
        except: comp=100.0
        G[gid]=dict(group=g,tip=r['tip'],species=r['species'],size=float(r['genome_size_Mb']),comp=comp,origin=r['origin'])
    return G
if __name__=='__main__':
    import collections
    G=load(); print(collections.Counter(v['group'] for v in G.values()))
