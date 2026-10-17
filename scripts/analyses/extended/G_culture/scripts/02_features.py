"""Per-genome NCBI evidence table for all GTDB-reference genomes in S1."""
import json, re, pandas as pd, xml.etree.ElementTree as ET
G='/path/to/project/analyses/extended/G_culture/'
S1='/path/to/project/scripts/inputs/tsv/S1_genome_metadata.tsv'
s=pd.read_csv(S1,sep='\t')
r=json.load(open(G+'cache/assembly_esummary.json')); b=json.load(open(G+'cache/biosample_xml.json'))
rows=[]
for _,x in s[s.origin=='GTDB reference'].iterrows():
    a=x.genome_id_or_accession; e=r[a]
    bs=e.get('biosampleaccn'); el=ET.fromstring(b[bs])
    attrs={at.get('harmonized_name') or at.get('attribute_name'):(at.text or '') for at in el.iter('Attribute')}
    pkg=(el.findtext('Package') or '')
    org=el.find('Description/Organism'); orgname=(org.get('taxonomy_name') if org is not None else '') or ''
    title=el.findtext('Description/Title') or ''
    ids=' '.join((i.text or '') for i in el.iter('Id'))
    blob=' '.join([title,orgname,ids]+[f'{k}={v}' for k,v in attrs.items()]).lower()
    inf=e['biosource'].get('infraspecieslist') or []
    rows.append(dict(tip=x.tip, accession=a, clade=x.clade, old_status=x.culture_status, s1_assembly_status=x.assembly_status,
        n_contigs=x.n_contigs, species=x.species, ncbi_organism=e['organism'], ncbi_assembly_level=e['assemblystatus'],
        excl=';'.join(e['exclfromrefseq']), refseq_category=e['refseq_category'], fromtype=e['fromtype'],
        strain=';'.join(i['sub_value'] for i in inf if i['sub_type']=='strain'), isolate_field=e['biosource'].get('isolate',''),
        partial=e['partialgenomerepresentation'], bioproject=';'.join(p['bioprojectaccn'] for p in e['gb_bioprojects']),
        submitter=e['submitterorganization'], biosample=bs, bs_package=pkg, bs_organism=orgname,
        bs_isolation_source=attrs.get('isolation_source',''), bs_env_medium=attrs.get('env_medium',''),
        bs_host=attrs.get('host',''),
        bs_metagenome_source=attrs.get('metagenome_source',''), bs_derived_from=attrs.get('derived_from',''),
        bs_sample_type=attrs.get('sample_type',''), bs_attr_names=';'.join(attrs),
        kw_metagenom=('metagenom' in blob), kw_mag=bool(re.search(r'metagenome[- ]assembled|\bmag\b|\bbin[_ .-]?\d',blob)),
        kw_single_cell=bool(re.search(r'single[- ]cell|\bsag\b|misag',blob)),
        kw_culture=bool(re.search(r'pure culture|cultur|isolat',blob))))
d=pd.DataFrame(rows); d.to_csv(G+'ncbi_evidence_gtdb.tsv',sep='\t',index=False); print(d.shape)
