r"""Corrected culture status for all 862 S1 genomes.
Rules (applied in order):
 1. origin == 'SPIRE MAG'                                   -> SPIRE MAG
 2. GTDB reference with any metagenome evidence             -> GTDB MAG
      (a) NCBI exclfromrefseq contains 'derived from metagenome', OR (b) BioSample package MIMAG.*,
      OR (c) the BioSample record mentions 'metagenom' anywhere (title, organism, attributes),
      OR (d) the BioSample title / sample name / isolate name is a MAG or bin name
      (\bMAG\d*\b, bin.N, metabat/concoct/maxbin), OR (e) the BioProject title/description
      mentions metagenomics AND the assembly carries no strain name AND the BioSample is not a
      pure-culture/pathogen/MIGS-bacterium record
 3. GTDB reference with BioSample package MISAG.*           -> SAG (none found)
 4. otherwise (no metagenome / single-cell evidence)        -> isolate
      closed/complete if NCBI assembly level is 'Complete Genome' or the S1 assembly is a single
      contig (S1 assembly_status == 'closed'); otherwise draft.
 'sample_type = single cell' in Microbe.1.0 BioSamples of named, type or deposited strains is a
 known mis-use of the field and is not taken as SAG evidence (flagged in evidence).
"""
import pandas as pd, json, re, xml.etree.ElementTree as ET
G='/path/to/project/analyses/extended/G_culture/'
S1='/path/to/project/scripts/inputs/tsv/S1_genome_metadata.tsv'
s=pd.read_csv(S1,sep='\t')
d=pd.read_csv(G+'ncbi_evidence_gtdb.tsv',sep='\t').fillna('').set_index('tip')
R=json.load(open(G+'cache/assembly_esummary.json')); B=json.load(open(G+'cache/biosample_xml.json')); BP=json.load(open(G+'cache/bioproject_esummary.json'))
def extra(a):
    e=R[a]; el=ET.fromstring(B[e['biosampleaccn']]); ev=[]
    names=[el.findtext('Description/Title') or '', e['biosource'].get('isolate','')]+[at.text or '' for at in el.iter('Attribute')
           if (at.get('harmonized_name') or at.get('attribute_name')) in ('sample_name','isolate','Submitter Id')]
    hit=[n for n in names if re.search(r'\bMAG\d*\b|\bbin[._-]?\d|metabat|concoct|maxbin',n,re.I)]
    if hit: ev.append(f'BioSample/isolate name is a MAG/bin name ("{hit[0][:60]}")')
    strain=any(i['sub_type']=='strain' for i in (e['biosource'].get('infraspecieslist') or []))
    pkg=el.findtext('Package') or ''; pure='pure culture' in B[e['biosampleaccn']].lower()
    for p in e.get('gb_bioprojects',[]):
        v=BP[p['bioprojectaccn']]; txt=(v['title']+' '+v['description']).lower()
        if 'metagenom' in txt and not strain and not pure and not pkg.startswith(('Pathogen','MIGS.ba')):
            ev.append(f'BioProject {p["bioprojectaccn"]} is a metagenomic study ("{v["title"][:70]}") and the assembly has no strain name')
    return ev
out=[]
for _,x in s.iterrows():
    r=dict(tip=x.tip, accession=x.genome_id_or_accession, clade=x.clade, origin=x.origin,
           old_status=x.culture_status, s1_assembly_status=x.assembly_status, n_contigs=x.n_contigs)
    if x.origin=='SPIRE MAG':
        r.update(new_status='SPIRE MAG', ncbi_assembly_level='', evidence='SPIRE v1 metagenome bin (this study)')
    else:
        e=d.loc[x.tip]; ev=[]; a=x.genome_id_or_accession
        meta='derived from metagenome' in e.excl
        if meta: ev.append("NCBI excluded-from-RefSeq: 'derived from metagenome'")
        if e.bs_package.startswith('MIMAG'): ev.append(f'BioSample package {e.bs_package}')
        if e.kw_metagenom and not meta: ev.append('BioSample text mentions metagenome')
        ev+= [z for z in extra(a) if not ev]   # (d)/(e) only decide when (a)-(c) are silent
        sag=e.bs_package.startswith('MISAG')
        r['ncbi_assembly_level']=e.ncbi_assembly_level
        if ev:
            r['new_status']='GTDB MAG'
        elif sag:
            r['new_status']='SAG'; ev.append(f'BioSample package {e.bs_package}')
        else:
            closed=(e.ncbi_assembly_level=='Complete Genome') or (x.assembly_status=='closed')
            r['new_status']='isolate (closed/complete)' if closed else 'isolate (cultured, draft)'
            ev.append('no metagenome/MIMAG/MISAG evidence in NCBI assembly or BioSample')
            if e.strain: ev.append(f'strain {e.strain}')
            if e.fromtype: ev.append(f'type material ({e.fromtype})')
            ev.append(f'BioSample package {e.bs_package}')
            if e.bs_isolation_source and e.bs_isolation_source.lower() not in ('missing','not collected','not applicable'):
                ev.append(f'isolation_source "{e.bs_isolation_source}"')
            ev.append(f'NCBI level {e.ncbi_assembly_level}')
            if e.kw_single_cell: ev.append("BioSample sample_type='single cell' (strain record; not a SAG)")
            if 'Candidatus' in e.ncbi_organism: ev.append('BORDERLINE: Candidatus name; kept as isolate (BioSample/BioProject describe a strain or pure culture)')
        if e.excl and not meta: ev.append(f'RefSeq exclusion: {e.excl}')
        ev.append(f'{e.biosample}; {e.bioproject}; {e.submitter}')
        r['evidence']=' | '.join(ev)
    out.append(r)
o=pd.DataFrame(out)
o['status_changed']=o.apply(lambda r: not ((r.old_status=='MAG / uncultured' and r.new_status in('SPIRE MAG','GTDB MAG')) or
     (r.old_status in('Cultured isolate','NTM (isolate)','Pathogen (isolate)') and r.new_status.startswith('isolate'))),axis=1)
o['status_changed']=o.status_changed.map({True:'yes',False:'no'}); o.loc[o.old_status=='Outgroup','status_changed']='n/a (old label was clade, not status)'
cols=['tip','accession','clade','origin','old_status','new_status','status_changed','s1_assembly_status','n_contigs','ncbi_assembly_level','evidence']
o[cols].to_csv(G+'corrected_culture_status.tsv',sep='\t',index=False)
pd.set_option('display.width',250)
print(pd.crosstab([o.old_status],o.new_status,margins=True))
print(pd.crosstab([o.clade],o.new_status,margins=True))
