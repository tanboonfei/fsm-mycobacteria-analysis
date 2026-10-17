"""Table 1 composition counts, old (published S1 labels) vs corrected."""
import pandas as pd
G='/path/to/project/analyses/extended/G_culture/'
s=pd.read_csv('/path/to/project/scripts/inputs/tsv/S1_genome_metadata.tsv',sep='\t')
c=pd.read_csv(G+'corrected_culture_status.tsv',sep='\t'); s=s.merge(c[['tip','new_status']],on='tip')
s['c3']=s.clade.str.split().str[0]
ISO=['Cultured isolate','NTM (isolate)','Pathogen (isolate)']
rows=[('FSM',['FSM']),('RGM',['RGM']),('SGM',['SGM']),('Ingroup (crown clade)',['FSM','RGM','SGM']),('Outside crown clade',['Outside']),('Outgroup (Tsukamurella)',['Outgroup']),('Total',None)]
out=[]
for name,cl in rows:
    x=s if cl is None else s[s.c3.isin(cl)]
    out.append({'Genome set':name,'Genomes':len(x),
        'Cultured isolates (old S1 label)':int(x.culture_status.isin(ISO).sum()),
        'Cultured isolates (corrected)':int(x.new_status.str.startswith('isolate').sum()),
        '  of which closed/complete (NCBI Complete Genome or single contig)':int((x.new_status=='isolate (closed/complete)').sum()),
        'SPIRE MAGs':int((x.new_status=='SPIRE MAG').sum()),
        'Uncultured reference genomes (old S1 label)':int(((x.origin=='GTDB reference')&(x.culture_status=='MAG / uncultured')).sum()),
        'Uncultured reference genomes = GTDB MAGs (corrected)':int((x.new_status=='GTDB MAG').sum()),
        'Closed assemblies (single contig; unchanged)':int((x.assembly_status=='closed').sum()),
        'Closed assemblies that are MAGs (corrected)':int(((x.assembly_status=='closed')&~x.new_status.str.startswith('isolate')).sum())})
o=pd.DataFrame(out); o.to_csv(G+'Table1_old_vs_corrected.tsv',sep='\t',index=False)
print(o.T.to_string())
