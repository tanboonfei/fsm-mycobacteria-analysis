"""Composition of the 59 raw-size and 33 completeness-corrected non-FSM small genomes (Table S13 notes),
the fate table, and the S13c 'cultured_isolates' column, under old and corrected labels."""
import pandas as pd
G='/path/to/project/analyses/extended/G_culture/'
T='/path/to/project/scripts/inputs/tsv/'
s=pd.read_csv(T+'S1_genome_metadata.tsv',sep='\t'); c=pd.read_csv(G+'new_completeness.tsv',sep='\t')
s=s.merge(c[['tip','new_status','completeness_used_new']],on='tip'); s['c3']=s.clade.str.split().str[0]
ISO=['Cultured isolate','NTM (isolate)','Pathogen (isolate)']
s['cls_old']=s.apply(lambda x:'isolate/closed' if (x.culture_status in ISO or x.assembly_status=='closed') else ('SPIRE MAG' if x.origin=='SPIRE MAG' else 'GTDB reference MAG'),axis=1)
s['cls_new']=s.new_status.map(lambda v:'isolate' if v.startswith('isolate') else ('SPIRE MAG' if v=='SPIRE MAG' else 'GTDB MAG'))
ns=s[s.c3.isin(['RGM','SGM'])]
raw=ns[ns.genome_size_Mb<4]
print('raw <4 Mb non-FSM:',len(raw)); print(' old',raw.cls_old.value_counts().to_dict()); print(' new',raw.cls_new.value_counts().to_dict())
out=[]
for lab,comp,cls in [('old','completeness_used_pct','cls_old'),('new','completeness_used_new','cls_new')]:
    b=ns[ns[comp]>=90].copy(); b['csize']=b.genome_size_Mb*100/b[comp]; sm=b[b.csize<4]
    print(lab,'corrected small',len(sm),sm[cls].value_counts().to_dict())
    fate=raw.apply(lambda x:'MAG, completeness <90%' if x[comp]<90 else ('MAG, corrected size >=4 Mb' if x.genome_size_Mb*100/x[comp]>=4 else
          ('isolate/closed (retained)' if x[cls].startswith('isolate') else 'MAG, survives')),axis=1)
    print(lab,'fate of 59',fate.value_counts().to_dict())
    for k,v in fate.value_counts().items(): out.append(dict(labels=lab,table='fate of 59 raw-size small non-FSM genomes',category=k,n=v))
    for k,v in raw[cls].value_counts().items(): out.append(dict(labels=lab,table='59 raw-size small non-FSM genomes',category=k,n=v))
    for k,v in sm[cls].value_counts().items(): out.append(dict(labels=lab,table='33 corrected small non-FSM genomes',category=k,n=v))
pd.DataFrame(out).to_csv(G+'S13_small_genome_composition_old_vs_new.tsv',sep='\t',index=False)
# S13 main table column and S13c cultured_isolates column
st=pd.read_csv(T+'S9_small_genomes.tsv',sep='\t',nrows=59)
m=st.merge(s[['genome_id_or_accession','culture_status','new_status','completeness_used_new','checkm2_this_study_completeness_pct']],on='genome_id_or_accession',how='left')
chg=m[m.culture_status.isin(['MAG / uncultured'])&m.new_status.str.startswith('isolate')]
print('\nS13 rows whose class changes MAG->isolate:'); print(chg[['genome_id_or_accession','species','genome_size_Mb','completeness_pct_or_genome_class','new_status']].to_string())
chg.to_csv(G+'S13_rows_relabelled.tsv',sep='\t',index=False)
sc=pd.read_csv(T+'S9c.tsv',sep='\t',nrows=43)
st_new=s.set_index('genome_id_or_accession').new_status
sc['cultured_isolates_new']=sc.accessions.map(lambda a:sum(st_new[x.strip()].startswith('isolate') for x in a.split(';')))
d=sc[sc.cultured_isolates!=sc.cultured_isolates_new]
print('\nS13c groups whose cultured_isolates count changes:'); print(d[['group','n_genomes','clade','cultured_isolates','cultured_isolates_new','species','accessions']].to_string())
print('S13c total isolates old',sc.cultured_isolates.sum(),'new',sc.cultured_isolates_new.sum(),'groups with >=1 isolate old',(sc.cultured_isolates>0).sum(),'new',(sc.cultured_isolates_new>0).sum())
sc[['group','n_genomes','clade','cultured_isolates','cultured_isolates_new','species','accessions']].to_csv(G+'S13c_cultured_isolates_old_vs_new.tsv',sep='\t',index=False)
