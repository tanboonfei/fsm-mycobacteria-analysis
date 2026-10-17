"""completeness_used under corrected labels: isolates (draft or closed) = 100 by definition (the
manuscript's rule); MAGs = CheckM2 v1.0.1 (this study; the 4 newly recognised GTDB MAGs scored here)."""
import pandas as pd
G='/path/to/project/analyses/extended/G_culture/'
s=pd.read_csv('/path/to/project/scripts/inputs/tsv/S1_genome_metadata.tsv',sep='\t')
c=pd.read_csv(G+'corrected_culture_status.tsv',sep='\t').set_index('tip')
q=pd.concat([pd.read_csv(G+f'{d}/quality_report.tsv',sep='\t') for d in ('checkm2_out','checkm2_out2')]).set_index('Name')
s['new_status']=s.tip.map(c.new_status); s['c3']=s.clade.str.split().str[0]
rows=[]
for _,x in s.iterrows():
    if x.c3 not in ('FSM','RGM','SGM'): v,src,ct=float('nan'),'not in ingroup; not scored',float('nan')
    elif x.new_status.startswith('isolate'): v,src,ct=100.0,'complete by definition (isolate)',float('nan')
    elif pd.notna(x.checkm2_this_study_completeness_pct): v,src,ct=x.checkm2_this_study_completeness_pct,'CheckM2 v1.0.1 (this study)',x.checkm2_this_study_contamination_pct
    else:
        a=x.genome_id_or_accession; v,src,ct=q.loc[a,'Completeness'],'CheckM2 v1.0.1 (this study)',q.loc[a,'Contamination']
    rows.append(dict(tip=x.tip,new_status=x.new_status,completeness_used_old=x.completeness_used_pct,completeness_source_old=x.completeness_source,
                     completeness_used_new=v,completeness_source_new=src,checkm2_contamination_new=ct))
o=pd.DataFrame(rows); o.to_csv(G+'new_completeness.tsv',sep='\t',index=False)
i=o[o.completeness_source_new!='not in ingroup; not scored']
print(i.completeness_source_new.value_counts()); print(i.completeness_source_old.value_counts())
print('ingroup <90 old',(i.completeness_used_old<90).sum(),'new',(i.completeness_used_new<90).sum())
ch=i[i.completeness_used_old!=i.completeness_used_new]
print('changed completeness:',len(ch)); print(ch.describe())
