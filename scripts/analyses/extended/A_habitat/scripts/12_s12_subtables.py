"""Assemble compact TSVs intended as new sub-tables of Additional file 2, Table S12."""
import pandas as pd
O='/path/to/project/analyses/extended/A_habitat/'
t=pd.read_csv(O+'T1_studies_per_clade.tsv',sep='\t')
keep=['SPIRE MAG','GTDB MAG (NCBI metagenome-derived)','GTDB genome not flagged metagenome-derived','true MAGs (SPIRE + NCBI metagenome-derived)','S1 biome annotated','S1 salinity annotated','salinity = low salinity','salinity = medium salinity','salinity = high salinity','salinity = very high salinity']
t=t[t.genome_set.isin(keep)][['clade','genome_set','n_genomes','n_species_units','n_studies','n_source_samples','largest_study','largest_study_n','pct_in_largest_study','n_studies_for_50pct']]
t['o']=t.genome_set.map({k:i for i,k in enumerate(keep)}); t=t.sort_values(['o','clade']).drop(columns='o')
t.to_csv(O+'S12f_source_studies_per_clade.tsv',sep='\t',index=False)
x=pd.read_csv(O+'T2_clade_habitat_tests.tsv',sep='\t')
x=x[~x.unit.str.contains('random draw')]
x['interpretable']=~x.note.fillna('').str.contains('DEGENERATE')
x.loc[~x.interpretable,['OR','CI_lo','CI_hi']]=float('nan')
x['note']=x.note.fillna('').str.replace(r'DEGENERATE \(quasi-separation by study; estimates not interpretable\); ','OR not estimable (outcome constant within studies); ',regex=True)
x[['test','unit','FSM_pos','FSM_n','RGM_pos','RGM_n','n_units','OR','CI_lo','CI_hi','p','note']].to_csv(O+'S12g_study_controlled_clade_tests.tsv',sep='\t',index=False)
a=pd.read_csv(O+'T3b_opsin_by_light_class.tsv',sep='\t'); b=pd.read_csv(O+'T3c_opsin_light_tests.tsv',sep='\t')
a.to_csv(O+'S12h_opsin_by_source_light_class.tsv',sep='\t',index=False)
with open(O+'S12h_opsin_by_source_light_class.tsv','a') as fh:
    fh.write('\n'); b.to_csv(fh,sep='\t',index=False)
s=pd.read_csv(O+'T4_sandpiper_habitat_summary.tsv',sep='\t'); s=s[s.run_set.str.startswith('all runs')]
p=pd.read_csv(O+'T4d_within_study_FSM_vs_RGM_prevalence.tsv',sep='\t').set_index('habitat')
cols=['habitat','n_runs','FSM_pos_runs','FSM_prev_pct','FSM_median_pos','FSM_p90_pos','FSM_max','FSM_studies_total','FSM_studies_pos','FSM_studies_pos_pct',
      'RGM_pos_runs','RGM_prev_pct','RGM_median_pos','RGM_p90_pos','RGM_max','RGM_studies_pos','RGM_studies_pos_pct']
s=s[cols].rename(columns={'FSM_studies_total':'studies_total'}).set_index('habitat')
s=s.join(p[['studies_with_either','FSM_higher','RGM_higher','ties','sign_test_p']].add_prefix('within_study_'))
s.reset_index().to_csv(O+'S12i_sandpiper_prevalence_abundance_by_habitat.tsv',sep='\t',index=False)
print(s.to_string())
