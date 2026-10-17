"""Compose Table S11k (per-lineage opsin phylogeny, species-tree mapping, synteny and prevalence) and S11l (per-locus)."""
import pandas as pd
O='/path/to/project/analyses/extended/B_opsin'
S=pd.read_csv(O+'/species/lineage_species_tree_summary.tsv',sep='\t').set_index('trait')
N=pd.read_csv(O+'/prevalence/prevalence_null_per_lineage.tsv',sep='\t').set_index('trait')
R=pd.read_csv(O+'/reconciliation/rf_congruence.tsv',sep='\t').set_index('set')
D=pd.read_csv(O+'/reconciliation/lca_dl_reconciliation.tsv',sep='\t'); D.index=D.lineage+'_'+D.set
C=pd.read_csv(O+'/reconciliation/control/control_rf.tsv',sep='\t').set_index('set')
f=lambda x,d=2: f"{x:.{d}f}"
def pv(p): return f"{p:.2g}"
rows=[]
def add(sec,met,a,b,allv=''): rows.append(dict(Section=sec,Metric=met,**{'Lineage A (actinorhodopsin-related)':a,'Lineage B (DTG-type)':b,'All opsins':allv}))
add('Screen','Opsin loci (FSM / RGM / SGM)','77 / 20 / 14','38 / 0 / 2','115 / 21 (incl. 1 lineage C, M. madagascariense) / 16')
add('Screen','Genomes carrying (FSM / RGM / SGM)','74 / 19 / 14','32 / 0 / 2','103 / 20 / 14 (FSM: 71 A only, 29 B only, 3 both)')
add('Screen','Motif at BR 85/89/96 (mycobacterial loci)','DTE 106, DTG 1 (basal FSM singleton), partial 4','DTG 38, partial 2','')
add('Opsin tree','Mycobacterial sequences in tree (>=150 aa)','107 (74 FSM, 19 RGM, 14 SGM)','38 (36 FSM, 2 SGM)','146 of 152 loci; 424 sequences, 334 columns, LG+F+R8')
add('Opsin tree','Clade and support (SH-aLRT/UFBoot)','114-sequence clade 95/95 (106 mycobacterial + 8 other Actinomycetota opsins: Micrococcales 5, Nakamurella, Rhodococcoides, Williamsia); one basal FSM DTG sequence joins at 79.3/90','100/100','')
add('Opsin tree','FSM opsins monophyletic?','No: of the 73 DTE FSM opsins, 63 form one FSM-only clade (94.4/84) and 10 fall in three clades (5, 3, 2) among RGM/SGM opsins; plus one basal DTG singleton','No: 36 FSM in clades of 28, 5, 2 and 1; the 2 SGM opsins nest inside (97/96 with an FSM opsin)','')
add('Opsin tree','Non-FSM mycobacterial opsins','RGM/SGM interleaved with FSM in the 92-sequence core (99.4/94; 68 FSM, 10 RGM, 14 SGM; 12 SGM in one clade 98.7/100); 9 RGM opsins (93.2/100) sister to 7 opsins of other Actinomycetota (91.5/98)','2 SGM (one pair, 100/100) nested within FSM','M. madagascariense opsin (lineage C) sister to Nakamurella opsin (96.1/100), outside A and B')
add('Opsin tree','Placement','Inside a 139-sequence actinorhodopsin clade (99.4/100) holding all Luna-cluster and acI (Nanopelagicales) actinorhodopsins; position relative to them unresolved (51/47, 69.9/42)','Outside the actinorhodopsin clade; sister (92/92) to DTG-type rhodopsins of Alphaproteobacteria (Hyphomicrobiales, Sphingomonadales, Acetobacterales), Enterobacterales, Pseudomonas and Microlunatus, within a clade with Anabaena sensory rhodopsin (97.5/93)','')
add('Opsin tree','Closest sequence in reference set (median identity)','Williamsia serinedens opsin (Mycobacteriales), 75% (FSM and non-FSM alike)','Rhodopsins of Kushneria/Dictyobacter/Lichenibacterium/Microlunatus, 42% (max 46%)','')
add('Opsin tree','Closest characterised reference (median identity)','acI and Luna-cluster actinorhodopsins: FSM 61% (Ca. Nanopelagicus limnae), 60% (Rhodoluna lacicola); non-FSM 57%','Archaerhodopsin-3 40%, cruxrhodopsin-3 38%, bacteriorhodopsin 35%, Anabaena sensory rhodopsin 31%','')
for k,lab in [('A','A'),('B','B')]:
    pass
def sp(k):
    s=S.loc[k]
    return s
a,b,an=sp('A'),sp('B'),sp('any')
add('Species tree (824-genome ingroup)','FSM carriers monophyletic / MRCA span / maximal all-carrier blocks (largest)',f"No / {int(a.FSM_carrier_MRCA_tips)} of 134 FSM / {a.FSM_pure_blocks} ({a.FSM_largest_block})",f"No / {int(b.FSM_carrier_MRCA_tips)} of 134 / {b.FSM_pure_blocks} ({b.FSM_largest_block})",f"No / {int(an.FSM_carrier_MRCA_tips)} / {an.FSM_pure_blocks} ({an.FSM_largest_block})")
add('Species tree (824-genome ingroup)','Sankoff minimum changes (gain=loss=1; gain=2,loss=1)',f"{a.sankoff_eq}; {a.sankoff_gain2}",f"{b.sankoff_eq}; {b.sankoff_gain2}",f"{an.sankoff_eq}; {an.sankoff_gain2}")
add('Species tree (824-genome ingroup)','Cost: gain on FSM stem / FSM crown absent / present at ingroup root / absent at root (equal weights)',f"{a.cost_stem_gain} / {a.cost_later_gain} / {a.cost_vertical_root} / {a.cost_root_absent}",f"{b.cost_stem_gain} / {b.cost_later_gain} / {b.cost_vertical_root} / {b.cost_root_absent}",f"{an.cost_stem_gain} / {an.cost_later_gain} / {an.cost_vertical_root} / {an.cost_root_absent}")
add('Species tree (824-genome ingroup)','Mk model (AIC ER, ARD); rates 0->1, 1->0',f"ARD ({f(a.AIC_ER,1)}, {f(a.AIC_ARD,1)}); {f(a.q01)}, {f(a.q10)}",f"ARD ({f(b.AIC_ER,1)}, {f(b.AIC_ARD,1)}); {f(b.q01)}, {f(b.q10)}",f"ARD ({f(an.AIC_ER,1)}, {f(an.AIC_ARD,1)}); {f(an.q01)}, {f(an.q10)}")
add('Species tree (824-genome ingroup)','P(present): ingroup root / FSM+sister node / FSM crown',f"{f(a.P_root)} / {f(a.P_FSMsister)} / {f(a.P_FSMcrown)}",f"{f(b.P_root)} / {f(b.P_FSMsister)} / {f(b.P_FSMcrown,3)}",f"{f(an.P_root)} / {f(an.P_FSMsister)} / {f(an.P_FSMcrown)}")
add('Species tree (824-genome ingroup)','lnL conditioned on: stem gain / FSM crown absent / root present / root absent (unconditioned)',f"{f(a.lnL_stem_gain,1)} / {f(a.lnL_later_gain,1)} / {f(a.lnL_vertical_root,1)} / {f(a.lnL_root_absent,1)} ({f(a.lnL,1)})",f"{f(b.lnL_stem_gain,1)} / {f(b.lnL_later_gain,1)} / {f(b.lnL_vertical_root,1)} / {f(b.lnL_root_absent,1)} ({f(b.lnL,1)})",f"{f(an.lnL_stem_gain,1)} / {f(an.lnL_later_gain,1)} / {f(an.lnL_vertical_root,1)} / {f(an.lnL_root_absent,1)} ({f(an.lnL,1)})")
add('Species tree (824-genome ingroup)','Stochastic node-state sampling (1,000): gains (in FSM / elsewhere); losses (in FSM); P(gain on FSM stem branch)',f"{f(a.sm_gains,1)} ({f(a.sm_gains_inFSM,1)} / {f(a.sm_gains_outFSM,1)}); {f(a.sm_losses,1)} ({f(a.sm_losses_inFSM,1)}); {f(a.sm_frac_gain_on_FSM_stem)}",f"{f(b.sm_gains,1)} ({f(b.sm_gains_inFSM,1)} / {f(b.sm_gains_outFSM,1)}); {f(b.sm_losses,1)} ({f(b.sm_losses_inFSM,1)}); {f(b.sm_frac_gain_on_FSM_stem)}",f"{f(an.sm_gains,1)} ({f(an.sm_gains_inFSM,1)} / {f(an.sm_gains_outFSM,1)}); {f(an.sm_losses,1)} ({f(an.sm_losses_inFSM,1)}); {f(an.sm_frac_gain_on_FSM_stem)}")
add('Gene tree vs species tree','nRF (FSM genomes only), label-permutation null minimum; p',f"{R.loc['A_all','nRF']} ({R.loc['A_FSM','nRF']}), {R.loc['A_all','perm_nRF_min']}; p < 0.001",f"{R.loc['B_all','nRF']} ({R.loc['B_FSM','nRF']}), {R.loc['B_all','perm_nRF_min']}; p < 0.001",'')
add('Gene tree vs species tree','Gene-tree splits with SH-aLRT>=80 and UFBoot>=95 that conflict with the species tree (FSM only)',f"{R.loc['A_all','strong_splits_conflicting']} of {R.loc['A_all','strong_splits']} ({R.loc['A_FSM','strong_splits_conflicting']} of {R.loc['A_FSM','strong_splits']})",f"{R.loc['B_all','strong_splits_conflicting']} of {R.loc['B_all','strong_splits']} ({R.loc['B_FSM','strong_splits_conflicting']} of {R.loc['B_FSM','strong_splits']})",'')
add('Gene tree vs species tree','AU test, constraint = species-tree splits with UFBoot>=95: dlnL (p-AU); FSM only','950.2 (7.7e-8); 716.8 (0.0035)','499.4 (4.0e-63); 487.1 (1.6e-65)','')
add('Gene tree vs species tree','Same test, control genes thrS / topA on the same genomes: dlnL; FSM only','641.5 / 228.3; 95.6 / 131.5','10.8 (n.s.) / 37.1; 8.8 (n.s.) / 38.3','per alignment column, opsin 1.5-2.8 vs thrS 0.01-0.94 and topA 0.04-0.25')
add('Gene tree vs species tree','nRF opsin vs thrS vs topA (common genomes)',f"{C.loc['A_all','nRF_opsin']} vs {C.loc['A_all','nRF_thrS']} vs {C.loc['A_all','nRF_topA']}",f"{C.loc['B_all','nRF_opsin']} vs {C.loc['B_all','nRF_thrS']} vs {C.loc['B_all','nRF_topA']}",'')
add('Gene tree vs species tree','AU test, FSM opsins of the lineage forced monophyletic (mycobacterial sequences only): dlnL (p-AU)','10.4 (0.34), not rejected','0.0 (0.50), not rejected; but the outgroup-rooted tree nests the 2 SGM opsins inside FSM (97/96)','')
add('Gene tree vs species tree','LCA duplication-loss reconciliation (D, L) vs species-congruent tree (D, L); FSM only',f"{D.loc['A_all','D']}, {D.loc['A_all','L']} vs 0, {D.loc['A_all','congruent_ref_L']}; {D.loc['A_FSM','D']}, {D.loc['A_FSM','L']} vs 0, {D.loc['A_FSM','congruent_ref_L']}",f"{D.loc['B_all','D']}, {D.loc['B_all','L']} vs 0, {D.loc['B_all','congruent_ref_L']}; {D.loc['B_FSM','D']}, {D.loc['B_FSM','L']} vs 0, {D.loc['B_FSM','congruent_ref_L']}",'')
add('Gene tree vs species tree','GeneRax 2.1.3 reconciliation lnL, fixed ML gene trees: UndatedDL vs UndatedDTL','-853.2 vs -620.3','-287.3 vs -207.9','sum -1140.5 vs -828.2 (one extra parameter); event counts not obtained (program fault)')
add('Neighbourhood (Bakta, +-5 genes)','Blh gene adjacent to opsin','FSM 0 of 73; non-FSM 26 of 34 (27 within 5 genes)','FSM 0 of 36; SGM 0 of 2','M. madagascariense: adjacent')
add('Neighbourhood (Bakta, +-5 genes)','Commonest flanks','FSM: TetR regulator (62%), SDR, membrane protein, GT2 (TET site); 15 FSM loci at the alpha-mannosidase/Rv0647c/AarF site; non-FSM: Blh (82%) and an unrelated TetR-type regulator (68%)','FSM: alpha-mannosidase/Rv0647c/AarF kinase site (24 of 25 with context), shared with 15 lineage-A FSM loci; SGM: F420 oxidoreductase, CzcO','')
add('Neighbourhood (Bakta, +-5 genes)','Mean Jaccard of neighbour families (pairs sharing any)','FSM-A vs non-FSM-A 0.007 (5%); within FSM-A 0.24 (42%); within non-FSM-A 0.15 (93%)','FSM-B vs non-FSM-A 0.0003 (0.5%); FSM-A vs FSM-B 0.21 (24%)','')
add('Neighbourhood (Bakta, +-5 genes)','Blh anywhere in genome (tBLASTn, mycobacterial Blh queries)','FSM genomes with A only: 40 of 71; non-FSM carriers: 30 of 34 (0 of 656 non-carriers)','FSM genomes with B only: 0 of 29','FSM: 44 of 134, never on an opsin contig; opsin+ 41/103 vs opsin- 3/31, Fisher p = 0.0019')
for k in ['A','B','any']:
    pass
add('Prevalence (FSM vs rest)','phyloglm FSM vs RGM+SGM (824): beta (SE), p',f"{f(a.glm_ingroup_beta)} ({f(a.glm_ingroup_se)}), {pv(a.glm_ingroup_p)}",f"{f(b.glm_ingroup_beta)} ({f(b.glm_ingroup_se)}), {pv(b.glm_ingroup_p)}",f"{f(an.glm_ingroup_beta)} ({f(an.glm_ingroup_se)}), {pv(an.glm_ingroup_p)}; Table S1 calls: {f(S.loc['any_S1','glm_ingroup_beta'])}, {pv(S.loc['any_S1','glm_ingroup_p'])}")
add('Prevalence (FSM vs rest)','phyloglm FSM vs RGM (525): beta (SE), p',f"{f(a.glm_FSMRGM_beta)} ({f(a.glm_FSMRGM_se)}), {pv(a.glm_FSMRGM_p)}",f"{f(b.glm_FSMRGM_beta)} ({f(b.glm_FSMRGM_se)}), {pv(b.glm_FSMRGM_p)} (no RGM carrier; quasi-separation)",f"{f(an.glm_FSMRGM_beta)} ({f(an.glm_FSMRGM_se)}), {pv(an.glm_FSMRGM_p)}; Table S1 calls: {f(S.loc['any_S1','glm_FSMRGM_beta'])}, {pv(S.loc['any_S1','glm_FSMRGM_p'])}")
add('Prevalence (FSM vs rest)','Prevalence-conditioned ARD Mk null, FSM vs RGM: observed OR; null 95th percentile; p',f"{f(N.loc['A','observed_OR'],1)}; {f(N.loc['A','null_OR_p95'],1)}; {pv(N.loc['A','p'])}",f"{f(N.loc['B','observed_OR'],1)}; {f(N.loc['B','null_OR_p95'],1)}; {pv(N.loc['B','p'])}",f"{f(N.loc['any','observed_OR'],1)}; {f(N.loc['any','null_OR_p95'],1)}; {pv(N.loc['any','p'])}; Table S1 calls: {f(N.loc['any_S1','observed_OR'],1)}; {f(N.loc['any_S1','null_OR_p95'],1)}; {pv(N.loc['any_S1','p'])}")
T=pd.DataFrame(rows); T.to_csv(O+'/Table_S11k_opsin_lineages.tsv',sep='\t',index=False)
# per-locus table
L=pd.read_csv(O+'/screen/opsin_loci_lineages.tsv',sep='\t')
nb=pd.read_csv(O+'/synteny/opsin_neighbourhoods.tsv',sep='\t')[['locus','blh_adjacent','g-1','g+1']]
cx=pd.read_csv(O+'/synteny/opsin_context.tsv',sep='\t')
br=pd.read_csv(O+'/tree/best_reference_per_locus.tsv',sep='\t')[['locus','ref','ref_label','ref_org','ident']]
X=L.merge(nb,on='locus',how='left').merge(cx,on='locus',how='left').merge(br,on='locus',how='left')
cols=['locus','genome','clade','species','origin','contig','locus_start','locus_end','dist_to_contig_end','aa_len','best_query','best_bits','motif','lineage','sublineage','lineage_basis','in_tree','old_id','old_lineage','context','blh_adjacent','g-1','g+1','ref','ref_org','ident']
X=X[cols].rename(columns={'old_id':'published_tree_id','old_lineage':'published_tree_group','ref':'closest_reference_accession','ref_org':'closest_reference_organism','ident':'closest_reference_identity_pct','g-1':'upstream_gene','g+1':'downstream_gene','context':'neighbourhood_context'})
X.to_csv(O+'/Table_S11l_opsin_loci.tsv',sep='\t',index=False)
print(T.to_string()[:3000]); print(X.shape)
