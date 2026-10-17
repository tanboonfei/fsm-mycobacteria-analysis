import re, csv, os
D='t1'
desc={'ML_LG+G4_r1':'Unconstrained ML tree (LG+G4; FSM sister to 42-genome clade)','ML_PMSF_C60':'Unconstrained ML tree, LG+C60+F+R4 PMSF (FSM sister to 42-genome clade)',
'ML_PMSF_C20_published':'Published LG+C20+F+G PMSF tree (FSM sister to 42-genome clade)','ML_fast30_removed':'ML tree from 30% fast-site-removed alignment (FSM inside 42-genome clade, sister to M. anyangense group)',
'H1_anyangense4':'FSM sister to M. anyangense group (4; the 30%-removal placement)','H2_S42_rhodesiae30':'FSM sister to the M. rhodesiae subclade (30) of the 42-genome clade',
'H3_S42_sp020438785_12':'FSM sister to the sp020438785/sp020438865 subclade (12) of the 42-genome clade','H4_palauense4':'FSM sister to the M. palauense/M. fallax group (4)',
'H5_S42plusPalauense46':'FSM sister to (42-genome clade + palauense group) (46)','H6_RGM_chubuense_phlei19':'FSM sister to the M. chubuense/M. phlei RGM clade (19)',
'H7_RGM_fortuitum_neoaurum18':'FSM sister to the M. fortuitum/M. neoaurum/M. smegmatis RGM clade (18)','H8_SGM30':'FSM sister to the slow growers (30)','H9_allRGM83':'FSM sister to all other RGM (83)'}
def au(tag):
    f=f'{D}/au_{tag}.iqtree'
    if not os.path.exists(f): return {}
    names=open(f'{D}/au_cand_names_{tag}.txt').read().split()
    rows=[l.split() for l in open(f) if re.match(r'^ +\d+ +-\d',l)]
    return {n:(r[2],r[-2]) for n,r in zip(names,rows)}
A={t:au(t) for t in ('lgg4b','lgfr7','c60')}
if not A['lgg4b']: A['lgg4b']=au('lgg4')
names=[n for n in desc if any(n in A[t] for t in A)]
with open('tables/S4h_AU_FSM_placement.tsv','w') as w:
    w.write('Hypothesis\tTopology\tdlnL_LG+G4\tpAU_LG+G4\tdlnL_LG+F+R7\tpAU_LG+F+R7\tdlnL_LG+C60+F+R4_PMSF\tpAU_LG+C60+F+R4_PMSF\n')
    for n in names:
        v=[]
        for t in ('lgg4b','lgfr7','c60'):
            v+=list(A[t].get(n,('NA','NA')))
        w.write('\t'.join([n,desc[n]]+v)+'\n')
print(open('tables/S4h_AU_FSM_placement.tsv').read())

# ---------- Task 3: shift detection ----------
def rd(f): return list(csv.DictReader(open(f),delimiter='\t')) if os.path.exists(f) else []
meta={'l1ou_824_log_pBIC_chronos_l1':('l1ou','824 ingroup','log Mb','chronos (λ=1) ultrametric','pBIC'),
'l1ou_824_log_BIC_chronos_l1':('l1ou','824 ingroup','log Mb','chronos','BIC'),
'l1ou_824_raw_pBIC_chronos_l1':('l1ou','824 ingroup','Mb','chronos','pBIC'),
'l1ou_525_log_pBIC_chronos_l1':('l1ou','525 FSM+RGM','log Mb','chronos','pBIC'),
'l1ou_824_log_pBIC_chronos_l1_derep':('l1ou','648 (one genome per species)','log Mb','chronos','pBIC'),
'l1ou_824_log_pBIC_extend':('l1ou','824 ingroup','log Mb','phylogram with terminals extended (force.ultrametric)','pBIC')}
sh=rd('t3/l1ou_shifts_all.tsv'); reg=rd('t3/l1ou_regimes.tsv'); drop=rd('t3/l1ou_drop_fsm_shift.tsv')
boot={}
for f in os.listdir('t3'):
    if '_boot' in f and f.endswith('.tsv'):
        run=f.split('_boot')[0]
        for x in rd('t3/'+f): boot[(run,x['edge'])]=x['bootstrap_detection']
rows=[]
for run,(m,data,scale,tree,crit) in meta.items():
    s=[x for x in sh if x['run']==run]
    if not s: continue
    f=[x for x in s if x['label']=='FSM stem']; l=[x for x in s if x['label'].startswith('leprosy')]
    rg={x['group']:x for x in reg if x['run']==run}
    d=[x for x in drop if x['run']==run and x['criterion']==crit]
    ratio=lambda x: x['tip_mean_ratio'] if x['tip_mean_ratio'] not in ('NA','') else ''
    rows.append([m,data,scale,tree,crit,s[0]['nShifts'],s[0]['alpha_unit_height'],
      'yes' if f else 'no', f[0]['shift_in_tip_mean'] if f else '', ratio(f[0]) if f else '',
      rg.get('FSM',{}).get('optimum_Mb_modal',''), rg.get('sister42',{}).get('optimum_Mb_modal',''),
      d[0]['delta'] if d else '', boot.get((run,f[0]['edge']),'') if f else '',
      'yes' if l else 'no', l[0]['shift_in_tip_mean'] if l else '',
      str(sum(1 for x in s if x['n_tips']=='1'))])
pem={'pem_824_log':('824','chronos','scOU','no'),'pem_824_log_BMphylogram':('824','phylogram','BM','no'),'pem_824_log_BM_all_me':('824','phylogram','BM','yes'),
'pem_824_log_BM_derep_me':('648','phylogram','BM','yes'),'pem_824_log_scOU_all_me':('824','chronos','scOU','yes'),'pem_824_log_scOU_derep_me':('648','chronos','scOU','yes')}
for x in rd('t3/pem_selection_summary.tsv'):
    n,tree,proc,me=pem[x['run']]
    ss=[y for y in rd('t3/'+x['run']+'_shifts.tsv') if y['method']==x['method']]
    fe=[y for y in ss if y['label']=='FSM stem']; le=[y for y in ss if y['label'].startswith('leprosy')]
    rows.append(['PhylogeneticEM ('+proc+')',n+(' (one genome per species)' if n=='648' else ' ingroup'),'log Mb',tree+(' + tip-error term' if me=='yes' else ''),x['method'],x['K_selected'],x['alpha'],
      'yes' if x['FSM_stem_shift']=='TRUE' else 'no',fe[0]['effective_tip_shift'] if fe else '',fe[0]['effective_ratio'] if fe else '','','','','', 'yes' if x['leprosy_clade_shift']=='TRUE' else 'no',le[0]['effective_tip_shift'] if le else '',x['n_single_tip_shifts']])
hdr=['Method','Tips','Trait scale','Tree','Selection criterion','Shifts selected','alpha (per tree height)','Shift on FSM stem','FSM shift (tip-mean, trait units)','FSM shift ratio (log scale)','FSM optimum (Mb)','42-genome sister-clade optimum (Mb)','Criterion cost of removing FSM shift','FSM-shift bootstrap detection rate','Shift on leprosy clade','Leprosy shift (tip-mean)','Single-tip shifts']
with open('tables/S4j_unsupervised_shifts_summary.tsv','w') as w:
    w.write('\t'.join(hdr)+'\n'); [w.write('\t'.join(map(str,r))+'\n') for r in rows]
# detail of l1ou shifts for the main run
with open('tables/S4k_l1ou_shift_configuration_main.tsv','w') as w:
    cols=['n_tips','composition','label','mean_size_Mb_below','shift_in_optimum','shift_in_tip_mean','tip_mean_ratio','species_examples']
    w.write('\t'.join(['run']+cols+['bootstrap_detection'])+'\n')
    for x in sh:
        if x['run'] in ('l1ou_824_log_pBIC_chronos_l1','l1ou_824_log_pBIC_chronos_l1_derep'):
            w.write('\t'.join([x['run']]+[x[c] for c in cols]+[boot.get((x['run'],x['edge']),'')])+'\n')
print(open('tables/S4j_unsupervised_shifts_summary.tsv').read())

# ---------- Task 2: concordance ----------
def cf(f): return {x['bipartition']:x for x in rd(f)}
g=cf('t2/concord_gcf_clades.tsv'); s=cf('t2/concord_scfl_clades.tsv'); gb=cf('t2/gcf_on_bac120tree_clades.tsv')
sup={x['bipartition']:x for x in rd('support_concat299_partitioned.tsv')}
sup60={x['bipartition']:x for x in rd('support_PMSF_C60.tsv')}
freq={x['bipartition']:x for x in rd('t2/genetree_bipartition_freq.tsv')}
order=['FSM','42-genome sister clade','FSM + 42-genome clade','FSM + 42 + palauense group','S42 subclade (rhodesiae, 30)','S42 subclade (sp020438785/865, 12)']
supkey={'42-genome sister clade':'42-genome clade','FSM + 42 + palauense group':'FSM + 42 + palauense','S42 subclade (rhodesiae, 30)':'rhodesiae subclade','S42 subclade (sp020438785/865, 12)':'sp020438785 subclade'}
with open('tables/S4i_orthologue_tree_concordance.tsv','w') as w:
    w.write('Bipartition\tTaxa\tConcatenated 299-gene tree: recovered / SH-aLRT\tgCF % (gDF1/gDF2/gDFP), N decisive\tsCFL % (sDF1/sDF2), mean informative sites\tgCF on bac120 LG+G4 tree\tbac120 LG+C60+F+R4 PMSF: SH-aLRT/UFBoot\n')
    for b in order:
        k=supkey.get(b,b); a=sup.get(k,{}); c=sup60.get(k,{})
        G=g.get(b,{}); S=s.get(b,{}); B=gb.get(b,{})
        w.write('\t'.join([b,G.get('n_taxa',''),('yes / '+a['support']) if a.get('recovered')=='TRUE' else 'no',
          f"{G.get('gCF')} ({G.get('gDF1')}/{G.get('gDF2')}/{G.get('gDFP')}), {G.get('gN')}", f"{S.get('sCF')} ({S.get('sDF1')}/{S.get('sDF2')}), {S.get('sN')}",
          B.get('gCF',''), c.get('support','') if c.get('recovered')=='TRUE' else 'not recovered'])+'\n')
with open('tables/S4i2_gene_tree_competing_bipartitions.tsv','w') as w:
    w.write('Bipartition\tGene trees containing it\tGene trees (n)\tPercent\n')
    for x in rd('t2/genetree_bipartition_freq.tsv'): w.write(f"{x['bipartition']}\t{x['n_present']}\t{x['n_gene_trees_informative']}\t{x['percent']}\n")
print(open('tables/S4i_orthologue_tree_concordance.tsv').read())
