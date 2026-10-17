# Build the TSV sub-tables from all outputs present
import json, os, re, glob, sys, statistics as st
sys.path.insert(0, os.path.dirname(__file__)); from relax_summary import summ
D='/path/to/project/analyses/extended/D_relax/'
T=D+'tables/'; os.makedirs(T,exist_ok=True)
cl=json.load(open(D+'controls/control_clades.json')); occ=json.load(open(D+'controls/occupancy.json'))
f=lambda x,n=4: ('' if x is None else (f'{x:.{n}g}' if isinstance(x,float) else str(x)))
# --- S9m control clades
rows=[]; Ks=[]
hdr=['run','group','n_genomes','species','placement_note','n_test_branches','n_reference_fitted','n_unclassified','summed_test_length_MG94','K','LRT','p','omega_test_MG94','omega_reference_MG94','lnL_alternative','RELAX_test_distribution','RELAX_reference_distribution','genomes','mean_OG_occupancy_of_52']
graft={'CC03_RGM':'grafted at bac120 position (not monophyletic in 52-gene ML tree)','CC11_SGM':'grafted at bac120 position (not monophyletic in 52-gene ML tree)'}
for c in cl:
    if c in ('CC08_RGM','CC09_RGM'): continue
    p=D+f'controls/{c}.relax.json'
    if not os.path.exists(p) or os.path.getsize(p)==0: continue
    s=summ(p); Ks.append((c,s['K'],s['len_test'],s['n_test']))
    rows.append([c,c[-3:],len(cl[c]['tips']),'; '.join(cl[c]['species']),graft.get(c,'constrained ML placement'),s['n_test'],s['n_ref'],s['n_uncl'],f(s['len_test']),f(s['K'],5),f(s['LRT']),f(s['p']),f(s['omega_test']),f(s['omega_ref']),f(s['lnL_alt'],10),s['alt_test'],s['alt_ref'],','.join(cl[c]['tips']),f(st.mean(occ[c].values()),3)])
prim={'FSM (published run)':D+'inputs/common_reference_FSM_tips_stem.json','FSM (rerun, this study)':D+'boot/real/FSM_tips_stem.json','Leprosy clade (published run)':D+'inputs/common_reference_LEP_tips_stem.json'}
for n,p in prim.items():
    s=summ(p); rows.append([n,'reduced','14' if 'FSM' in n else '3','','primary 37-taxon tree',s['n_test'],s['n_ref'],s['n_uncl'],f(s['len_test']),f(s['K'],5),f(s['LRT']),f(s['p']),f(s['omega_test']),f(s['omega_ref']),f(s['lnL_alt'],10),s['alt_test'],s['alt_ref'],'',''])
with open(T+'S9m_control_clade_null.tsv','w') as o:
    o.write('\t'.join(hdr)+'\n'); [o.write('\t'.join(map(str,r))+'\n') for r in rows]
k=sorted(x[1] for x in Ks)
summary={'n_controls':len(k),'K_min':min(k) if k else None,'K_max':max(k) if k else None,'K_median':st.median(k) if k else None}
for name,val in [('FSM_published',1.043835267381168),('FSM_rerun',summ(prim['FSM (rerun, this study)'])['K']),('LEP_published',0.7579124)]:
    summary[name+'_n_controls_le']=sum(1 for x in k if x<=val); summary[name+'_empirical_lower_tail_p']=(1+sum(1 for x in k if x<=val))/(len(k)+1)
json.dump({'controls':Ks,'summary':summary},open(T+'control_null_summary.json','w'),indent=1); print(summary)
# --- bootstrap
res=[]
for p in sorted(glob.glob(D+'boot/fast2/rep*.log'), key=lambda x:int(re.findall(r'rep(\d+)',x)[0])):
    for l in open(p):
        if l.startswith('RESULT'):
            a=l.rstrip().split('\t'); res.append((os.path.basename(p)[:-4],float(a[6]),float(a[8]),float(a[2]),float(a[4])))
obs=None
for l in open(D+'boot/fast2/REAL.log') if os.path.exists(D+'boot/fast2/REAL.log') else []:
    if l.startswith('RESULT'): a=l.rstrip().split('\t'); obs=(float(a[6]),float(a[8]))
with open(T+'S9p_parametric_bootstrap.tsv','w') as o:
    o.write('replicate\tLRT_conditional\tK_conditional\tlnL_null\tlnL_alt\n')
    if obs: o.write(f'OBSERVED (real alignment)\t{obs[0]:.4f}\t{obs[1]:.5f}\t\t\n')
    for r in res: o.write(f'{r[0]}\t{r[1]:.4f}\t{r[2]:.5f}\t{r[3]:.3f}\t{r[4]:.3f}\n')
if res:
    Kb=sorted(r[2] for r in res); Lb=sorted(r[1] for r in res)
    q=lambda v,pq: v[min(len(v)-1,int(round(pq*(len(v)-1))))]
    bs={'n':len(res),'K_median':st.median(Kb),'K_2.5%':q(Kb,0.025),'K_97.5%':q(Kb,0.975),'K_min':Kb[0],'K_max':Kb[-1],'LRT_median':st.median(Lb),'LRT_95%':q(Lb,0.95),'LRT_max':Lb[-1],'frac_LRT_gt_3.84':sum(1 for x in Lb if x>3.84)/len(Lb)}
    if obs: bs.update({'obs_LRT':obs[0],'obs_K':obs[1],'p_LRT_ge_obs':(1+sum(1 for x in Lb if x>=obs[0]))/(len(Lb)+1),'p_absdev_K_ge_obs':(1+sum(1 for x in Kb if abs(x-1)>=abs(obs[1]-1)))/(len(Kb)+1)})
    json.dump(bs,open(T+'bootstrap_summary.json','w'),indent=1); print(bs)
# full refits on simulated data
for p in sorted(glob.glob(D+'boot/full/rep*.json')):
    if os.path.getsize(p)>0: s=summ(p); print('FULL',os.path.basename(p),round(s['K'],5),round(s['LRT'],3),s['alt_test'])
