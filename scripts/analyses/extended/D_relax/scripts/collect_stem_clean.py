import json, sys, os, hashlib
sys.path.insert(0, os.path.dirname(__file__)); from relax_summary import summ
D='/path/to/project/analyses/extended/D_relax/'; T=D+'tables/'; os.makedirs(T,exist_ok=True)
def prof(p):
    pts=[];mle=None;ll=None
    for l in open(p):
        if l.startswith('PROFILE'):
            a=l.rstrip().split('\t')
            if a[1]=='MLE': mle=float(a[2])
            else: pts.append((float(a[1]),float(a[3])))
    lo=[x for x in pts if x[0]<mle]; hi=[x for x in pts if x[0]>mle]
    def cross(s):
        for (x1,y1),(x2,y2) in zip(s,s[1:]):
            if (y1-3.8415)*(y2-3.8415)<=0: return x1+(3.8415-y1)*(x2-x1)/(y2-y1)
    return mle,cross(lo),cross(hi)
def hyci(p):
    for l in open(p):
        if 'rate ratio for *test*' in l and '[' in l:
            a=l.split('=')[1].replace('[',' ').replace(']',' ').replace('-',' ').split(); return float(a[0]),float(a[1]),float(a[2])
ab=json.load(open(D+'absrel/absrel_stems.json'))['branch attributes']['0']
hc=json.load(open(D+'absrel/absrel_highclass_sites_by_OG.json'))
rows=[]
for nm,b,tag in [('FSM stem','IntNode024','FSMstem'),('Leprosy-clade stem','IntNode017','LEPstem')]:
    x=ab[b]; rd=x['Rate Distributions']
    bu=json.load(open(D+f'bustedE/bustedE_{tag}.json')); ut=bu['fits']['Unconstrained model']['Rate Distributions']['Test']
    m,l,h=prof(D+f'ci/profile_{tag}.log'); hm,hl,hh=hyci(D+f'ci/singleomega_{tag}.log')
    mA,lA,hA=prof(D+f'ci/profile_{tag}cleanA.log')
    rows.append([nm,b,f"{x['Baseline MG94xREV']:.5f}",
      f"{m:.4f}",f"{l:.4f}",f"{h:.4f}",f"{hl:.4f}-{hh:.4f}",f"{mA:.4f}",f"{lA:.4f}",f"{hA:.4f}",
      f"{x['Baseline MG94xREV omega ratio']:.4f}",x['Rate classes'],'; '.join(f"{r[0]:.4g}@{100*r[1]:.2f}%" for r in rd),f"{x['LRT']:.2f}",f"{x['Uncorrected P-value']:.3g}",f"{x['Corrected P-value']:.3g}",
      hc[b]['n_sites'],hc[b]['n_OG'],
      '; '.join(f"{v['omega']:.4g}@{100*v['proportion']:.2f}%" for k,v in sorted(ut.items())),f"{bu['test results']['LRT']:.2f}",f"{bu['test results']['p-value']:.3g}"])
hdr=['branch','node','branch_length_MG94_baseline','omega_stem_MG94_two_ratio','omega_95CI_lower_profile','omega_95CI_upper_profile','omega_95CI_HyPhy_conditional','omega_stem_cleanA','omega_cleanA_CI_lower','omega_cleanA_CI_upper','aBSREL_baseline_omega','aBSREL_rate_classes','aBSREL_rate_distribution','aBSREL_LRT','aBSREL_p_uncorrected','aBSREL_p_Holm','aBSREL_sites_P_highclass_gt0.9','aBSREL_OGs_with_such_sites','BUSTED-E_stem_distribution(class0=error sink)','BUSTED-E_LRT','BUSTED-E_p']
with open(T+'S9n_stem_tests.tsv','w') as o:
    o.write('\t'.join(hdr)+'\n'); [o.write('\t'.join(map(str,r))+'\n') for r in rows]
# cleaning
ca=json.load(open(D+'clean/cleanA_stats.json')); cb=json.load(open(D+'clean/cleanB_stats.json'))
runs=[('none (published)','FSM',D+'inputs/common_reference_FSM_tips_stem.json','21280'),('none (published)','Leprosy clade',D+'inputs/common_reference_LEP_tips_stem.json','21280'),
      ('A: BUSTED-E error-sink masking','FSM',D+'clean/cleanA_FSM_tips_stem.json','21280'),('A: BUSTED-E error-sink masking','Leprosy clade',D+'clean/cleanA_LEP_tips_stem.json','21280'),
      ('B: segment masking + trimAl automated1','FSM',D+'clean/cleanB_FSM_tips_stem.json',str(cb['cols_kept'])),('B: segment masking + trimAl automated1','Leprosy clade',D+'clean/cleanB_LEP_tips_stem.json',str(cb['cols_kept']))]
with open(T+'S9o_alignment_cleaning.tsv','w') as o:
    o.write('cleaning\ttest_set\tcodons\tcells_masked\tK\tLRT\tp\tomega_test_MG94\tomega_reference_MG94\tRELAX_alt_test\tRELAX_alt_reference\tRELAX_null\tlnL_alt\tsha256\n')
    for c,t,p,n in runs:
        s=summ(p); m={'none (published)':'0','A: BUSTED-E error-sink masking':str(ca['n_codon_cells_masked']),'B: segment masking + trimAl automated1':str(cb['masked_residues'])+' residues + '+str(21280-cb['cols_kept'])+' columns removed'}[c]
        o.write('\t'.join([c,t,n,m,f"{s['K']:.4f}",f"{s['LRT']:.2f}",f"{s['p']:.3g}",f"{s['omega_test']:.4f}",f"{s['omega_ref']:.4f}",s['alt_test'],s['alt_ref'],s['null'],f"{s['lnL_alt']:.2f}",s['sha256']])+'\n')
print(open(T+'S9n_stem_tests.tsv').read()); print(open(T+'S9o_alignment_cleaning.tsv').read())
