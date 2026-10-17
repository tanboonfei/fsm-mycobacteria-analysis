"""Table S3f (old S4f) completeness-floor and extrapolation, old vs corrected labels.
Reimplemented; the old-label run must reproduce the published table."""
import numpy as np, pandas as pd, statsmodels.api as sm
from scipy import stats
G='/path/to/project/analyses/extended/G_culture/'
s=pd.read_csv('/path/to/project/scripts/inputs/tsv/S1_genome_metadata.tsv',sep='\t')
c=pd.read_csv(G+'new_completeness.tsv',sep='\t'); s=s.merge(c[['tip','new_status','completeness_used_new']],on='tip')
s['c3']=s.clade.str.split().str[0]
rows=[]
for lab in ['old','corrected']:
    if lab=='old': mag=s.culture_status=='MAG / uncultured'; comp='completeness_used_pct'
    else: mag=~s.new_status.str.startswith('isolate'); comp='completeness_used_new'
    F=s[mag&(s.c3=='FSM')]; R=s[mag&(s.c3=='RGM')]
    for fl in [None,95,98,99,99.5]:
        f=F if fl is None else F[F[comp]>=fl]; r=R if fl is None else R[R[comp]>=fl]
        rows.append(dict(labels=lab,floor='none (all MAGs)' if fl is None else f'>= {fl}%',FSM_MAGs=len(f),FSM_median=round(f.genome_size_Mb.median(),2),
            RGM_MAGs=len(r),RGM_median=round(r.genome_size_Mb.median(),2),gap=round(r.genome_size_Mb.median()-f.genome_size_Mb.median(),2),
            p=f'{stats.mannwhitneyu(f.genome_size_Mb,r.genome_size_Mb).pvalue:.2g}'))
    for nm,x in [('FSM',F),('RGM',R)]:
        X=sm.add_constant(x[comp].values); m=sm.OLS(x.genome_size_Mb.values,X).fit()
        pr=m.get_prediction(np.array([[1,100.0]])).summary_frame(alpha=0.05)
        rows.append(dict(labels=lab,floor=f'extrapolation {nm}',FSM_MAGs=len(x),FSM_median=round(m.params[1],4),RGM_MAGs=f'{m.pvalues[1]:.2g}',
            RGM_median=round(pr['mean'].iloc[0],2),gap=f"{pr['mean_ci_lower'].iloc[0]:.2f} - {pr['mean_ci_upper'].iloc[0]:.2f}",p=''))
o=pd.DataFrame(rows); o.to_csv(G+'S3f_floor_old_vs_corrected.tsv',sep='\t',index=False)
print('(extrapolation rows: FSM_MAGs=n, FSM_median=slope, RGM_MAGs=slope p, RGM_median=predicted Mb at 100%, gap=95% CI)')
print(o.to_string())
# composition of the RGM >=99.5% set under old labels
R=s[(s.culture_status=='MAG / uncultured')&(s.c3=='RGM')&(s.completeness_used_pct>=99.5)]
print('\nold RGM MAGs >=99.5%:',len(R),'of which now isolates:',int(R.new_status.str.startswith('isolate').sum()))
