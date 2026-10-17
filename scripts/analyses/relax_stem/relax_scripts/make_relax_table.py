"""Build the RELAX convergence sub-table (S9r) from the collected run tables."""
import csv
from scipy.stats import chi2
T = '/path/to/project/analyses/relax_stem/tables/'
runs = list(csv.DictReader(open(T + 'relax_all_runs.tsv'), delimiter='\t'))
jobs = list(csv.DictReader(open(T + 'relax_polish_profile_jobs.tsv'), delimiter='\t'))
f = lambda x: float(x) if x not in ('', 'None', None) else None
rows = []
label = {'published (RELAX 3.1.1, 1 start)': 'Published run (RELAX 3.1.1, 1 starting point)',
         'extended rerun, identical command': 'Rerun of the published command (extended analysis)',
         'R_FSM_s1_sp5': 'RELAX, --starting-points 5, seed 201', 'R_FSM_s2_sp5': 'RELAX, --starting-points 5, seed 202', 'R_FSM_s3_sp5': 'RELAX, --starting-points 5, seed 203', 'R_FSM_s4_sp1': 'RELAX, published command (1 starting point), seed 204',
         'R_LEP_s1_sp1': 'RELAX rerun of the published command, seed 101', 'R_LEP_s1_sp5': 'RELAX, --starting-points 5, seed 301',
         'R_LEP_s2_sp5': 'RELAX, --starting-points 5, seed 302',
         'W_FSM_cond': 'Warm start from conditional optimum: null, then alternative (HyPhy default optimiser)',
         'W_LEP_cond': 'Warm start (FSM conditional-optimum omegas, relabelled tree): null, then alternative',
         'W_FSM_xLEPbest': 'Warm start from best leprosy null (relabelled): null, then alternative',
         'W_LEP_xFSMbest': 'Warm start from best FSM null (relabelled): null, then alternative'}
for r in runs:
    if r['run'] not in label: continue
    n, a = f(r['lnL_null']), f(r['lnL_alt'])
    rows.append(dict(test_set='FSM terminals+stem' if r['test'] == 'FSM' else 'Leprosy clade terminals+stem', procedure=label[r['run']],
        lnL_null='%.2f' % n if n else '', lnL_alt='%.2f' % a if a else '', LRT='%.2f' % (2 * (a - n)) if n and a else '',
        p_chi2_1df=('%.2g' % chi2.sf(max(2 * (a - n), 0), 1)) if n and a else '', K='%.4f' % f(r['K']) if f(r['K']) else '',
        ref_omega_alt=r['ref_alt'], test_omega_alt=r['test_alt'], runtime_min='%.0f' % (f(r['runtime_s']) / 60) if f(r['runtime_s']) else '',
        note=('alternative below null: optimiser failure' if n and a and a < n else '')))
jl = {'L_FSM_sp5null_cw': 'Coordinate-wise re-optimisation of the seed-201 null', 'L_FSM_sp5alt_cw': 'Coordinate-wise re-optimisation of the seed-201 alternative',
      'X_FSM_alt_from_sp5null_cw': 'Alternative re-optimised from the seed-201 null (K start 1.04, coordinate-wise)',
      'X_FSM_null_from_sp5alt_cw': 'Null re-optimised from the seed-201 alternative (coordinate-wise)',
      'L_FSM_null': 'Null re-optimised from the K = 1.04 profile solution (second basin)', 'P_FSM_K1.00_fromK1.04': 'Null from the K = 1.04 profile solution (default optimiser, 1 pass)', 'L_FSM_alt': 'Alternative re-optimised (second basin)',
      'L_LEP_rnull_cw': 'Coordinate-wise re-optimisation of the seed-101 null', 'L_LEP_ralt_cw': 'Coordinate-wise re-optimisation of the seed-101 alternative',
      'A_LEP_free_fromK0.76': 'Alternative re-optimised from the K = 0.76 profile solution'}
for j in jobs:
    if j['job'] not in jl: continue
    rows.append(dict(test_set='FSM terminals+stem' if j['test'] == 'FSM' else 'Leprosy clade terminals+stem', procedure=jl[j['job']],
        lnL_null='%.2f' % f(j['best_lnL']) if j['model'].startswith('null') else '', lnL_alt='%.2f' % f(j['best_lnL']) if j['model'].startswith('alt') else '',
        LRT='', p_chi2_1df='', K='%.4f' % f(j['K']), ref_omega_alt='', test_omega_alt='', runtime_min='%.0f' % (f(j['runtime_s']) / 60) if f(j['runtime_s']) else '',
        note='successive-pass lnL changes: ' + (j['pass_deltas'] or 'none') + ('; ' + j['done'] if j['done'] != 'finished' else '')))
with open(T + 'S9r_relax_convergence.tsv', 'w', newline='') as o:
    w = csv.DictWriter(o, fieldnames=list(rows[0].keys()), delimiter='\t'); w.writeheader(); w.writerows(rows)
print(len(rows), 'rows')
