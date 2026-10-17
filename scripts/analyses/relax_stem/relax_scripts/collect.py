"""Collect all relax_stem RELAX runs (warm-start HBL fits and seeded RELAX runs) into one TSV + print."""
import json, glob, os, re, sys
from math import log
from scipy.stats import chi2
R = '/path/to/project/analyses/relax_stem/relax/runs'
A8 = '/path/to/project/analyses/extended/D_relax'
def dist(om, aux0, aux1):
    w1 = aux0; w2 = (1 - aux0) * aux1; w3 = 1 - w1 - w2
    return sorted(zip(om, (w1, w2, w3)))
def fmt(d): return '; '.join('%.4g (%.1f%%)' % (o, 100 * w) for o, w in d)
def wall(prefix):
    f = prefix + '.time'
    if os.path.exists(f):
        m = re.search(r'Elapsed \(wall clock\) time.*: (.*)', open(f).read())
        if m:
            p = [float(x) for x in m.group(1).split(':')]
            return int(p[-1] + 60 * p[-2] + (3600 * p[-3] if len(p) > 2 else 0))
    return None
rows = []
def add(test, run, kind, start, lnN, lnA, K, refN, refA, testA, unc, secs, note=''):
    lrt = 2 * (lnA - lnN) if lnN is not None and lnA is not None else None
    p = chi2.sf(max(lrt, 0), 1) if lrt is not None else None
    rows.append(dict(test=test, run=run, kind=kind, start=start, lnL_null=lnN, lnL_alt=lnA, LRT=lrt, p=p, K=K,
                     ref_null=refN, ref_alt=refA, test_alt=testA, unclassified=unc, runtime_s=secs, note=note))
# published + extended rerun
for test, f, run in [('FSM', A8 + '/inputs/common_reference_FSM_tips_stem.json', 'published (RELAX 3.1.1, 1 start)'),
                     ('LEP', A8 + '/inputs/common_reference_LEP_tips_stem.json', 'published (RELAX 3.1.1, 1 start)'),
                     ('FSM', A8 + '/boot/real/FSM_tips_stem.json', 'extended rerun, identical command')]:
    j = json.load(open(f)); fa = j['fits']['RELAX alternative']; fn = j['fits']['RELAX null']
    g = lambda d: fmt(sorted((v['omega'], v['proportion']) for v in d.values()))
    add(test, run, 'RELAX', 'default', fn['Log Likelihood'], fa['Log Likelihood'], j['test results']['relaxation or intensification parameter'],
        g(fn['Rate Distributions']['Reference']), g(fa['Rate Distributions']['Reference']), g(fa['Rate Distributions']['Test']), '', j['timers']['Overall']['timer'])
add('FSM', 'extended conditional refit (branch lengths fixed)', 'conditional', 'extended rerun null', -694918.7973694124, None, None,
    fmt(dist((0.03165687374094876, 0.9898442458781159, 17.57427538892248), 0.9185622357422768, 0.748553643135527)), '', '', '', None, 'null only; same objective')
# seeded RELAX runs
for f in sorted(glob.glob(R + '/R_*.json')):
    try:
        j = json.load(open(f)); j['test results']['LRT']
    except Exception:
        continue
    name = os.path.basename(f)[:-5]
    test = name.split('_')[1]; fa = j['fits']['RELAX alternative']; fn = j['fits']['RELAX null']
    g = lambda d: fmt(sorted((v['omega'], v['proportion']) for v in d.values()))
    seed, sp = re.match(r'R_\w+?_s(\d+)_sp(\d+)', name).groups()
    add(test, name, 'RELAX', 'seed %s, starting-points %s' % (seed, sp), fn['Log Likelihood'], fa['Log Likelihood'],
        j['test results']['relaxation or intensification parameter'], g(fn['Rate Distributions']['Reference']), g(fa['Rate Distributions']['Reference']),
        g(fa['Rate Distributions']['Test']), '', j['timers']['Overall']['timer'])
# warm-start fits
for f in sorted(glob.glob(R + '/W_*.log')):
    name = os.path.basename(f)[:-4]; test = name.split('_')[1]; res = {}
    for l in open(f):
        if l.startswith('RESULT'):
            t = l.rstrip('\n').split('\t'); res[t[1]] = dict(zip(t[2::2], t[3::2]))
    def parse(tag):
        if tag not in res: return None
        d = res[tag]; ref = [float(x) for x in d['ref'].split(',')]; tst = [float(x) for x in d['test'].split(',')]; unc = [float(x) for x in d['unc'].split(',')]
        return dict(lnL=float(d['lnL']), K=float(d['K']), secs=float(d['secs']), ref=fmt(dist(ref[:3], ref[3], ref[4])),
                    test=fmt(dist(tst, ref[3], ref[4])), unc=fmt(dist(unc[:3], unc[3], unc[4])))
    n, a, ap, npol = parse('null_fit'), parse('alt_fit'), parse('alt_repolish'), parse('null_repolish')
    first = 'null' if 'cond' in name or 'pubnull' in name else 'alt'
    if npol and n and npol['lnL'] > n['lnL']: n = npol
    if ap and a and ap['lnL'] > a['lnL']: a = ap
    secs = sum(x['secs'] for x in (parse('null_fit'), parse('alt_fit'), ap, npol) if x)
    add(test, name, 'warm start (all parameters free)', name.split('_', 2)[2] + ' (%s first)' % first,
        n and n['lnL'], a and a['lnL'], a and a['K'], n and n['ref'], a and a['ref'], a and a['test'], (a or n or {}).get('unc', ''), wall(R + '/' + name) or secs,
        note='' if (n and a) else 'incomplete')
import csv
out = '/path/to/project/analyses/relax_stem/tables/relax_all_runs.tsv'
with open(out, 'w', newline='') as o:
    w = csv.DictWriter(o, fieldnames=list(rows[0].keys()), delimiter='\t'); w.writeheader(); w.writerows(rows)
for r in rows:
    print('%-4s %-48s null %s alt %s LRT %s K %s t %s %s' % (r['test'], r['run'][:48], r['lnL_null'] and '%.2f' % r['lnL_null'], r['lnL_alt'] and '%.2f' % r['lnL_alt'],
          r['LRT'] is not None and '%.2f' % r['LRT'], r['K'] and '%.4f' % r['K'], r['runtime_s'], r['note']))

# ---- single-model jobs run with profile.bf / profile_cw.bf (P_, A_, L_, X_, Q_) ----
prof = []
for f in sorted(glob.glob(R + '/[PALXQ]_*.log')):
    name = os.path.basename(f)[:-4]; test = name.split('_')[1]
    lines = open(f).read().split('\n')
    st = [l for l in lines if l.startswith('START')]
    passes = [l.split('\t') for l in lines if l.startswith('PASS')]
    resl = [l.split('\t') for l in lines if l.startswith('RESULT')]
    if not st: continue
    k0 = float(st[0].split()[-1]); ln0 = float(st[0].split()[2])
    free = name.startswith('A_') or 'alt' in name.split('_from')[0]
    # best lnL seen: start value (a valid evaluated point) or any pass
    vals = [(ln0, k0, 'start')] + [(float(p[3]), float(p[7]), 'pass%s' % p[1]) for p in passes]
    if resl and not passes:
        r = resl[0]; vals.append((float(r[3]), float(r[7]), 'pass1'))
    best = max(vals)
    model = 'alt (K free)' if free else ('null (K=1)' if abs(k0 - 1) < 1e-9 else 'profile K=%g' % k0)
    cmd = [l for l in open(R + '/jobs.tsv') if '\t' + name + '\t' in l]
    method = 'coordinate-wise' if cmd and 'profile_cw' in cmd[0] else 'HyPhy default'
    deltas = ' / '.join('%+.2f' % float(p[5]) for p in passes)
    prof.append(dict(test=test, job=name, model=model, method=method, start_lnL=ln0, best_lnL=best[0], K=best[1], best_at=best[2],
                     pass_deltas=deltas, done=('finished' if resl else ('stopped after pass 1 (slow later pass)' if any(l.startswith('END\t' + name + '\t') for l in open(R + '/sched_status.tsv')) else 'RUNNING')), runtime_s=wall(R + '/' + name)))
out2 = '/path/to/project/analyses/relax_stem/tables/relax_polish_profile_jobs.tsv'
with open(out2, 'w', newline='') as o:
    w = csv.DictWriter(o, fieldnames=list(prof[0].keys()), delimiter='\t'); w.writeheader(); w.writerows(prof)
print()
for p in prof:
    print('%-4s %-28s %-13s %-15s start %.2f best %.2f K %.4f (%s) deltas %s %s' % (p['test'], p['job'], p['model'], p['method'], p['start_lnL'], p['best_lnL'], p['K'], p['best_at'], p['pass_deltas'], p['done']))
