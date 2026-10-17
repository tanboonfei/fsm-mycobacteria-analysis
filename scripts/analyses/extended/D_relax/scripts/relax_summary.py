# Summarise a HyPhy RELAX json: K, LRT, p, branch counts, summed MG94 test length, branch-set omegas, alt/null distributions
import json, sys, hashlib
def summ(p):
    d=json.load(open(p)); t=d['test results']; tested=d['tested']['0']; ba=d['branch attributes']['0']
    mg='MG94xREV with separate rates for branch sets'
    lens={g:sum(ba[b][mg] for b in tested if tested[b]==g and b in ba) for g in ('Test','Reference','Unclassified')}
    n={g:sum(1 for b in tested if tested[b]==g) for g in ('Test','Reference','Unclassified')}
    rd=d['fits'][mg]['Rate Distributions']
    om={g:rd[k][0][0] for k in rd for g in ('Test','Reference','Unclassified') if '*'+g+'*' in k}
    alt=d['fits']['RELAX alternative']['Rate Distributions']; nul=d['fits']['RELAX null']['Rate Distributions']
    fmt=lambda x:';'.join(f"{v['omega']:.4g}@{100*v['proportion']:.2f}%" for k,v in sorted(x.items()))
    return dict(K=t['relaxation or intensification parameter'], LRT=t['LRT'], p=t['p-value'], n_test=n['Test'], n_ref=n['Reference'],
                n_uncl=n['Unclassified'], len_test=lens['Test'], len_ref=lens['Reference'], omega_test=om.get('Test'), omega_ref=om.get('Reference'),
                alt_test=fmt(alt['Test']), alt_ref=fmt(alt['Reference']), null=fmt(nul['Test']),
                lnL_alt=d['fits']['RELAX alternative']['Log Likelihood'], lnL_null=d['fits']['RELAX null']['Log Likelihood'],
                runtime_s=d['timers']['Overall']['timer'], sha256=hashlib.sha256(open(p,'rb').read()).hexdigest())
if __name__=='__main__':
    for p in sys.argv[1:]:
        s=summ(p); print(p.split('/')[-1], {k:(round(v,5) if isinstance(v,float) else v) for k,v in s.items() if k!='sha256'})
