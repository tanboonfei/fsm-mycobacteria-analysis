#!/usr/bin/env python3
"""Locus-level checks: (1) decay rates within bins of identity to the best
reference, FSM vs matched RGM (divergence control; Mantel-Haenszel);
(2) anatomy of short-ORF calls; (3) per-OG overlap of decayed loci across FSM."""
import os, numpy as np, pandas as pd
from statsmodels.stats.contingency_tables import StratifiedTable
OUT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
g = pd.read_csv(OUT + '/results/per_genome_decay_screen_grouped.tsv', sep='\t')
L = pd.read_csv(OUT + '/results/per_locus_calls.tsv.gz', sep='\t')
L = L.merge(g[['genome', 'group']], on='genome')
L['ibin'] = pd.cut(L.ident, [0.35, 0.6, 0.7, 0.8, 0.9, 1.01], right=False)
F, R = 'FSM SPIRE MAG (>=95%)', 'RGM SPIRE MAG (>=95%), matched'
sub = L[L.group.isin([F, R])].copy(); sub['FSM'] = sub.group == F
rows = []
for met in ['disrupt', 'short80', 'split', 'any_decay']:
    tabs = []
    for b, d in sub.groupby('ibin', observed=True):
        a = d[d.FSM][met]; r = d[~d.FSM][met]
        rows.append(dict(metric=met, identity_bin=str(b), FSM_loci=len(a), FSM_pct=100 * a.mean(),
                         RGM_loci=len(r), RGM_pct=100 * r.mean()))
        t = np.array([[a.sum(), len(a) - a.sum()], [r.sum(), len(r) - r.sum()]])
        if (t.sum(1) > 0).all(): tabs.append(t)
    st = StratifiedTable(tabs)
    rows.append(dict(metric=met, identity_bin='Mantel-Haenszel OR (FSM/RGM)', FSM_pct=st.oddsratio_pooled,
                     RGM_pct=st.test_null_odds().pvalue))
T = pd.DataFrame(rows); T.to_csv(OUT + '/results/identity_stratified.tsv', sep='\t', index=False, float_format='%.4g')
print(T.to_string())
# anatomy of short-ORF calls
def anat(d):
    s = d[d.short80]
    return pd.Series({'short80_loci': len(s), 'with_fs_or_stop_%': 100 * s.disrupt.mean(),
                      'split_%': 100 * s.split.mean(), 'no_ORF_%': 100 * (s.n_orf == 0).mean(),
                      'alignment_cov<0.8_%': 100 * (s.qcov < 0.8).mean(),
                      'intact_full_alignment_%': 100 * ((~s.disrupt) & (s.qcov >= 0.95) & (s.n_orf == 1)).mean(),
                      'median_orf_ratio': s.orf_ratio.median()})
A = L[L.group != ''].groupby('group').apply(anat)
A.to_csv(OUT + '/results/short_orf_anatomy.tsv', sep='\t', float_format='%.4g'); print(A.to_string())
# are FSM decayed loci shared (lineage-level) or genome-specific?
f = L[L.group == F]
per = f.groupby('og').agg(n=('any_decay', 'size'), dec=('any_decay', 'sum'))
print('FSM: OGs with any decay call in >=1 genome', (per.dec > 0).sum(), 'of', len(per),
      '; in >=50% of genomes where assessed', ((per.dec / per.n) >= 0.5).sum(),
      '; in >=25%', ((per.dec / per.n) >= 0.25).sum())
lep = L[L.group == 'Leprosy clade, intact']
pl = lep.groupby('og').any_decay.mean()
print('leprosy: OGs decayed in all 3:', (pl == 1).sum())
pass  # named list: scripts/08_name_shared.py
# identity-matched closed-RGM comparison: share of short80 loci among high-identity loci
# genome-specific vs lineage-shared decay: a locus call is 'shared' when its OG is
# called decayed in >=50% of the genomes of the same panel in which it was assessed
from scipy.stats import mannwhitneyu
rec = []
for grp in [F, R, 'Closed RGM, intact', 'Closed RGM, shredded + 4.4% deleted', 'Leprosy clade, intact']:
    d = L[L.group == grp]
    fr = d.groupby('og').any_decay.mean(); shared = set(fr[fr >= 0.5].index)
    for gname, x in d.groupby('genome'):
        sh = x.og.isin(shared)
        rec.append(dict(group=grp, genome=gname, n_shared_ogs_in_panel=len(shared),
                        pct_decay_shared=100 * (x.any_decay & sh).sum() / len(x),
                        pct_decay_genome_specific=100 * (x.any_decay & ~sh).sum() / len(x),
                        pct_fs_stop_genome_specific=100 * (x.disrupt & ~sh).sum() / len(x)))
S = pd.DataFrame(rec); S.to_csv(OUT + '/results/shared_vs_specific_per_genome.tsv', sep='\t', index=False, float_format='%.4g')
Q = S.groupby('group').agg(n=('genome', 'size'), shared_OGs=('n_shared_ogs_in_panel', 'first'),
    shared_median=('pct_decay_shared', 'median'), specific_median=('pct_decay_genome_specific', 'median'),
    specific_fs_stop_median=('pct_fs_stop_genome_specific', 'median'))
print(Q.to_string())
for c in ['pct_decay_genome_specific', 'pct_fs_stop_genome_specific', 'pct_decay_shared']:
    a = S[S.group == F][c]; b = S[S.group == R][c]
    print(c, 'FSM vs matched RGM: medians %.3g vs %.3g, MW p = %.3g' % (a.median(), b.median(), mannwhitneyu(a, b).pvalue))
Q.to_csv(OUT + '/results/shared_vs_specific_summary.tsv', sep='\t', float_format='%.4g')
