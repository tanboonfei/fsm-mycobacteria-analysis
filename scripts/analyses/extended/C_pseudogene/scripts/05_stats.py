#!/usr/bin/env python3
"""Panel definition (fragmentation/completeness matching), group summaries,
effect sizes and the Table S6 sub-table."""
import os, numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
from scipy.optimize import linear_sum_assignment
OUT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
g = pd.read_csv(OUT + '/results/per_genome_decay_screen.tsv', sep='\t')
g['logN50'] = np.log10(g.N50); g['logcpm'] = np.log10(g.contigs_per_Mb)
mag = g[(g.panel == 'SPIRE_MAG') & (g.completeness >= 95)]
fsm = mag[mag.clade == 'FSM']; rgm = mag[mag.clade == 'RGM']
# optimal 1:1 matching on standardised log N50 and log contigs per Mb
Z = pd.concat([fsm, rgm])[['logN50', 'logcpm']]; mu, sd = Z.mean(), Z.std()
A = ((fsm[['logN50', 'logcpm']] - mu) / sd).values; B = ((rgm[['logN50', 'logcpm']] - mu) / sd).values
D = np.sqrt(((A[:, None, :] - B[None, :, :]) ** 2).sum(-1))
i, j = linear_sum_assignment(D)
rgm_m = rgm.iloc[j]
pairs = pd.DataFrame({'FSM': fsm.genome.values[i], 'RGM_match': rgm_m.genome.values, 'distance': D[i, j]})
pairs.to_csv(OUT + '/results/matched_pairs.tsv', sep='\t', index=False, float_format='%.3f')
g['group'] = ''
g.loc[g.genome.isin(fsm.genome), 'group'] = 'FSM SPIRE MAG (>=95%)'
g.loc[g.genome.isin(rgm.genome), 'group'] = 'RGM SPIRE MAG (>=95%), unmatched pool'
g.loc[g.genome.isin(rgm_m.genome), 'group'] = 'RGM SPIRE MAG (>=95%), matched'
lab = {'closed_RGM_intact': 'Closed RGM, intact', 'closed_RGM_shredded': 'Closed RGM, shredded',
       'closed_RGM_shredded_incomplete': 'Closed RGM, shredded + 4.4% deleted',
       'leprosy_intact': 'Leprosy clade, intact', 'leprosy_shredded': 'Leprosy clade, shredded'}
for k, v in lab.items(): g.loc[g.panel == k, 'group'] = v
g['sister_clade'] = g.sister
g.to_csv(OUT + '/results/per_genome_decay_screen_grouped.tsv', sep='\t', index=False, float_format='%.4g')

MET = [('pct_fs_or_stop', 'loci with frameshift or internal stop (%)'), ('pct_frameshift', 'frameshift (%)'),
       ('pct_internal_stop', 'internal stop (%)'), ('pct_orf_lt80', 'ORF <80% of reference length (%)'),
       ('n_split', 'core genes split into >=2 ORFs (n)'), ('pct_split', 'core genes split (%)'),
       ('pct_any_decay', 'any decay signal (%)'), ('pct_hit', 'core groups with a hit (%)'),
       ('hits_at_contig_end', 'hits within 200 bp of a contig end (n)'), ('assessed', 'loci assessed (n)'),
       ('median_identity', 'median identity to best reference'),
       ('spacer_median_bp', 'intergenic spacer median (bp)'), ('pct_spacer_gt300', 'spacers >300 bp (%)'),
       ('N50', 'N50 (bp)'), ('contigs_per_Mb', 'contigs per Mb'), ('completeness', 'completeness (%)')]
ORDER = ['FSM SPIRE MAG (>=95%)', 'RGM SPIRE MAG (>=95%), matched', 'Closed RGM, intact', 'Closed RGM, shredded',
         'Closed RGM, shredded + 4.4% deleted', 'Leprosy clade, intact', 'Leprosy clade, shredded']
groups = {k: g[g.group == k] for k in ORDER}
groups['RGM SPIRE MAG (>=95%), all'] = g[g.genome.isin(rgm.genome)]
groups['RGM sister-clade SPIRE MAG (>=95%)'] = g[g.genome.isin(rgm.genome) & g.sister]
groups['FSM GTDB-deposited MAG (>=95%)'] = g[(g.panel == 'GTDB_FSM_MAG') & (g.completeness >= 95)]
rows = []
for k, d in groups.items():
    r = {'group': k, 'n': len(d)}
    for c, _ in MET:
        r[c + '_median'] = d[c].median(); r[c + '_range'] = f'{d[c].min():.3g}-{d[c].max():.3g}'
    rows.append(r)
S = pd.DataFrame(rows); S.to_csv(OUT + '/results/group_summary.tsv', sep='\t', index=False, float_format='%.4g')

def cliff(a, b):
    a = np.asarray(a); b = np.asarray(b)
    return ((a[:, None] > b[None, :]).sum() - (a[:, None] < b[None, :]).sum()) / (len(a) * len(b))
def boot_md(a, b, n=5000, seed=1):
    rng = np.random.default_rng(seed); a = np.asarray(a); b = np.asarray(b)
    v = [np.median(rng.choice(a, len(a))) - np.median(rng.choice(b, len(b))) for _ in range(n)]
    return np.percentile(v, [2.5, 97.5])
comp = [('FSM SPIRE MAG (>=95%)', 'RGM SPIRE MAG (>=95%), matched'),
        ('FSM SPIRE MAG (>=95%)', 'RGM SPIRE MAG (>=95%), all'),
        ('FSM SPIRE MAG (>=95%)', 'RGM sister-clade SPIRE MAG (>=95%)'),
        ('FSM SPIRE MAG (>=95%)', 'Closed RGM, shredded + 4.4% deleted'),
        ('RGM SPIRE MAG (>=95%), matched', 'Closed RGM, shredded + 4.4% deleted'),
        ('Closed RGM, shredded', 'Closed RGM, intact'),
        ('FSM SPIRE MAG (>=95%)', 'Leprosy clade, intact')]
out = []
for a, b in comp:
    for c, desc in MET[:13]:
        x, y = groups[a][c].dropna(), groups[b][c].dropna()
        p = mannwhitneyu(x, y).pvalue if len(y) > 2 and len(x) > 2 else np.nan
        lo, hi = boot_md(x, y)
        out.append(dict(group_a=a, group_b=b, metric=desc, n_a=len(x), median_a=x.median(), n_b=len(y),
                        median_b=y.median(), difference=x.median() - y.median(), diff_CI95=f'{lo:.3g} to {hi:.3g}',
                        cliffs_delta=cliff(x, y), MannWhitney_p=p))
E = pd.DataFrame(out); E.to_csv(OUT + '/results/effect_sizes.tsv', sep='\t', index=False, float_format='%.4g')
# matching balance
bal = []
for c in ['N50', 'contigs', 'contigs_per_Mb', 'completeness', 'size_bp']:
    bal.append(dict(variable=c, FSM_median=fsm[c].median(), RGM_matched_median=rgm_m[c].median(),
                    RGM_all_median=rgm[c].median(), p_FSM_vs_matched=mannwhitneyu(fsm[c], rgm_m[c]).pvalue))
pd.DataFrame(bal).to_csv(OUT + '/results/matching_balance.tsv', sep='\t', index=False, float_format='%.4g')
# Table S6 sub-table: one row per panel
keep = [('pct_hit', 'Core groups with a hit (%)'), ('assessed', 'Loci assessed (n, >200 bp from contig ends)'),
        ('median_identity', 'Median identity to best reference'),
        ('pct_frameshift', 'Frameshift (%)'), ('pct_internal_stop', 'Internal stop (%)'),
        ('pct_fs_or_stop', 'Frameshift or internal stop (%)'), ('pct_orf_lt80', 'Prodigal ORF <80% of reference (%)'),
        ('n_split', 'Core genes split into >=2 ORFs (n)'), ('pct_any_decay', 'Any decay signal (%)'),
        ('spacer_median_bp', 'Intergenic spacer median (bp)'), ('pct_spacer_gt300', 'Spacers >300 bp (%)')]
T = []
for k in ORDER + ['RGM sister-clade SPIRE MAG (>=95%)', 'FSM GTDB-deposited MAG (>=95%)']:
    d = groups[k]
    r = {'Panel': k, 'Genomes': len(d), 'N50 median (kb)': round(d.N50.median() / 1e3, 1),
         'Contigs per Mb median': round(d.contigs_per_Mb.median(), 1)}
    for c, nm in keep:
        if len(d) <= 3:
            r[nm] = '; '.join(f'{v:.3g}' for v in d[c])
        else:
            r[nm] = f'{d[c].median():.3g} ({d[c].quantile(.25):.3g}-{d[c].quantile(.75):.3g})'
    T.append(r)
pd.DataFrame(T).to_csv(OUT + '/results/Table_S6i_orthologue_decay_screen.tsv', sep='\t', index=False)
print(pd.DataFrame(T).T.to_string())
print(pd.DataFrame(bal).to_string())
print(E[E.metric.isin(['loci with frameshift or internal stop (%)', 'ORF <80% of reference length (%)',
     'core genes split into >=2 ORFs (n)', 'any decay signal (%)', 'intergenic spacer median (bp)'])].to_string())
