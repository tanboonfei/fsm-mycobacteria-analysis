#!/usr/bin/env python3
"""Orthologue groups called decayed in >=50% of the FSM SPIRE MAGs (>=95%) in
which they were assessed, with the modal non-hypothetical PGAP product name."""
import os, pandas as pd
OUT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = '/path/to/project'
L = pd.read_csv(OUT + '/results/per_locus_calls.tsv.gz', sep='\t')
g = pd.read_csv(OUT + '/results/per_genome_decay_screen_grouped.tsv', sep='\t')
L = L.merge(g[['genome', 'group']], on='genome'); f = L[L.group == 'FSM SPIRE MAG (>=95%)']
per = f.groupby('og').agg(n=('any_decay', 'size'), dec=('any_decay', 'sum'), fs=('disrupt', 'mean'),
                          short=('short80', 'mean'), split=('split', 'mean'), qcov=('qcov', 'median'), ident=('ident', 'median'))
per['frac'] = per.dec / per.n; per = per[per.frac >= 0.5].sort_values('frac', ascending=False)
og = pd.read_csv(OUT + '/results/orthogroups.tsv', sep='\t'); names = {}
for gg in og.genome.unique():
    for line in open(f'{ROOT}/refprot/ncbi_dataset/data/{gg}/protein.faa'):
        if line[0] == '>': names[(gg, line[1:].split()[0])] = line.split(' ', 1)[1].split('[')[0].strip()
def nm(o):
    c = [names.get((r.genome, r.protein)) for r in og[og.OG == o].itertuples()]
    c = [x for x in c if x and 'hypothetical' not in x]
    return pd.Series(c).mode()[0] if c else 'hypothetical protein'
per['product'] = [nm(o) for o in per.index]
per.to_csv(OUT + '/results/fsm_lineage_shared_decay_ogs.tsv', sep='\t', float_format='%.3g')
print(per.to_string())
