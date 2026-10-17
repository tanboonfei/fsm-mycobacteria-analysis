#!/usr/bin/env python3
"""Build the genome manifest and the shredded fragmentation controls.
Shredding re-implements publication/analysis_fragmentation_sim/shred_genomes.py
(the S3b-S3e simulation; its FASTA outputs were not retained): exponential
fragment lengths, mean 16.7 kb (N50 ~28 kb), fragments <500 bp dropped, seed
20260910; FI additionally deletes random fragments until 95.6% remains."""
import os, gzip, random, pandas as pd
ROOT = '/path/to/project'
OUT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FNA = OUT + '/work/fna'; os.makedirs(FNA, exist_ok=True)
MEAN, KEEP = 16700, 0.956
S3 = """GCF_010731895.1 GCF_010726245.1 GCA_002007745.1 GCF_010731855.1 GCF_015710575.1
GCF_010725885.1 GCF_900637205.1 GCF_010729895.1 GCF_010728925.1 GCF_000230895.2 GCF_010726085.1
GCA_900637135.1 GCF_010731795.1 GCF_001583415.1 GCF_900108565.1 GCF_010725725.1 GCF_019645835.1
GCF_900603025.1 GCF_003855255.1 GCF_022179545.1""".split()
LEP = ['GCF_000195855.1', 'GCF_000975265.2', 'GCF_003408705.1']
paths = dict(l.rstrip('\n').split('\t') for l in open(ROOT + '/Publication/scripts/analyses/protist/A/genome_paths.tsv'))
s1 = pd.read_csv(ROOT + '/Publication/scripts/inputs/tsv/S1_genome_metadata.tsv', sep='\t')
s1['cl'] = s1.clade.str.split(' ').str[0]
sis = set(open(ROOT + '/analyses/sister_tips.txt').read().split())
asm = pd.read_csv(ROOT + '/analysis_gaps/assembly_stats.tsv', sep='\t').set_index('genome_id_or_accession')

def read(p):
    op = gzip.open if p.endswith('.gz') else open
    seqs = []; cur = []
    for line in op(p, 'rt'):
        if line.startswith('>'):
            if cur: seqs.append(''.join(cur)); cur = []
        else: cur.append(line.strip())
    if cur: seqs.append(''.join(cur))
    return seqs
def write(p, fr):
    with open(p, 'w') as f:
        for k, s in enumerate(fr, 1):
            f.write('>ctg%d\n' % k)
            for j in range(0, len(s), 80): f.write(s[j:j + 80] + '\n')
def shred(seqs, rng):
    fr = []
    for s in seqs:
        i = 0
        while i < len(s):
            L = max(500, int(rng.expovariate(1.0 / MEAN))); fr.append(s[i:i + L]); i += L
    return [f for f in fr if len(f) >= 500]
def n50(lens):
    lens = sorted(lens, reverse=True); t = sum(lens) / 2; c = 0
    for L in lens:
        c += L
        if c >= t: return L

rows = []
rng = random.Random(20260910)
for acc in S3 + LEP:
    tip = s1.loc[s1.genome_id_or_accession == acc].iloc[0]
    seqs = read(paths[tip.tip]); tot = sum(map(len, seqs))
    grp = 'closed_RGM' if acc in S3 else 'leprosy'
    rows.append(dict(genome=acc + '__C', source=acc, panel=grp + '_intact', clade=tip.cl, species=tip.species,
                     fasta=paths[tip.tip], prodigal='', contigs=len(seqs), N50=n50(map(len, seqs)), size_bp=tot,
                     completeness=100.0, sister=tip.tip in sis))
    fr = shred(seqs, rng); p = f'{FNA}/{acc}__F.fna'; write(p, fr)
    rows.append(dict(genome=acc + '__F', source=acc, panel=grp + '_shredded', clade=tip.cl, species=tip.species,
                     fasta=p, prodigal='', contigs=len(fr), N50=n50(map(len, fr)), size_bp=sum(map(len, fr)),
                     completeness=100.0, sister=tip.tip in sis))
    order = list(range(len(fr))); rng.shuffle(order); target = tot * KEEP; kept = []; a = 0
    for i in order:
        if a + len(fr[i]) <= target: kept.append(i); a += len(fr[i])
    if grp == 'closed_RGM':  # FI only for the RGM controls (the leprosy FI adds nothing to a positive control)
        kept.sort(); ff = [fr[i] for i in kept]; p = f'{FNA}/{acc}__FI.fna'; write(p, ff)
        rows.append(dict(genome=acc + '__FI', source=acc, panel=grp + '_shredded_incomplete', clade=tip.cl,
                         species=tip.species, fasta=p, prodigal='', contigs=len(ff), N50=n50(map(len, ff)),
                         size_bp=sum(map(len, ff)), completeness=95.6, sister=tip.tip in sis))
# MAGs: all SPIRE MAGs, plus the GTDB-deposited FSM MAGs
for _, t in s1.iterrows():
    gid = t.genome_id_or_accession
    if t.origin == 'SPIRE MAG':
        a = asm.loc[gid]
        rows.append(dict(genome=gid, source=gid, panel='SPIRE_MAG', clade=t.cl, species=t.species,
                         fasta=paths.get(t.tip, paths.get(gid)), prodigal=f'{ROOT}/proteins/{gid}.faa', contigs=int(a.n_contigs),
                         N50=int(a.N50), size_bp=int(a.asm_len), completeness=t.checkm2_this_study_completeness_pct,
                         sister=t.tip in sis))
    elif t.cl == 'FSM':
        seqs = read(paths[t.tip])
        rows.append(dict(genome=gid, source=gid, panel='GTDB_FSM_MAG', clade='FSM', species=t.species,
                         fasta=paths[t.tip], prodigal='', contigs=len(seqs), N50=n50(map(len, seqs)),
                         size_bp=sum(map(len, seqs)), completeness=t.completeness_used_pct, sister=False))
m = pd.DataFrame(rows)
m['tip'] = m.source.map(dict(zip(s1.genome_id_or_accession, s1.tip)))
m.to_csv(OUT + '/results/manifest.tsv', sep='\t', index=False)
print(m.groupby(['panel', 'clade']).agg(n=('genome', 'size'), N50=('N50', 'median'), contigs=('contigs', 'median')))
