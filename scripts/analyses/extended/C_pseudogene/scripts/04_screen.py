#!/usr/bin/env python3
"""Orthologue-guided decay screen (Pseudofinder-like, caller-independent).
For each genome and each of the 930 single-copy core orthologue groups, the
best-scoring miniprot alignment among the 12 reference proteins of the group
is taken as the target locus. Loci within 200 bp of a contig end are excluded.
Per locus: frameshifts (miniprot fs), in-frame stop codons (st), intron-like
gaps (N/U/V CIGAR ops), reference coverage, and the Prodigal ORF structure of
the locus (longest same-strand ORF / reference length; number of ORFs)."""
import os, re, gzip, sys, collections
import numpy as np, pandas as pd
OUT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EDGE = 200          # bp from either contig end
MIN_ID = 0.35       # miniprot Identity for a locus to count as a hit
MIN_QCOV = 0.30     # aligned fraction of the reference protein for a hit
ORF_OV = 90         # nt overlap for a Prodigal ORF to be assigned to a locus
ORF_FRAC = 0.20     # ... and at least this fraction of the ORF inside the locus
m = pd.read_csv(OUT + '/results/manifest.tsv', sep='\t')
m = m[m.genome.map(lambda x: os.path.exists(f'{OUT}/work/mp/{x}.paf.gz'))]
SUFFIX = sys.argv[1] if len(sys.argv) > 1 else ''
MINLEN = int(os.environ.get('MINLEN', '100'))  # exclude groups whose median reference length is < MINLEN aa
_og = pd.read_csv(OUT + '/results/orthogroups.tsv', sep='\t').groupby('OG').length_aa.median()
OGKEEP = set(_og[_og >= MINLEN].index)

def prodigal_genes(row):
    genes = collections.defaultdict(list)  # contig -> [(start,end,strand,partial)]
    if isinstance(row.prodigal, str) and row.prodigal:
        for l in open(row.prodigal):
            if l[0] != '>': continue
            f = l[1:].split(' # ')
            ctg = f[0].rsplit('_', 1)[0]
            part = re.search(r'partial=(\d\d)', l).group(1)
            genes[ctg].append((int(f[1]), int(f[2]), '+' if f[3] == '1' else '-', part))
    else:
        for l in gzip.open(f'{OUT}/work/prodigal/{row.genome}.gff.gz', 'rt'):
            if l[0] == '#': continue
            f = l.split('\t')
            if len(f) < 9 or f[2] != 'CDS': continue
            part = re.search(r'partial=(\d\d)', f[8]).group(1)
            genes[f[0]].append((int(f[3]), int(f[4]), f[6], part))
    for c in genes: genes[c].sort()
    return genes

def miniprot(genome):
    best = {}
    paf = None
    for l in gzip.open(f'{OUT}/work/mp/{genome}.paf.gz', 'rt'):
        if l.startswith('##PAF'):
            paf = l.rstrip('\n').split('\t')[1:]
            continue
        f = l.rstrip('\n').split('\t')
        attr = dict(x.split('=', 1) for x in f[8].split(';'))
        q = paf
        tags = dict((t[:2], t[5:]) for t in q[12:])
        og = q[0].split('|')[0]
        rec = dict(og=og, ref=q[0], qlen=int(q[1]), qs=int(q[2]), qe=int(q[3]), strand=q[4], ctg=q[5],
                   tlen=int(q[6]), ts=int(q[7]), te=int(q[8]), AS=int(tags['AS']), fs=int(tags['fs']),
                   st=int(tags['st']), intron=len(re.findall(r'\d+[NUV]', tags.get('cg', ''))),
                   ident=float(attr['Identity']), rank=int(attr['Rank']))
        if og not in best or rec['AS'] > best[og]['AS']: best[og] = rec
    return best

nOG = len(OGKEEP)
print('core groups analysed', nOG, 'min median length', MINLEN, file=sys.stderr)
loci_rows, gsum, ig_rows = [], [], []
for _, row in m.iterrows():
    genes = prodigal_genes(row)
    # intergenic spacers between consecutive Prodigal genes on a contig (contig ends excluded)
    sp = []
    for c, gl in genes.items():
        for a, b in zip(gl, gl[1:]): sp.append(b[0] - a[1] - 1)
    sp = np.array(sp)
    ncds = sum(len(v) for v in genes.values())
    best = miniprot(row.genome)
    n_hit = n_edge = 0; recs = []
    for og, r in best.items():
        if og not in OGKEEP: continue
        if r['ident'] < MIN_ID or (r['qe'] - r['qs']) / r['qlen'] < MIN_QCOV: continue
        n_hit += 1
        if r['ts'] < EDGE or r['tlen'] - r['te'] < EDGE: n_edge += 1; continue
        s, e = r['ts'] + 1, r['te']
        orfs = []
        for (a, b, strand, part) in genes.get(r['ctg'], []):
            if b < s: continue
            if a > e: break
            ov = min(b, e) - max(a, s) + 1
            if strand == r['strand'] and ov >= ORF_OV and ov >= ORF_FRAC * (b - a + 1): orfs.append((a, b, part))
        L = max([(b - a + 1) / 3 - (1 if p[1 if r['strand'] == '+' else 0] == '0' else 0) for a, b, p in orfs], default=0)
        r.update(n_orf=len(orfs), orf_ratio=L / r['qlen'], qcov=(r['qe'] - r['qs']) / r['qlen'])
        recs.append(r)
    d = pd.DataFrame(recs)
    d['genome'] = row.genome
    d['disrupt'] = (d.fs > 0) | (d.st > 0)
    d['short80'] = d.orf_ratio < 0.8
    d['split'] = d.n_orf >= 2
    d['any_decay'] = d.disrupt | d.short80 | d.split
    loci_rows.append(d)
    n = len(d)
    gsum.append(dict(genome=row.genome, core_groups=nOG, hits=n_hit, hits_at_contig_end=n_edge, assessed=n,
        pct_hit=100 * n_hit / nOG, median_identity=d.ident.median(),
        pct_frameshift=100 * (d.fs > 0).mean(), pct_internal_stop=100 * (d.st > 0).mean(),
        pct_fs_or_stop=100 * d.disrupt.mean(), pct_intron_like=100 * (d.intron > 0).mean(),
        pct_orf_lt80=100 * d.short80.mean(), median_orf_ratio=d.orf_ratio.median(),
        n_split=int(d.split.sum()), pct_split=100 * d.split.mean(),
        pct_any_decay=100 * d.any_decay.mean(), median_ref_coverage=d.qcov.median(),
        prodigal_cds=ncds, n_spacers=len(sp), spacer_median_bp=float(np.median(sp)),
        spacer_mean_bp=float(sp.mean()), pct_spacer_gt300=100 * (sp > 300).mean(),
        pct_spacer_gt500=100 * (sp > 500).mean(), pct_overlap=100 * (sp < 0).mean()))
    print(row.genome, n, round(gsum[-1]['pct_fs_or_stop'], 2), round(gsum[-1]['pct_any_decay'], 2), file=sys.stderr)
g = m.merge(pd.DataFrame(gsum), on='genome')
g['contigs_per_Mb'] = g.contigs / (g.size_bp / 1e6)
g.to_csv(OUT + '/results/per_genome_decay_screen' + SUFFIX + '.tsv', sep='\t', index=False, float_format='%.4g')
pd.concat(loci_rows).to_csv(OUT + '/results/per_locus_calls' + SUFFIX + '.tsv.gz', sep='\t', index=False, float_format='%.4g')
