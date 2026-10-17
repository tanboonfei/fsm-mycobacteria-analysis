#!/usr/bin/env python3
"""Define single-copy core orthologue groups from the 862-genome MMseqs2 protein
clusters (easy-cluster --min-seq-id 0.5 -c 0.8 --cluster-mode 2; mmclust862/),
restricted to 12 closed reference genomes (6 RGM, 6 SGM; none of them a negative
control, none in the leprosy clade). A group is kept when every reference genome
contributes exactly one PGAP protein. Writes the reference proteins (one per
reference per group) with headers OG|genome|protein."""
import collections, sys, os
ROOT = '/path/to/project'
OUT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REFS = {  # accession: (clade, species)
 'GCF_010731875.1': ('RGM', 'M. sarraceniae (42-genome sister clade)'),
 'GCF_001457595.1': ('RGM', 'M. smegmatis'),
 'GCF_900637195.1': ('RGM', 'M. aurum'),
 'GCF_900187065.1': ('RGM', 'M. thermoresistibile'),
 'GCF_010731575.1': ('RGM', 'M. monacense'),
 'GCF_016313185.1': ('RGM', 'M. baixiangningiae'),
 'GCF_000340435.2': ('SGM', 'M. haemophilum'),
 'GCF_003391395.1': ('SGM', 'M. marinum'),
 'GCF_009741445.1': ('SGM', 'M. avium'),
 'GCF_010727125.1': ('SGM', 'M. terrae'),
 'GCF_010727605.1': ('SGM', 'M. simiae'),
 'GCF_009936235.1': ('SGM', 'M. xenopi'),
}
cnt = collections.defaultdict(lambda: collections.defaultdict(list))
with open(ROOT + '/mmclust862/c862_cluster.tsv') as f:
    for line in f:
        rep, mem = line.rstrip('\n').split('\t')
        g, p = mem.split('|', 1)
        if g in REFS:
            cnt[rep][g].append(p)
keep = {rep: d for rep, d in cnt.items()
        if len(d) == len(REFS) and all(len(v) == 1 for v in d.values())}
print('clusters touching refs', len(cnt), 'single-copy core', len(keep), file=sys.stderr)
want = collections.defaultdict(dict)  # genome -> protein -> OG
ogs = sorted(keep)
ogname = {rep: 'OG%04d' % (i + 1) for i, rep in enumerate(ogs)}
rep_of = {v: k for k, v in ogname.items()}
for rep in ogs:
    for g, v in keep[rep].items():
        want[g][v[0]] = ogname[rep]
n = 0
with open(OUT + '/work/core_ref_proteins.faa', 'w') as out, open(OUT + '/results/orthogroups.tsv', 'w') as tab:
    tab.write('OG\tmmseqs_rep\tgenome\tprotein\tlength_aa\n')
    for g in REFS:
        seqs = {}; cur = None
        for line in open(f'{ROOT}/refprot/ncbi_dataset/data/{g}/protein.faa'):
            if line.startswith('>'):
                cur = line[1:].split()[0]; seqs[cur] = [] if cur in want[g] else None
            elif seqs.get(cur) is not None:
                seqs[cur].append(line.strip())
        for p, og in want[g].items():
            s = ''.join(seqs[p]).rstrip('*')
            out.write(f'>{og}|{g}|{p}\n{s}\n'); n += 1
            tab.write(f'{og}\t{rep_of[og]}\t{g}\t{p}\t{len(s)}\n')
print('wrote', n, 'proteins', file=sys.stderr)
with open(OUT + '/results/reference_genomes.tsv', 'w') as f:
    f.write('accession\tclade\tspecies\n')
    for g, (c, sp) in REFS.items(): f.write(f'{g}\t{c}\t{sp}\n')
