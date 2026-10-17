"""Make a saved RELAX LF file for another labelling by rewriting the per-branch model annotations in the Tree line.
usage: relabel.py <saved fit .bf> <labelled newick> <out .bf>"""
import re, sys, collections
fit, nwk, out = sys.argv[1:]
lab = dict(re.findall(r'([A-Za-z0-9_]+)\{(Reference|Test|Unclassified)\}', open(nwk).read()))
lines = open(fit).read().split('\n')
idx = [i for i, l in enumerate(lines) if l.startswith('Tree ')]
assert len(idx) == 1
seen = collections.Counter(); missing = []
def sub(m):
    n = m.group(1)
    if n not in lab: missing.append(n); return m.group(0)
    seen[lab[n]] += 1
    return '%s{relax.%s}' % (n, lab[n].lower())
lines[idx[0]] = re.sub(r'([A-Za-z0-9_]+)\{relax\.[a-z]+\}', sub, lines[idx[0]])
print('relabelled', dict(seen), 'missing', missing, file=sys.stderr)
open(out, 'w').write('\n'.join(lines))
