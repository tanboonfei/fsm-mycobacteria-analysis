"""Print K and omega distributions (reference/test/unclassified) stored in a spooled HyPhy RELAX fit (.bf)."""
import re, sys
for f in sys.argv[1:]:
    s = open(f).read(); v = {}
    for m in re.finditer(r'^global (relax\.[A-Za-z0-9_.]+)(:?=)([^;]+);', s, re.M):
        v[m.group(1)] = (m.group(2), m.group(3))
    def val(n):
        op, x = v[n]
        try: return float(x)
        except ValueError: return x
    K = val('relax.K'); K = float(K) if not isinstance(K, str) else K
    def dist(ns, Kexp=None):
        om = [val('relax.%s.omega%d' % (ns if ns != 'test' else 'reference', i)) for i in (1, 2, 3)]
        a0 = val('relax.%s.bsrel_mixture_aux_0' % (ns if ns != 'test' else 'reference')); a1 = val('relax.%s.bsrel_mixture_aux_1' % (ns if ns != 'test' else 'reference'))
        if Kexp is not None: om = [o ** Kexp for o in om]
        w = [a0, (1 - a0) * a1, (1 - a0) * (1 - a1)]
        return '; '.join('%.4g (%.2f%%)' % (o, 100 * x) for o, x in sorted(zip(om, w)))
    k = K if isinstance(K, float) else 1.0
    print('%s\tK=%s\tref: %s\ttest: %s\tunclassified: %s' % (f.split('/')[-1], K, dist('reference'), dist('test', k), dist('unclassified')))
