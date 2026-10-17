#!/usr/bin/env python3
"""Extract the best multi-start RELAX fits drawn in Figure 3a (Table S8r, rows 'BEST').

FSM  alternative: X_FSM_alt_from_sp5null_cw.fit.bf (lnL -694848.84, spooled HyPhy fit,
                  coordinate-wise re-optimisation started from the seed-201 5-start null)
     null       : L_FSM_sp5null_cw.fit.bf          (lnL -694852.34)
LEP  alternative/null: R_LEP_s1_sp1.json (RELAX 3.1.1 JSON; lnL -694939.57 / -695030.99)

Output: inputs/relax_best_fits.tsv  (one row per rate class and branch set)
Usage : python3 extract_relax_best_fits.py [RUNS_DIR]
"""
import json
import re
import sys
from pathlib import Path

import pandas as pd

RUNS = Path(sys.argv[1] if len(sys.argv) > 1 else
            "/path/to/project/analyses/relax_stem/relax/runs")
OUT = Path(__file__).resolve().parents[1] / "inputs" / "relax_best_fits.tsv"


def bf_params(path):
    """K, reference and test classes from a spooled HyPhy RELAX fit (.bf)."""
    s = open(path).read()
    v = {m.group(1): m.group(3) for m in
         re.finditer(r"^global (relax\.[A-Za-z0-9_.]+)(:?=)([^;]+);", s, re.M)}
    f = lambda n: float(v["relax." + n])
    K = f("K")
    om = [f(f"reference.omega{i}") for i in (1, 2, 3)]
    a0, a1 = f("reference.bsrel_mixture_aux_0"), f("reference.bsrel_mixture_aux_1")
    w = [a0, (1 - a0) * a1, (1 - a0) * (1 - a1)]
    cls = sorted(zip(om, [o ** K for o in om], w))
    return K, cls


def lnl_from_log(path, pat):
    txt = open(path).read()
    return float(re.findall(pat, txt)[-1])


rows = []
# ---- FSM
K, cls = bf_params(RUNS / "X_FSM_alt_from_sp5null_cw.fit.bf")
ln_alt = lnl_from_log(RUNS / "X_FSM_alt_from_sp5null_cw.log", r"PASS\t1\tlnL\s+(-?[\d.]+)")
ln_null = lnl_from_log(RUNS / "L_FSM_sp5null_cw.log", r"PASS\t1\tlnL\s+(-?[\d.]+)")
for i, (w0, w1, p) in enumerate(cls):
    rows.append(dict(test="FSM", fit="X_FSM_alt_from_sp5null_cw.fit.bf", cls=i, omega_ref=w0,
                     omega_test=w1, proportion=p, K=K, LRT=2 * (ln_alt - ln_null),
                     lnL_alt=ln_alt, lnL_null=ln_null))
# ---- LEP
j = json.load(open(RUNS / "R_LEP_s1_sp1.json"))
ref = j["fits"]["RELAX alternative"]["Rate Distributions"]["Reference"]
tst = j["fits"]["RELAX alternative"]["Rate Distributions"]["Test"]
for i, c in enumerate(sorted(ref, key=lambda c: ref[c]["omega"])):
    rows.append(dict(test="LEP", fit="R_LEP_s1_sp1.json", cls=i, omega_ref=ref[c]["omega"],
                     omega_test=tst[c]["omega"], proportion=ref[c]["proportion"],
                     K=j["test results"]["relaxation or intensification parameter"],
                     LRT=j["test results"]["LRT"],
                     lnL_alt=j["fits"]["RELAX alternative"]["Log Likelihood"],
                     lnL_null=j["fits"]["RELAX null"]["Log Likelihood"]))
df = pd.DataFrame(rows)
df.to_csv(OUT, sep="\t", index=False)
print(df.round(5).to_string())
