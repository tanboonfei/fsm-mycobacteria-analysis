#!/usr/bin/env python3
"""Test whether MAG binning could manufacture the FSM reduction signal.

The attack: binning preferentially discards low-coverage repetitive regions,
which are disproportionately non-coding and pseudogene-rich. One process would
then make genomes look smaller, denser and cleaner at once, which is three of
the four headline observations.

The obvious control is unavailable. FSM contains no cultured isolate: all 134
genomes are metagenome-assembled, 96 recovered here from SPIRE and 38 already
deposited in GTDB. The `origin` column separates the two catalogues, not
binned from unbinned, and `culture_status` is "MAG / uncultured" for every
FSM genome.

So the artefact is measured where it can be measured. RGM and SGM each contain
both cultured isolates and MAGs, so the isolate-minus-MAG offset within those
clades estimates the size and direction of any binning effect on the same
metrics, under the same pipeline. That offset is then compared with the
FSM-versus-RGM difference the manuscript reports.

Two comparisons are reported per clade:
    isolate vs MAG      the binning estimate (RGM, SGM)
    GTDB vs SPIRE MAG   catalogue/assembly effect among MAGs (all three clades)

Genome size and GC are annotation-independent. Coding density is not: references
are PGAP-annotated and MAGs Prodigal-annotated, and the two agree to about one
point outside pseudogene-rich genomes (Table S24). Coding-density offsets below
roughly one point therefore cannot be separated from annotation method, and are
flagged in the output.

Writes supplementary_text/tsv/S33_binning_control.tsv.

Usage: python3 binning_artefact_control.py
"""
import numpy as np
import pandas as pd
from scipy import stats

from paths import ORDER, RUNOUT, TSV, load_s1

# (column, label, decimals, annotation-independent?, MAG subset for the offset)
METRICS = [("genome_size_Mb", "Genome size (Mb)", 2, True, None),
           ("GC_pct", "GC content (%)", 2, True, None),
           ("coding_density_pct", "Coding density (%)", 2, False, None)]

# Cultured isolates are annotated by PGAP and MAGs by Prodigal, so the coding
# density of an isolate is not comparable with that of a GTDB-deposited MAG.
# The isolate-vs-MAG offset therefore reads the uniform Prodigal recall for the
# isolates and restricts the MAG side to the SPIRE bins, which Prodigal
# annotated natively. Genome size and GC need no such care.
METRICS_ISO = [("genome_size_Mb", "Genome size (Mb)", 2, True, None),
               ("GC_pct", "GC content (%)", 2, True, None),
               ("coding_density_pct_uniform_prodigal",
                "Coding density (%), uniform Prodigal", 2, False, "SPIRE MAG")]
ANNOTATION_FLOOR = 1.0   # points of coding density; see Table S24
SHORT = {"genome_size_Mb": "genome size", "GC_pct": "GC",
         "coding_density_pct": "coding density",
         "coding_density_pct_uniform_prodigal": "coding density"}

s1 = load_s1(pd)
s1["is_isolate"] = s1.culture_status.isin(
    ["Cultured isolate", "NTM (isolate)", "Pathogen (isolate)"])

rows = []


def compare(label, a, b, name_a, name_b, note="", metrics=None):
    for col, pretty, dp, clean, subset in (metrics or METRICS):
        bb = b[b.origin == subset] if subset else b
        x = a[col].dropna().values
        y = bb[col].dropna().values
        if len(x) < 3 or len(y) < 3:
            continue
        p = stats.mannwhitneyu(x, y).pvalue
        diff = float(np.median(x) - np.median(y))
        flag = note
        if not clean and abs(diff) < ANNOTATION_FLOOR:
            flag = (flag + "; " if flag else "") + \
                "within the PGAP/Prodigal annotation floor, not separable from method"
        if subset:
            flag += ("; the isolate column is the uniform Prodigal recall and the MAG "
                     f"column is restricted to the {subset.split()[0]} bins, which Prodigal "
                     "annotated natively, so the two sides carry the same gene caller")
        rows.append({
            "comparison": label, "metric": pretty,
            f"n_{name_a}": len(x), f"median_{name_a}": round(float(np.median(x)), dp),
            f"n_{name_b}": len(y), f"median_{name_b}": round(float(np.median(y)), dp),
            "difference": round(diff, dp), "MannWhitney_p": f"{p:.2g}", "note": flag,
        })


print("=" * 78)
print("COMPOSITION")
print("=" * 78)
comp = pd.crosstab(s1[s1.clade3.isin(ORDER)].clade3,
                   np.where(s1[s1.clade3.isin(ORDER)].is_isolate, "cultured isolate", "MAG"))
print(comp.to_string())
print("\nFSM cultured isolates:", int(s1[(s1.clade3 == 'FSM')].is_isolate.sum()),
      "-> the within-FSM isolate control does not exist")

# ---- 1. the binning estimate, where isolates and MAGs coexist
print("\n" + "=" * 78)
print("1. ISOLATE vs MAG WITHIN CLADE  (estimates the binning offset)")
print("=" * 78)
for c in ["RGM", "SGM"]:
    sub = s1[s1.clade3 == c]
    iso, mag = sub[sub.is_isolate], sub[~sub.is_isolate]
    compare(f"{c}: cultured isolate vs MAG", iso, mag, "isolate", "MAG",
            metrics=METRICS_ISO)
    for col, pretty, dp, _, subset in METRICS_ISO:
        m = mag[mag.origin == subset] if subset else mag
        d = np.median(iso[col].dropna()) - np.median(m[col].dropna())
        print(f"  {c} {pretty:22s} isolate {np.median(iso[col].dropna()):7.2f}  "
              f"MAG {np.median(m[col].dropna()):7.2f}  diff {d:+7.2f}")

# ---- 2. catalogue effect among MAGs, including within FSM
print("\n" + "=" * 78)
print("2. GTDB-DEPOSITED vs SPIRE MAGs  (both binned; assembly/pipeline effect)")
print("=" * 78)
for c in ORDER:
    sub = s1[(s1.clade3 == c) & (~s1.is_isolate)]
    g, sp = sub[sub.origin == "GTDB reference"], sub[sub.origin == "SPIRE MAG"]
    if len(g) < 3:
        print(f"  {c}: only {len(g)} GTDB-deposited MAGs, skipped")
        continue
    compare(f"{c}: GTDB-deposited MAG vs SPIRE MAG", g, sp, "gtdb", "spire",
            note="both sets are metagenome-assembled")
    for col, pretty, dp, _, _sub in METRICS:
        d = np.median(g[col].dropna()) - np.median(sp[col].dropna())
        print(f"  {c} {pretty:22s} GTDB {np.median(g[col].dropna()):7.2f}  "
              f"SPIRE {np.median(sp[col].dropna()):7.2f}  diff {d:+7.2f}  (n={len(g)}/{len(sp)})")

# ---- 3. does the artefact scale with how much is missing?
print("\n" + "=" * 78)
print("3. WITHIN FSM MAGs: does the signal track assembly quality?")
print("=" * 78)
fsm = s1[(s1.clade3 == "FSM") & s1.completeness_pct.notna()]
for col, pretty, dp, _, _sub in METRICS:
    v = fsm[[col, "completeness_pct"]].dropna()
    r, p = stats.spearmanr(v.completeness_pct, v[col])
    print(f"  {pretty:22s} vs completeness  rho={r:+.2f}  p={p:.2g}  (n={len(v)})")
    rows.append({"comparison": "FSM MAGs: metric vs estimated completeness",
                 "metric": pretty, "n_isolate": len(v), "median_isolate": "",
                 "n_MAG": "", "median_MAG": "", "difference": f"rho = {r:+.2f}",
                 "MannWhitney_p": f"{p:.2g}",
                 "note": "binning loss should make less complete MAGs smaller and denser"})

hi = fsm[fsm.completeness_pct >= 99]
print(f"\n  FSM MAGs at >=99% complete (n={len(hi)}): "
      f"size {hi.genome_size_Mb.median():.2f} Mb, "
      f"coding {hi.coding_density_pct.median():.2f}%, GC {hi.GC_pct.median():.2f}%")
print(f"  all FSM MAGs           (n={len(fsm)}): "
      f"size {fsm.genome_size_Mb.median():.2f} Mb, "
      f"coding {fsm.coding_density_pct.median():.2f}%, GC {fsm.GC_pct.median():.2f}%")
rows.append({"comparison": "FSM MAGs ≥99% complete vs all FSM MAGs",
             "metric": "Genome size (Mb)", "n_isolate": len(hi),
             "median_isolate": round(float(hi.genome_size_Mb.median()), 2),
             "n_MAG": len(fsm), "median_MAG": round(float(fsm.genome_size_Mb.median()), 2),
             "difference": round(float(hi.genome_size_Mb.median() - fsm.genome_size_Mb.median()), 2),
             "MannWhitney_p": "", "note": "the most complete FSM MAGs are not larger"})

# ---- 4. the comparison that matters
print("\n" + "=" * 78)
print("4. ARTEFACT SIZE AGAINST THE REPORTED EFFECT")
print("=" * 78)
rgm = s1[s1.clade3 == "RGM"]
fsm_all = s1[s1.clade3 == "FSM"]
for col, pretty, dp, _, subset in METRICS_ISO:
    effect = float(np.median(fsm_all[col].dropna()) - np.median(rgm[col].dropna()))
    art = []
    for c in ["RGM", "SGM"]:
        sub = s1[s1.clade3 == c]
        m = sub[~sub.is_isolate]
        if subset:
            m = m[m.origin == subset]
        art.append(float(np.median(sub[sub.is_isolate][col].dropna())
                         - np.median(m[col].dropna())))
    worst = max(abs(a) for a in art)
    ratio = abs(effect) / worst if worst else float("inf")
    print(f"  {pretty:22s} FSM-RGM {effect:+7.2f} | largest binning offset "
          f"{worst:6.2f} | effect is {ratio:5.1f}x larger")
    rows.append({"comparison": "FSM vs RGM effect against the largest binning offset",
                 "metric": pretty, "n_isolate": "", "median_isolate": "",
                 "n_MAG": "", "median_MAG": "",
                 "difference": f"{effect:+.2f} vs {worst:.2f}",
                 "MannWhitney_p": "",
                 "note": (f"reported effect is {ratio:.1f}x the largest "
                          f"isolate-vs-MAG offset") if ratio >= 1 else
                         (f"FSM sits {abs(effect):.2f} points below RGM on exact "
                          f"{SHORT[col]}; the largest isolate-vs-MAG offset is "
                          f"{worst:.2f}, and no {SHORT[col]} shift is resolvable")})

out = pd.DataFrame(rows)
dest = RUNOUT / "S4_binning_control.tsv"
out.to_csv(dest, sep="\t", index=False)
print(f"\nwrote {dest}  ({len(out)} rows)")
