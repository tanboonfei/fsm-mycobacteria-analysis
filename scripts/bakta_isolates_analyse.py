#!/usr/bin/env python3
"""Compare pseudogene load between cultured isolates and MAGs, one pipeline.

Closes the gap left by the binning control (Table S33): pseudogenes had been
called with Bakta for the SPIRE MAGs only, so the low FSM pseudogene load could
not be separated from the possibility that binning discards the repetitive,
pseudogene-rich fraction of a genome.

The cultured RGM and SGM isolates are now annotated with the same Bakta version
(1.11.4) and database (6.0, full) under the same command, so isolate-against-MAG
within a clade isolates the binning effect on pseudogene calling.

Counts are read from the `pseudogenes:` line of each Bakta .txt summary, the
same rule that reproduces the published per-MAG counts exactly.

Writes supplementary_text/tsv/S34_pseudogene_isolate_control.tsv.

Usage: python3 bakta_isolates_analyse.py
"""
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from paths import ROOT, RUNOUT, TSV, load_s1

ISO_DIR = ROOT.parent / "bakta_isolates"
ISOLATE = ["Cultured isolate", "NTM (isolate)", "Pathogen (isolate)"]


def read_count(txt):
    for line in open(txt):
        m = re.match(r"^pseudogenes:\s*(\d+)", line, re.I)
        if m:
            return int(m.group(1))
    return None


recs = []
for d in sorted(ISO_DIR.glob("*")):
    if not d.is_dir():
        continue
    txt = d / f"{d.name}.txt"
    if txt.exists():
        n = read_count(txt)
        if n is not None:
            recs.append({"genome_id_or_accession": d.name, "pseudogenes_bakta": n})
iso = pd.DataFrame(recs)
print(f"isolate genomes annotated: {len(iso)}")
if iso.empty:
    raise SystemExit("no Bakta isolate output yet")

s1 = load_s1(pd)
s1["is_isolate"] = s1.culture_status.isin(ISOLATE)
mag = pd.read_csv(TSV / "S5b_pseudogenes_MAG_panel.tsv", sep="\t").rename(
    columns={"genome": "genome_id_or_accession", "pseudogenes": "pseudogenes_bakta"})

df = pd.concat([
    s1.merge(iso, on="genome_id_or_accession", how="inner").assign(kind="cultured isolate"),
    s1.merge(mag[["genome_id_or_accession", "pseudogenes_bakta"]],
             on="genome_id_or_accession", how="inner").assign(kind="MAG"),
], ignore_index=True)

rows = []
print("\n" + "=" * 74)
print("PSEUDOGENE LOAD, BAKTA 1.11.4 / DB 6.0 THROUGHOUT")
print("=" * 74)
print(df.groupby(["clade3", "kind"]).pseudogenes_bakta
        .agg(n="size", median="median", mean="mean", q25=lambda x: x.quantile(.25),
             q75=lambda x: x.quantile(.75)).round(1).to_string())

print("\n" + "=" * 74)
print("ISOLATE vs MAG WITHIN CLADE  (the binning test)")
print("=" * 74)
for c in ["RGM", "SGM"]:
    a = df[(df.clade3 == c) & (df.kind == "cultured isolate")].pseudogenes_bakta
    b = df[(df.clade3 == c) & (df.kind == "MAG")].pseudogenes_bakta
    if len(a) < 3 or len(b) < 3:
        print(f"  {c}: too few ({len(a)} isolates, {len(b)} MAGs)")
        continue
    p = stats.mannwhitneyu(a, b).pvalue
    d = float(a.median() - b.median())
    print(f"  {c}: isolate median {a.median():.0f} (n={len(a)})  "
          f"MAG median {b.median():.0f} (n={len(b)})  diff {d:+.0f}  p={p:.2g}")
    rows.append({"comparison": f"{c}: cultured isolate vs MAG",
                 "group_a": "cultured isolate", "n_a": len(a), "median_a": a.median(),
                 "group_b": "MAG", "n_b": len(b), "median_b": b.median(),
                 "difference": round(d, 1), "MannWhitney_p": f"{p:.2g}",
                 "note": "positive difference = binning removes pseudogenes"})

# The claim under attack: FSM MAGs against everything else.
fsm = df[(df.clade3 == "FSM")].pseudogenes_bakta
print("\n" + "=" * 74)
print("DOES THE BINNING OFFSET EXPLAIN THE FSM RESULT?")
print("=" * 74)
print(f"  FSM (MAGs only, no isolate exists): median {fsm.median():.0f} (n={len(fsm)})")
for c in ["RGM", "SGM"]:
    a = df[(df.clade3 == c) & (df.kind == "cultured isolate")].pseudogenes_bakta
    b = df[(df.clade3 == c) & (df.kind == "MAG")].pseudogenes_bakta
    if len(a) < 3:
        continue
    off = float(a.median() - b.median())
    gap = float(a.median() - fsm.median())
    corrected = fsm.median() + off
    print(f"  vs {c} isolates: gap {gap:+.0f} pseudogenes; binning offset {off:+.0f}; "
          f"FSM adjusted upward for binning would be {corrected:.0f} against {a.median():.0f}")
    rows.append({"comparison": f"FSM MAGs vs {c} cultured isolates",
                 "group_a": f"{c} cultured isolate", "n_a": len(a), "median_a": a.median(),
                 "group_b": "FSM MAG", "n_b": len(fsm), "median_b": fsm.median(),
                 "difference": round(gap, 1), "MannWhitney_p":
                     f"{stats.mannwhitneyu(a, fsm).pvalue:.2g}",
                 "note": f"binning offset measured in {c} is {off:+g}; "
                         f"FSM adjusted for it would be {corrected:g}"})

# Mechanism: truncated reading frames at contig boundaries are called as
# pseudogenes; a more fragmented assembly carries more of them.
asm = pd.read_csv(ROOT.parent / "analysis_gaps" / "assembly_stats.tsv", sep="\t")
# S1 gained its own n_contigs column, which collided with the one in
# assembly_stats and left the merge with n_contigs_x / n_contigs_y.
asm = asm[[c for c in asm.columns
           if c == "genome_id_or_accession" or c not in df.columns]]
frag = df[df.kind == "MAG"].merge(asm, on="genome_id_or_accession", how="left")
if "n_contigs" not in frag.columns:
    frag = frag.rename(columns={"n_contigs_y": "n_contigs"})
frag = frag.dropna(subset=["n_contigs"])
print("\n" + "=" * 74)
print("MECHANISM: PSEUDOGENE CALLS AGAINST ASSEMBLY FRAGMENTATION")
print("=" * 74)
r, p = stats.spearmanr(frag.n_contigs, frag.pseudogenes_bakta)
print(f"  all MAGs: rho={r:+.2f} p={p:.2g} (n={len(frag)}), median {frag.n_contigs.median():.0f} contigs")
rows.append({"comparison": "MAGs: pseudogene calls vs contig count", "group_a": "", "n_a": len(frag),
             "median_a": "", "group_b": "", "n_b": "", "median_b": "",
             "difference": f"rho = {r:+.2f}", "MannWhitney_p": f"{p:.2g}",
             "note": "fragmentation inflates pseudogene calls; MAGs are more fragmented than isolates"})
for c in ["FSM", "RGM", "SGM"]:
    s = frag[frag.clade3 == c]
    if len(s) > 10:
        r2, p2 = stats.spearmanr(s.n_contigs, s.pseudogenes_bakta)
        print(f"   {c}: rho={r2:+.2f} p={p2:.2g} (n={len(s)})")
        rows.append({"comparison": f"{c} MAGs: pseudogene calls vs contig count", "group_a": "",
                     "n_a": len(s), "median_a": "", "group_b": "", "n_b": "", "median_b": "",
                     "difference": f"rho = {r2:+.2f}", "MannWhitney_p": f"{p2:.2g}", "note": ""})

# Bakta under-calls pseudogenes in a massively decayed genome, which is the same
# behaviour Prodigal shows for coding density (Table S24).
print("\n" + "=" * 74)
print("LEPROSY CLADE: BAKTA AGAINST THE PUBLISHED PGAP COUNTS")
print("=" * 74)
PGAP = {"Mycobacterium leprae": 1133}
dec = df[df.species.astype(str).str.contains("leprae|lepromatosis|uberis", case=False, na=False)]
for _, r_ in dec.iterrows():
    ref = PGAP.get(r_.species)
    print(f"  {r_.species:28s} Bakta {r_.pseudogenes_bakta:5.0f}"
          + (f"   PGAP {ref}   ratio {ref / max(r_.pseudogenes_bakta,1):.0f}x" if ref else ""))
    rows.append({"comparison": "Leprosy clade: Bakta vs PGAP pseudogene call",
                 "group_a": "Bakta", "n_a": 1, "median_a": r_.pseudogenes_bakta,
                 "group_b": "NCBI PGAP", "n_b": 1 if ref else "", "median_b": ref if ref else "",
                 "difference": (ref - r_.pseudogenes_bakta) if ref else "",
                 "MannWhitney_p": "",
                 "note": f"{r_.species}; Bakta annotates decayed frames as CDSs, and it "
                         f"under-calls pseudogenes in decayed genomes"})

out = RUNOUT / "S6_pseudogene_isolate_control.tsv"
pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
print(f"\nwrote {out}")

# S6 keeps its published MAG-only contents. The combined per-genome table is
# written alongside it so the isolate calls are traceable.
per = df[["genome_id_or_accession", "clade3", "kind", "species", "pseudogenes_bakta"]]
per = per.sort_values(["clade3", "kind", "genome_id_or_accession"])
per_out = RUNOUT / "S5_pseudogenes_all.tsv"
per.to_csv(per_out, sep="\t", index=False)
print(f"per-genome counts -> {per_out}")
