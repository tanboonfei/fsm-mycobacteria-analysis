#!/usr/bin/env python3
"""Additional checks (Table S14 and S11 notes).
1. Table S14: re-run the FSM enrichment without the four FSM MAGs carrying a non-mycobacterial 16S rRNA gene,
   and count named FSM species carrying each enriched family.
2. Table S11: opsin prevalence counted per GTDB species instead of per genome.
Inputs: gene_presence_by_genome.pkl (per-genome eggNOG gene names, from the gene-family step),
Additional_file_2.xlsx (Table S1) and scripts/inputs/tsv/S18_FSM_enriched_genes.tsv.
Run from the project root: python3 table_checks_S11_S14.py <gene_presence.pkl> <Additional_file_2.xlsx> <S18.tsv>"""
import sys, pickle, collections, openpyxl, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
pkl, xlsx, s18 = sys.argv[1:4]
gname = pickle.load(open(pkl, "rb"))["gname"]
rows = list(openpyxl.load_workbook(xlsx, read_only=True)["S1"].iter_rows(values_only=True))
I = {k: i for i, k in enumerate(rows[0])}
sp = {r[I["genome_id_or_accession"]]: r[I["species"]] for r in rows[1:] if r[0]}
cl = {r[I["genome_id_or_accession"]]: str(r[I["clade"]])[:3] for r in rows[1:] if r[0]}
enr = pd.read_csv(s18, sep="\t"); enr = enr[enr.pct_FSM.notna()]
bad = {"spire_mag_01108461", "spire_mag_01753944", "spire_mag_01960731", "spire_mag_02102612"}
fsm = [g for g in gname if cl.get(g) == "FSM" and g not in bad]
oth = [g for g in gname if cl.get(g) in ("RGM", "SGM")]
res = []
for gene in set().union(*gname.values()):
    a = sum(gene in gname[g] for g in fsm); b = sum(gene in gname[g] for g in oth)
    res.append((gene, stats.fisher_exact([[a, len(fsm) - a], [b, len(oth) - b]], alternative="greater")[1],
                100 * a / len(fsm), 100 * b / len(oth)))
d = pd.DataFrame(res, columns=["gene", "p", "pct_FSM", "pct_other"])
d["q"] = multipletests(d.p, method="fdr_bh")[1]
sel = d[(d.q < 0.05) & (d.pct_FSM >= 50) & (d.pct_FSM - d.pct_other >= 30)]
print("enriched without the four MAGs:", len(sel), "| kept:", len(set(sel.gene) & set(enr.gene)), "of", len(enr),
      "| lost:", set(enr.gene) - set(sel.gene), "| gained:", set(sel.gene) - set(enr.gene))
fsp = collections.defaultdict(list)
for g in gname:
    if cl.get(g) == "FSM" and sp[g] != "novel sp.": fsp[sp[g]].append(g)
k = [sum(any(gene in gname[g] for g in v) for v in fsp.values()) for gene in enr.gene]
print(f"named FSM species in panel: {len(fsp)}; species per enriched family min {min(k)}, median {sorted(k)[len(k)//2]}")
for c in ("FSM", "RGM", "SGM"):
    s = collections.defaultdict(list)
    for r in rows[1:]:
        if r[0] and str(r[I["clade"]]).startswith(c) and r[I["species"]] != "novel sp.":
            s[r[I["species"]]].append(r[I["rhodopsin"]] == "present")
    print(f"{c}: opsin in {sum(any(v) for v in s.values())} of {len(s)} species")
