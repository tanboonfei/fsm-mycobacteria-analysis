#!/usr/bin/env python3
"""Supplementary Figure S3: candidate FSM-associated gene families by functional category.

Inputs : inputs/tsv/S13_FSM_enriched_genes.tsv   (Table S13: 99 candidate named families + COG category)
         inputs/tsv/S19_discriminating_families.tsv (2,477 clade-discriminating MMseqs2 clusters,
                                                  prevalence over 134 FSM / 389 RGM / 299 SGM genomes;
                                                  old table number S19, shipped with the code)
         inputs/cog_corrections.tsv   (columns gene, old_category, new_category, evidence;
                                       COG-category corrections; if absent a built-in
                                       fallback is used, see COG_FALLBACK, with a warning)
Output : supplementary_figures/FigureS3_gene_signature.{png,pdf}
Usage  : python3 build_figureS3_gene_signature.py

Selection rule (as implemented)
  1. Table S13 families are matched to S19 clusters BY GENE NAME (inner merge). Only 59 of
     the 99 names occur in S19; the other 40 are never evaluated.
  2. Several clusters can carry one name (cfp2, lpqC, asnC, nicT, cysE); the per-clade
     prevalence of a name is the MAXIMUM over its clusters, taken separately for SGM, RGM
     and FSM. None of the 25 families drawn is affected.
  3. A family qualifies if prev_FSM - max(prev_RGM, prev_SGM) >= 30 points (33 of the 59).
  4. The three largest gaps per COG category are drawn (25 families, 14 categories).
Selection and display use the same genome panel; Table S13 differs (96 SPIRE FSM MAGs, RGM+SGM
pooled as the comparator, one named family at a time).

Labels. The displayed descriptions are short labels checked against the representative's
Pfam / KEGG-KO / COG evidence in the eggNOG-mapper output (disc.emapper.annotations) and the
KEGG KO definitions; they replace raw eggNOG free text (LABELS below, one entry per family).
"""
import os
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from figstyle import INK, INPUTS, TSV, WIDTH_FULL, save, seq_cmap

TOP_PER_CAT = 3
STRICT_GAP = 30.0

s7 = pd.read_csv(TSV / "S13_FSM_enriched_genes.tsv", sep="\t")
s7 = s7[~s7.gene.astype(str).str.contains(" ")].copy()      # drop the trailing "Test universe" note row
s4 = pd.read_csv(TSV / "S19_discriminating_families.tsv", sep="\t")

# Per-clade prevalence for every named family (max over same-named clusters).
prev = (s4[s4.gene_name.notna()]
        .assign(gene=lambda d: d.gene_name.astype(str))
        .groupby("gene")[["prev_SGM", "prev_RGM", "prev_FSM"]].max() * 100)
desc = (s4[s4.gene_name.notna()]
        .assign(gene=lambda d: d.gene_name.astype(str))
        .drop_duplicates("gene").set_index("gene")["description"])

# Short verified labels (eggNOG symbol kept as the row name). Evidence per label:
#   Pfam / KO of the cluster representative (disc.emapper.annotations); KO names from KEGG.
LABELS = {
    "pnuC": "nicotinamide mononucleotide transporter (KO K03811)",
    "pap": "PPK2-family polyphosphate kinase (Pfam PPK2; symbol 'pap' is eggNOG's)",
    "mpt": "fasciclin-domain protein (Pfam Fasciclin)",
    "ptsG": "glucose-family PTS permease, single-chain EIICBA (KO K02757)",
    "cutC": "copper homeostasis protein CutC (KO K06201)",
    "crtD": "carotenoid desaturase, CrtI/CrtD family (FAD-dependent; KO K09845)",
    "deoD_1": "nucleoside phosphorylase, PNP/UDP family (KO K00757, uridine phosphorylase)",
    "glpF": "MIP-family glycerol/water channel (KO K02440, glycerol uptake facilitator)",
    "ppaC": "Mn-dependent inorganic pyrophosphatase (KO K15986)",
    "galE": "UDP-glucose 4-epimerase (KO K01784)",
    "yeeR": "membrane methyltransferase, PEMT/ICMT-like (Pfam PEMT; KO K21310 mddA)",
    "sadh": "zinc-dependent alcohol dehydrogenase (Pfam ADH_N + ADH_zinc_N)",
    "lppN": "uncharacterised protein (symbol lppN from eggNOG; no domain, KO or COG hit)",
    "nhaP": "Na+/H+ exchanger family (Pfam Na_H_Exchanger)",
    "cyp147G1": "cytochrome P450 (Pfam p450)",
    "citB": "CitB-family two-component response regulator (KO K02475)",
    "mftF": "family-2 glycosyltransferase, mycofactocin-associated (Pfam Glycos_transf_2)",
    "citE_2": "HpcH/HpaI-family aldolase, citrate lyase beta-subunit type (KO K01644)",
    "crtL": "lycopene beta-cyclase (Pfam Lycopene_cycl; KO K06443)",
    "lys1": "saccharopine dehydrogenase (KO K00290; Pfam Sacchrp_dh)",
    "tesA": "thioesterase (Pfam Thioesterase)",
    "fadD11_1": "long-chain-fatty-acid--CoA ligase (KO K01897; Pfam AMP-binding)",
    "chaA": "Ca2+:H+ (Na+/Ca2+) antiporter (KO K07300)",
    "aguA": "agmatine deiminase (KO K10536)",
    "nagB": "glucosamine-6-phosphate deaminase (KO K02564)",
    # qualifying but not drawn (kept so a rebuild with another top-N still has a label)
    "PE36": "PE-family protein (Pfam PE)",
    "crtT": "methyltransferase (Pfam Methyltransf_25 / Ubie_methyltran)",
    "cry1": "photolyase/cryptochrome family (Pfam DPRP + FAD_binding_7)",
    "eccD2": "T7SS ubiquitin-like protein (Pfam YukD)",
    "lppR_1": "PknH-like extracellular-domain protein (Pfam PknH_C)",
    "lpqC": "esterase (KO K03932; Pfam Esterase_phd)",
    "papA5": "PapA5-type acyltransferase (Pfam Condensation + PapA_C)",
    "ybbC": "DUF1343 domain protein",
}

# COG-category corrections: read from cog_corrections.tsv, else this fallback.
COG_FALLBACK = {"crtL": "Secondary metabolites"}
COG_TSV = Path(os.environ.get("MYCO_COG_TSV", INPUTS / "cog_corrections.tsv"))
if COG_TSV.exists():
    _c = pd.read_csv(COG_TSV, sep="\t")
    CAT_FIX = dict(zip(_c.gene, _c.new_category))
    CAT_FIX = {**COG_FALLBACK, **CAT_FIX}   # crtL is not in the TSV
    print(f"COG corrections: {COG_TSV} ({len(CAT_FIX)} genes)")
else:
    CAT_FIX = COG_FALLBACK
    print(f"!! {COG_TSV} not found; using built-in fallback {COG_FALLBACK}", file=sys.stderr)
s7["COG_category"] = s7.apply(
    lambda r: CAT_FIX.get(r.gene, r.COG_category), axis=1)

g = s7.merge(prev, left_on="gene", right_index=True, how="inner")
N_TABLE, N_MATCHED = s7.gene.nunique(), len(g)      # 99 named families (plus a note row), 59 matched
g["strict_gap"] = g["prev_FSM"] - g[["prev_RGM", "prev_SGM"]].max(axis=1)
g = g[g.strict_gap >= STRICT_GAP].copy()

# Three most FSM-specific genes per functional category.
g = g.sort_values("strict_gap", ascending=False)
sel = g.groupby("COG_category", sort=False).head(TOP_PER_CAT).copy()
cat_rank = sel.groupby("COG_category").strict_gap.max().sort_values(ascending=False)
sel["cat_rank"] = sel.COG_category.map({c: i for i, c in enumerate(cat_rank.index)})
sel = sel.sort_values(["cat_rank", "strict_gap"], ascending=[True, False])


def short_desc(gene):
    if gene in LABELS:
        return LABELS[gene]
    d = desc.get(gene)
    if not isinstance(d, str):
        return ""
    d = d.strip().rstrip(".")
    return d if len(d) <= 64 else d[:61].rstrip() + "\u2026"   # raw eggNOG text (no verified label)


fig = plt.figure(figsize=(WIDTH_FULL, 5.05))
gs = fig.add_gridspec(1, 1, left=0.052, right=0.985, top=0.905, bottom=0.175)

# ------------------------------------------------------------------ b
ax = fig.add_subplot(gs[0])
M = sel[["prev_SGM", "prev_RGM", "prev_FSM"]].values
n = len(sel)
im = ax.imshow(M, aspect="auto", cmap=seq_cmap("GnBu"), vmin=0, vmax=100)
ax.set_xticks([0, 1, 2])
ax.set_xticklabels(["SGM", "RGM", "FSM"], fontsize=8, fontweight="bold")
ax.xaxis.set_ticks_position("top")
ax.set_yticks(range(n))
ax.set_yticklabels([f"{r.gene}" for r in sel.itertuples()], fontsize=7,
                   fontstyle="italic")
for i in range(n):
    for j in range(3):
        v = M[i, j]
        ax.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=6.4,
                color="white" if v > 55 else INK)
for i, r in enumerate(sel.itertuples()):
    ax.text(5.05, i, short_desc(r.gene), va="center", fontsize=6.4, color=INK)
# Category brackets.
start = 0
for cat, grp in sel.groupby("cat_rank", sort=True):
    k = len(grp)
    ax.plot([2.60, 2.60], [start - 0.42, start + k - 0.58], color=INK, lw=0.9,
            clip_on=False)
    ax.text(2.70, start + (k - 1) / 2, grp.COG_category.iloc[0], va="center",
            fontsize=6.6, fontweight="bold", color=INK)
    start += k
ax.set_xlim(-0.5, 10.2)
ax.set_ylim(n - 0.5, -0.5)
for s in ax.spines.values():
    s.set_visible(False)
ax.tick_params(length=0)
cb = fig.colorbar(im, ax=ax, fraction=0.022, pad=0.045, location="bottom",
                  shrink=0.26, aspect=26, anchor=(0.02, 1.0))
cb.set_label("% of genomes encoding the family", fontsize=7)
cb.ax.tick_params(labelsize=6.6)
cb.outline.set_linewidth(0.5)
ax.set_title("Candidate FSM-associated gene families by functional category", loc="left",
             fontsize=8, pad=16)
fig.text(0.052, 0.012,
         f"A selected subset, not a validated functional signature. Of the 99 candidate families in Table S13, {N_MATCHED} map to a discriminating MMseqs2 cluster\n"
         f"(families file); {len(g)} of these {N_MATCHED} have an FSM prevalence at least 30 points above the higher of RGM and SGM, and the three largest gaps per COG\n"
         f"category are drawn ({n} families). Prevalence is over 134 FSM, 389 RGM and 299 SGM genomes (Table S13 uses 96 SPIRE FSM MAGs and pools RGM + SGM).\n"
         "Labels are short descriptions checked against Pfam / KEGG-KO evidence of the cluster representative.",
         fontsize=6.0, va="bottom", color=INK, linespacing=1.45)

save(fig, "FigureS3_gene_signature", supp=True)
print(f"genes passing strict gap >= {STRICT_GAP}: {len(g)} of {len(s7)} significant")
print("shown:", n, "genes in", sel.COG_category.nunique(), "categories")
print(f"Table S13 families {N_TABLE}; matched to S19 clusters {N_MATCHED}; qualifying {len(g)}")
print("qualifying but not drawn:", sorted(set(g.gene) - set(sel.gene)))
