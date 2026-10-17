#!/usr/bin/env python3
"""Figure 4: catabolic and osmotic modules rare or undetected in FSM, and habitat.

  a  prevalence of five catabolic / osmotic modules in RGM, FSM and the leprosy
     clade                                                      (Table S7)
  b  habitat of clade-dominant metagenome grid cells (SingleM/Sandpiper)  (Table S11d)
  c  detection in public metagenome runs by habitat, FSM vs RGM (Table S11i)
  d  aquatic against other biomes, from genome metadata          (Table S1)

Inputs : inputs/fig4a_values.tsv   panel a values (columns: module, gene_label,
                                   group, n_genomes, n_positive, pct, rule);
                                   if absent, the older inputs/tsv/S7_functional_modules.tsv
                                   (module, group, pct, order) is used, with a warning
         inputs/tsv/S11d_mag_global_grid_both.tsv, inputs/tsv/S11i_sandpiper_...tsv,
         inputs/S1_genome_metadata.tsv
Output : figures/Figure4_niche.{png,pdf}   (set MYCO_FIG_OUT to write elsewhere)
Usage  : python3 build_figure4_niche.py

Fisher p-values and odds ratios live in the manuscript legend, not in the
panel titles, which overran the panel boundary and collided with each other.
"""
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import figstyle
from figstyle import (AQUATIC_C, FSM_C, HABITAT_C, INK, INPUTS, LEPRAE_C,
                      ORDER, OTHER_BIOME_C, RGM_C, TSV, WIDTH_FULL, load_s1,
                      panel_tag, save)

if os.environ.get("MYCO_FIG_OUT"):          # e.g. for a preview build
    figstyle.OUT = Path(os.environ["MYCO_FIG_OUT"])

# ---- panel a values: one TSV path variable, with a fallback to the S7 export
FIG4A_TSV = Path(os.environ.get("MYCO_FIG4A_TSV", INPUTS / "fig4a_values.tsv"))
FIG4A_FALLBACK = TSV / "S7_functional_modules.tsv"
if FIG4A_TSV.exists():
    s3 = pd.read_csv(FIG4A_TSV, sep="\t")
    print(f"panel a values: {FIG4A_TSV}")
else:
    print(f"!! {FIG4A_TSV} not found; panel a falls back to {FIG4A_FALLBACK.name} "
          "(older Bakta-product values, no numerators)", file=sys.stderr)
    s3 = pd.read_csv(FIG4A_FALLBACK, sep="\t")
    s3["gene_label"] = s3["module"]
    s3["n_genomes"] = np.nan
    s3["n_positive"] = np.nan
    s3["rule"] = "fallback (S7 export)"
s3["group"] = s3["group"].replace({"Leprosy clade": "Decay"})   # internal key
if "order" not in s3.columns:
    s3["order"] = s3.groupby("module", sort=False).ngroup()

HAB_ORDER = ["freshwater", "soil/sediment", "wastewater", "drinking-water",
             "marine/saline", "host-assoc.", "aquaculture", "other"]

grid = pd.read_csv(TSV / "S11d_mag_global_grid_both.tsv", sep="\t")
# The grid file predates the final clade names: RGM1_saline_large = RGM,
# RGM2_fresh_small = FSM. Two further cells are dominated by other lineages
# and are not plotted (542 cells in the file, 540 plotted).
grid["dom_clade"] = grid["dom_clade"].replace({"RGM1_saline_large": "RGM",
                                               "RGM2_fresh_small": "FSM"})
s1 = load_s1(pd)

fig = plt.figure(figsize=(WIDTH_FULL, 5.00))
gs = fig.add_gridspec(2, 3, height_ratios=[0.66, 1.0],
                      width_ratios=[1.0, 1.15, 0.80],
                      hspace=0.58, wspace=0.75, left=0.075, right=0.985,
                      top=0.905, bottom=0.175)

# ------------------------------------------------------------------ a
ax = fig.add_subplot(gs[0, :])
s3 = s3[~s3["module"].str.contains("hodopsin", case=False)]      # rhodopsin is not a Fig. 4a module
mods = list(dict.fromkeys(s3.sort_values("order")["module"]))
x = np.arange(len(mods))
# The host-restricted decayers are carried alongside FSM: both reduced
# lineages have lost these modules, and the two modes separate on the other
# signatures, not on this panel.
w = 0.27


def rows(group):
    d = s3[s3.group == group].set_index("module").loc[mods]
    return d


def label_of(m):
    g = rows("RGM").loc[m, "gene_label"] if "gene_label" in s3 else m
    g = m if (not isinstance(g, str) or g == m) else g
    return g.replace(" (", "\n(").replace("K+", "K$^{+}$")


grp = {g: rows(g) for g in ("RGM", "FSM", "Decay")}
n_of = {g: grp[g].n_genomes.max() for g in grp}
def _n(g):
    return f" (n = {int(n_of[g])})" if np.isfinite(n_of[g]) else ""
ax.bar(x - w, grp["RGM"].pct, w, color=RGM_C, label="RGM" + _n("RGM"), edgecolor="white", lw=0.5)
ax.bar(x, grp["FSM"].pct, w, color=FSM_C, label="FSM" + _n("FSM"), edgecolor="white", lw=0.5)
ax.bar(x + w, grp["Decay"].pct, w, color=LEPRAE_C,
       label="Leprosy clade, three tissue-derived genomes", edgecolor="white", lw=0.5)
for off, g in ((-w, "RGM"), (0, "FSM"), (w, "Decay")):
    for xi, v in zip(x + off, grp[g].pct):
        txt = "0%" if v == 0 else ("<1%" if v < 0.5 else f"{v:.0f}%")
        ax.text(xi, v + 2.5, txt, ha="center", fontsize=6.4, color=INK)
for k, m in enumerate(mods):
    if "ctC" in m or "kdp" in m:
        ax.annotate("osmotic", xy=(k, 95), ha="center", fontsize=6.6, color=INK,
                    style="italic")
ax.set_xticks(x)
ax.set_xticklabels([label_of(m) for m in mods])
ax.set_ylabel("% of genomes encoding")
ax.set_ylim(0, 108)
ax.legend(frameon=False, ncol=1, loc="upper center", fontsize=7.2,
          bbox_to_anchor=(0.43, 1.02), handlelength=1.1, labelspacing=0.3)
ax.spines[["top", "right"]].set_visible(False)
panel_tag(ax, "a", "Catabolic and osmotic modules rare or undetected in FSM", dx=-0.052, dy=1.10)

# ------------------------------------------------------------------ b
ax = fig.add_subplot(gs[1, 0])
comp = {c: grid[grid.dom_clade == c].dom_biome.value_counts(normalize=True) * 100
        for c in ORDER}
bottoms = {c: 0.0 for c in ORDER}
for hab in HAB_ORDER:
    vals = [float(comp[c].get(hab, 0.0)) for c in ORDER]
    ax.bar(range(len(ORDER)), vals, bottom=[bottoms[c] for c in ORDER],
           color=HABITAT_C[hab], width=0.62, edgecolor="white", lw=0.6,
           label=hab if max(vals) > 0 else None)
    for i, c in enumerate(ORDER):
        if vals[i] >= 10:
            ax.text(i, bottoms[c] + vals[i] / 2, f"{vals[i]:.0f}%", ha="center",
                    va="center", fontsize=6.4, color="white", fontweight="bold")
        bottoms[c] += vals[i]
ax.set_xticks(range(len(ORDER)))
ax.set_xticklabels([f"{c}\n{(grid.dom_clade == c).sum()} cells" for c in ORDER], fontsize=6.6)
ax.set_ylabel("Clade-dominant grid cells (%)")
ax.set_ylim(0, 100)
ax.spines[["top", "right"]].set_visible(False)
ax.legend(frameon=False, fontsize=6.4, ncol=2, loc="upper center",
          bbox_to_anchor=(0.5, -0.20), handlelength=1.0, columnspacing=0.8,
          handletextpad=0.5, labelspacing=0.35)
panel_tag(ax, "b", "Habitat of clade-dominant\nmetagenomes", dx=-0.20, dy=1.06)

# ------------------------------------------------------------------ c
# Read-based detection in public metagenome runs (Sandpiper, GTDB R232), the
# replacement for the MAG salinity panel, which did not survive control for
# source study (Table S12g).
ax = fig.add_subplot(gs[1, 1])
sp = pd.read_csv(TSV / "S11i_sandpiper_prevalence_abundance_by_habitat.tsv", sep="\t").set_index("habitat")
# Nine of the 13 habitat classes of Table S11i are drawn. Omitted: saline/hypersaline lake,
# freshwater (unspecified), aquatic (unspecified), other/unspecified. The FSM-RGM
# comparison for those is in Table S11i. RGM here includes the sister RGM lineages
# (295 species against 42 FSM species); panels a, b and d use RGM only.
SP_ROWS = [("lake/reservoir", "Lake/reservoir"), ("river/stream", "River/stream"),
           ("estuary/brackish", "Estuary/brackish"), ("marine", "Marine"),
           ("groundwater/drinking water", "Ground/drinking water"), ("sediment", "Sediment"),
           ("soil/peat/wetland", "Soil/peat/wetland"), ("wastewater/engineered", "Wastewater"),
           ("host-associated", "Host-associated")]
y = np.arange(len(SP_ROWS))[::-1]
hh = 0.38
fv = [float(sp.loc[k, "FSM_prev_pct"]) for k, _ in SP_ROWS]
rv = [float(sp.loc[k, "RGM_prev_pct"]) for k, _ in SP_ROWS]
ax.barh(y + hh / 2, fv, hh, color=FSM_C, label="FSM", edgecolor="white", lw=0.4)
ax.barh(y - hh / 2, rv, hh, color=RGM_C, label="RGM (+ sister RGM)", edgecolor="white", lw=0.4)
for yy, v in zip(y + hh / 2, fv):          # label bars too short to see
    if v < 1:
        ax.text(v + 0.3, yy, f"{v:.2g}", va="center", fontsize=6.0, color=INK)
for yy, v in zip(y - hh / 2, rv):
    if v < 1:
        ax.text(v + 0.3, yy, f"{v:.2g}", va="center", fontsize=6.0, color=INK)
ax.set_yticks(y)
ax.set_yticklabels([lab for _, lab in SP_ROWS], fontsize=6.4)
ax.set_xlabel("Runs with the clade detected (%)")
ax.set_xlim(0, 32)
ax.spines[["top", "right"]].set_visible(False)
ax.legend(frameon=False, fontsize=6.2, loc="center right", handlelength=1.0,
          handletextpad=0.5, labelspacing=0.35, bbox_to_anchor=(1.02, 0.45))
ax.text(0.0, -0.30, "9 of 13 habitat classes shown (Table S11i)", transform=ax.transAxes,
        fontsize=6.2, color=INK, va="top")
panel_tag(ax, "c", "Detection in public\nmetagenome runs", dx=-0.62, dy=1.06)

# ------------------------------------------------------------------ d
ax = fig.add_subplot(gs[1, 2])
bio = s1[s1.biome.notna() & s1.clade3.isin(ORDER)].copy()
bio["aq"] = bio.biome == "aquatic"
AQ_C = [AQUATIC_C, OTHER_BIOME_C]
bottoms = {c: 0.0 for c in ORDER}
for j, (lab, mask) in enumerate([("aquatic", True), ("other biomes (terrestrial,\nanthropogenic, host)", False)]):
    vals, ns = [], []
    for c in ORDER:
        sub = bio[bio.clade3 == c]
        vals.append(100 * (sub.aq == mask).mean())
        ns.append((sub.aq == mask).sum())
    ax.bar(range(len(ORDER)), vals, bottom=[bottoms[c] for c in ORDER],
           color=AQ_C[j], width=0.60, edgecolor="white", lw=0.6, label=lab)
    for i, c in enumerate(ORDER):
        if vals[i] >= 12:
            ax.text(i, bottoms[c] + vals[i] / 2, f"{ns[i]}", ha="center",
                    va="center", fontsize=6.8,
                    color="white" if j == 0 else INK, fontweight="bold")
        bottoms[c] += vals[i]
ax.set_xticks(range(len(ORDER)))
ax.set_xticklabels([f"{c}\nn={(bio.clade3 == c).sum()}" for c in ORDER], fontsize=6.8)
ax.set_ylabel("Genomes with biome metadata (%)")
ax.set_ylim(0, 100)
ax.spines[["top", "right"]].set_visible(False)
ax.legend(frameon=False, fontsize=6.4, loc="upper center",
          bbox_to_anchor=(0.5, -0.20), handlelength=1.0, handletextpad=0.5,
          labelspacing=0.35)
panel_tag(ax, "d", "Aquatic against\nother biomes", dx=-0.30, dy=1.06)

save(fig, "Figure4_niche")
for c in ORDER:
    print(f"{c}: {(grid.dom_clade == c).sum():3d} cells, "
          f"freshwater {float(comp[c].get('freshwater', 0.0)):.0f}%")
print("aquatic counts:", {c: int(((bio.clade3 == c) & bio.aq).sum()) for c in ORDER},
      "of", {c: int((bio.clade3 == c).sum()) for c in ORDER})
print("sandpiper prevalence FSM/RGM:", list(zip([k for k, _ in SP_ROWS], fv, rv)))
