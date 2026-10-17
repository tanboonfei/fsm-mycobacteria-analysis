#!/usr/bin/env python3
"""Figure 3: selection under the two reduction modes.

  a  RELAX dN/dS by site class, FSM and the leprosy clade     (Table S8r, best multi-start fits)
  b  codon usage: genome-wide Nc, and the bulk-minus-ribosomal gap  (Table S9, S9d)
  c  piN/piS among conspecific MAG consensus genomes          (Table S9e)

Inputs : inputs/relax_best_fits.tsv   (made by extract_relax_best_fits.py from the
                                       multi-start HyPhy runs; the BEST fits of Table S8r)
         inputs/tsv/S9_codon_usage.tsv, inputs/tsv/S9e_within_pop_pNpS.tsv
Output : figures/Figure3_selection.{png,pdf}
Usage  : python3 build_figure3_selection.py

Panel a draws the rate classes of the same fits whose K and LRT it prints (the
best null and alternative over all multi-start fits, Table S8r), so panel and
numbers agree. The published single-start fits (Table S8) stopped ~310 lnL
units short of this optimum and are not drawn.

Statistics stay in the manuscript legend; panels carry noun phrases.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from figstyle import (FSM_C, INK, INPUTS, LEPRAE_C, RGM_C, RULE, TSV,
                      WIDTH_FULL, panel_tag, save, violin_box)

s14 = pd.read_csv(TSV / "S9_codon_usage.tsv", sep="\t")
s15 = pd.read_csv(TSV / "S9e_within_pop_pNpS.tsv", sep="\t")
relax = pd.read_csv(INPUTS / "relax_best_fits.tsv", sep="\t")

PAIR = ["RGM", "FSM"]
PAIR_C = {"RGM": RGM_C, "FSM": FSM_C}
# Panel b carries the host-restricted decay clade as a third group; panel c
# cannot, because within-population pi needs several genomes per population and
# the decay clade contributes one genome per species.
TRIO = PAIR + ["Decay"]
TRIO_C = {**PAIR_C, "Decay": LEPRAE_C}

fig = plt.figure(figsize=(WIDTH_FULL, 5.35))
gs = fig.add_gridspec(2, 1, height_ratios=[1.00, 0.86], hspace=0.52,
                      left=0.078, right=0.985, top=0.895, bottom=0.095)

# ------------------------------------------------------------------ a
# Corrected design: RGM+SGM are the only reference branches and the other
# reduced lineage is unclassified, and therefore it enters neither set.
# Text colours are darker than the line colours (gold on white is low contrast).
MODES = [("FSM", "FSM clade", FSM_C, "#8A5A00"),
         ("LEP", "leprosy clade", LEPRAE_C, LEPRAE_C)]
gs_a = gs[0].subgridspec(1, 2, wspace=0.14)
for k, (tag, name, col, tcol) in enumerate(MODES):
    ax = fig.add_subplot(gs_a[k])
    d = relax[relax.test == tag].sort_values("cls")
    K, lrt = float(d.K.iloc[0]), float(d.LRT.iloc[0])
    for _, row in d.iterrows():
        w0, w1, prop = row.omega_ref, row.omega_test, row.proportion
        # line width is affine in the codon fraction (0.8 + 4.0 x fraction pt)
        ax.plot([0, 1], [w0, w1], color=col, lw=0.8 + 4.0 * prop, alpha=0.85,
                solid_capstyle="round", zorder=3)
        ax.scatter([0, 1], [w0, w1], s=8, color=col, zorder=4)
        ax.text(-0.06, w0, f"{w0:.3g}", ha="right", va="center", fontsize=6.6,
                color=INK)
        ax.text(1.06, w1, f"{w1:.3g}", ha="left", va="center", fontsize=6.6,
                color=tcol, fontweight="bold" if prop > 0.5 else "normal")
        if prop > 0.5:
            ax.text(0.5, (w0 * w1) ** 0.5 * 1.6, f"{prop*100:.0f}% of codons",
                    ha="center", fontsize=6.6, color=INK)
    # neutral line stops short of the right-hand value labels; its name sits
    # above the line in the empty middle of the panel
    ax.hlines(1.0, -0.01, 1.0, ls="--", lw=0.7, color=RULE, zorder=1)
    ax.text(0.5, 1.28, "neutral (\u03c9 = 1)", fontsize=6.4, color=INK, va="bottom",
            ha="center",
            bbox=dict(facecolor="white", edgecolor="none", pad=0.6))
    ax.set_yscale("log")
    ax.set_ylim(4e-3, 200)
    ax.set_xlim(-0.44, 1.46)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["reference\n(RGM + SGM)", "test\n(terminals + stem)"])
    ax.spines[["top", "right"]].set_visible(False)
    if k == 0:
        ax.set_ylabel("d$_N$/d$_S$ (\u03c9), log scale")
    else:
        ax.set_yticklabels([])
    # K and the LRT of the fit drawn (best multi-start fit, Table S8r) go on the
    # title row, away from the class labels; verdicts and p-values are in the legend.
    ax.set_title(name, loc="left", fontsize=8, color=tcol, pad=4)
    ax.set_title(f"K = {K:.2f}, LRT = {lrt:.1f}" if tag == "FSM" else f"K = {K:.2f}, LRT = {lrt:.0f}",
                 loc="right", fontsize=6.8, color=tcol, pad=4)

panel_tag(fig.axes[0], "a", dx=-0.155, dy=1.14)

# ------------------------------------------------------------------ b and c
gs_b = gs[1].subgridspec(1, 3, wspace=0.52, width_ratios=[1.0, 1.0, 0.86])

ax = fig.add_subplot(gs_b[0])
nc = {c: pd.to_numeric(s14.loc[s14.clade == c, "bulk_Nc"], errors="coerce").dropna().values
      for c in TRIO}
violin_box(ax, nc, "Wright's Nc\n(higher = weaker bias)", groups=TRIO, colours=TRIO_C)
ax.set_xticklabels([f"{'Leprosy' if c == 'Decay' else c}\n(n = {len(nc[c])})\nmed {np.median(nc[c]):.1f}"
                    for c in TRIO], fontsize=7)
ax.set_xlim(-0.62, 2.62)
panel_tag(ax, "b", "Genome-wide codon bias", dx=-0.30, dy=1.10)

ax = fig.add_subplot(gs_b[1])
gap = {c: pd.to_numeric(s14.loc[s14.clade == c, "gap_bulk_minus_ribo"],
                        errors="coerce").dropna().values for c in TRIO}
violin_box(ax, gap, "Bulk − ribosomal Nc\n(codons)", groups=TRIO, colours=TRIO_C, seed=3)
ax.set_xticklabels([f"{'Leprosy' if c == 'Decay' else c}\n(n = {len(gap[c])})\nmed {np.median(gap[c]):.1f}"
                    for c in TRIO], fontsize=7)
ax.set_xlim(-0.62, 2.62)
ax.set_title("Translational-selection control", loc="left", fontsize=8,
             color=INK, pad=4)

# c: one point per population, x-range matched to the data so 11 points fill
# the panel instead of floating in an empty axis.
ax = fig.add_subplot(gs_b[2])
rng = np.random.default_rng(7)
med = {}
for i2, c in enumerate(PAIR):
    v = s15.loc[s15.clade == c, "piN_piS"].astype(float).values
    med[c] = float(np.median(v))
    ax.scatter(i2 + rng.normal(0, 0.055, len(v)), v, s=13, color=PAIR_C[c],
               alpha=0.9, lw=0, zorder=3)
    ax.plot([i2 - 0.26, i2 + 0.26], [med[c]] * 2, color=INK, lw=1.3, zorder=4)
ax.set_xticks(range(len(PAIR)))
# medians on their own line under each group; the previous single-line form
# ran the two labels together at print width
ax.set_xticklabels([f"{c}\n({int((s15.clade == c).sum())} species)\nmed {med[c]:.2f}"
                    for c in PAIR], fontsize=6.8)
ax.tick_params(axis="x", pad=2)
ax.set_xlim(-0.55, 1.55)
ax.set_ylim(0.03, 0.19)
ax.set_ylabel("$\\pi_N$/$\\pi_S$ among conspecific\nMAG consensus genomes")
ax.spines[["top", "right"]].set_visible(False)
ax.set_title("Species with \u2265 4 conspecific MAGs", loc="left", fontsize=8,
             color=INK, pad=4)
panel_tag(ax, "c", dx=-0.42, dy=1.10)           # y-axis already names it

save(fig, "Figure3_selection")
print("Nc medians :", {c: round(float(np.median(nc[c])), 1) for c in PAIR})
print("gap medians:", {c: round(float(np.median(gap[c])), 1) for c in PAIR})
print("piN/piS    :", {c: round(med[c], 3) for c in PAIR})
print("populations:", {c: int((s15.clade == c).sum()) for c in PAIR})
print("piN/piS range:", round(float(s15.piN_piS.min()), 3), "-",
      round(float(s15.piN_piS.max()), 3))
