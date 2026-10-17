#!/usr/bin/env python3
"""Supplementary Figure S1: the reduction metrics survive a MAG-completeness control.

Pseudogene count, coding density and genome size against estimated completeness
for the 243 ingroup MAGs (FSM 96, RGM 97, SGM 50), with an analysis of
covariance (OLS on the raw values; completeness, contamination and clade, RGM as
reference) giving the FSM-minus-RGM term printed in each panel title.

Inputs : inputs/S1_genome_metadata.tsv           (Table S1: clade, completeness_used_pct,
                                                  checkm2_this_study_contamination_pct,
                                                  coding_density_pct_uniform_prodigal, genome_size_Mb)
         inputs/tsv/S5h_pseudogenes_MAG_panel.tsv (Table S5h: Bakta pseudogene count per MAG)
Output : supplementary_figures/FigureS1_completeness_control.{png,pdf}
Usage  : python3 build_figureS1_completeness.py

The p-values are uncorrected for phylogeny (genomes are treated as independent)
and show direction only; the phylogenetically corrected analysis in Table S3
governs the main-text claims. This is stated on the image itself.
"""
import matplotlib.pyplot as plt
import pandas as pd
import statsmodels.formula.api as smf

from figstyle import (CLADE_C, INK, ORDER, TSV, WIDTH_FULL, load_s1, panel_tag,
                      pfmt, save)

# (column, axis label, panel title, symlog y axis, effect-size format)
METRICS = [
    ("pseudogenes",        "Pseudogenes per genome (symlog)", "Pseudogene load", True,
     "{:+.1f} per genome"),
    ("coding_density_pct", "Coding density (%)",     "Coding density",  False,
     "{:+.2f} points"),
    ("genome_size_Mb",     "Genome size (Mb)",       "Genome size",     False,
     "{:+.2f} Mb"),
]

s1 = load_s1(pd)
s6 = pd.read_csv(TSV / "S5h_pseudogenes_MAG_panel.tsv", sep="\t")[["genome", "pseudogenes"]]
df = s1.merge(s6, left_on="genome_id_or_accession", right_on="genome", how="left")
mag = df[df.completeness_pct.notna() & df.clade3.isin(ORDER)].copy()
# Completeness and contamination are taken from this study's single CheckM2
# v1.0.1 run, not from the values distributed with SPIRE; the completeness
# control uses the same estimates as every other analysis in the paper.
mag["completeness_pct"] = pd.to_numeric(mag["completeness_used_pct"],
                                        errors="coerce")
mag["contamination_pct"] = pd.to_numeric(
    mag["checkm2_this_study_contamination_pct"], errors="coerce")
# coding density from the single uniform Prodigal call applied to every genome
mag["coding_density_pct"] = pd.to_numeric(
    mag["coding_density_pct_uniform_prodigal"], errors="coerce").fillna(
    pd.to_numeric(mag["coding_density_pct"], errors="coerce"))

fig = plt.figure(figsize=(WIDTH_FULL, 2.75))
gs = fig.add_gridspec(1, 3, wspace=0.36, left=0.068, right=0.988, top=0.725,
                      bottom=0.250)

stats_out = []
for j, (col, ylab, short, logy, efmt) in enumerate(METRICS):
    ax = fig.add_subplot(gs[j])
    sub = mag[mag[col].notna()]
    model = smf.ols(
        f"Q('{col}') ~ completeness_pct + contamination_pct + "
        "C(clade3, Treatment('RGM'))", data=sub).fit()
    coef_name = [k for k in model.params.index if "FSM" in k][0]
    p = model.pvalues[coef_name]
    beta = model.params[coef_name]
    stats_out.append((col, beta, p, len(sub)))

    for c in ORDER:
        s = sub[sub.clade3 == c]
        ax.scatter(s.completeness_pct, s[col], s=4.5, color=CLADE_C[c],
                   alpha=0.65, lw=0, label=c)
    if logy:
        ax.set_yscale("symlog", linthresh=1)
    ax.set_xlabel("Estimated completeness (%)")
    ax.set_ylabel(ylab)
    ax.spines[["top", "right"]].set_visible(False)
    # panel_tag would overwrite the title, and the letter is drawn without one
    # and the title is set afterwards. Second line is the FSM-against-RGM clade
    # term after controlling completeness and contamination: direction only.
    panel_tag(ax, "abc"[j], dx=-0.22, dy=1.30)
    ax.set_title(f"{short}\nFSM \u2212 RGM {efmt.format(beta).replace("-", chr(0x2212))}; {pfmt(p)}*", loc="left",
                 fontsize=7, color=INK, pad=4, linespacing=1.5)

handles = [plt.Line2D([], [], marker="o", ls="", ms=3.4, color=CLADE_C[c],
                      label=f"{c} (n = {int((mag.clade3 == c).sum())})")
           for c in ORDER]
fig.legend(handles=handles, loc="upper right", ncol=3, frameon=False, fontsize=7,
           bbox_to_anchor=(0.988, 1.005), handletextpad=0.3, columnspacing=1.2)

fig.text(0.068, 0.012,
         "* Uncorrected for phylogeny: OLS on raw values, genomes treated as independent; value ~ completeness + contamination + clade\n"
         "(reference RGM), 243 ingroup MAGs; effect = FSM minus RGM coefficient. Direction only; phylogenetically corrected tests are in Table S3.",
         fontsize=6.2, color=INK, va="bottom", ha="left", linespacing=1.3)

save(fig, "FigureS1_completeness_control", supp=True)
for col, beta, p, n in stats_out:
    print(f"{col:20s} n={n:3d}  FSM beta={beta:+9.3f}  clade p={p:.2e}  -> {pfmt(p)}")
