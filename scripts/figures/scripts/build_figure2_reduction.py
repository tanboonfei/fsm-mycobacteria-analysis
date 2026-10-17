#!/usr/bin/env python3
"""Figure 2: genome reduction in FSM is non-degenerative, and differs from the
sporadic small genomes found elsewhere in the ingroup.

  a  genome architecture: size, coding density, GC           (Table S1; leprosy clade Table S4b)
  b  orthologue-guided core-gene decay screen                (Table S5i)
  c  size against coding density and GC                      (Tables S1, S4b)

Inputs : inputs/S1_genome_metadata.tsv, inputs/tsv/S4b_decayer_prodigal.tsv,
         inputs/tsv/S5i_per_genome_decay_screen_grouped.tsv
Output : figures/Figure2_reduction.{png,pdf}
Usage  : python3 build_figure2_reduction.py

Findings live in the manuscript legend. Panel titles are deliberately mixed:
some noun phrases, some none, so the set does not read in one register.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from figstyle import (CLADE_C, FSM_C, INK, LARGE_C, LEPRAE_C, ORDER, RULE,
                      SMALL_C, TSV, WIDTH_FULL, load_s1, panel_tag, save,
                      violin_box)

INGROUP_ANC, FSM_ANC, LEPRAE_PSEUDO = 5.1, 4.8, 1133

s1 = load_s1(pd)
# Coding density is the single uniform Prodigal call, applied to the whole frame
# so that panel c scatters the same quantity panel a summarises. The deposited
# column is superseded (Table S2) and must not reach either panel; it is left
# NaN for the two ingroup assemblies not re-called and for the outgroup; the
# coding-density panels carry 822 genomes where the others carry 824.
s1["coding_density_pct"] = pd.to_numeric(
    s1["coding_density_pct_uniform_prodigal"], errors="coerce")
ing = s1[s1.clade3.isin(ORDER)].copy()
prod = pd.read_csv(TSV / "S4b_decayer_prodigal.tsv", sep="\t")

fig = plt.figure(figsize=(WIDTH_FULL, 7.1))
gs = fig.add_gridspec(3, 1, height_ratios=[1.00, 0.76, 1.04], hspace=0.80,
                      left=0.085, right=0.985, top=0.905, bottom=0.120)

# ---------------- a: genome architecture ----------------------------------
gs_a = gs[0].subgridspec(1, 3, wspace=0.46, width_ratios=[1.32, 1.0, 1.0])
metrics = [("genome_size_Mb", "Genome size (Mb)", "{:.2f}"),
           ("coding_density_pct", "Coding density (%)", "{:.1f}"),
           ("GC_pct", "GC content (%)", "{:.1f}")]
# The host-restricted decay clade is carried as a fourth group so that the
# three architecture metrics can be read against a known degenerative
# reduction. Its coding density is the Prodigal call (Table S4b), matching
# the uniform Prodigal column used for the MAGs.
DECAY_KEY = {"genome_size_Mb": "genome_size_Mb",
             "coding_density_pct": "coding_density_prodigal",
             "GC_pct": "GC_pct"}
GROUPS_A = ORDER + ["Leprosy"]
LAB_A = {"SGM": "SGM", "RGM": "RGM", "FSM": "FSM", "Leprosy": "Leprosy\n(SGM subset)"}
COLOURS_A = {**CLADE_C, "Leprosy": LEPRAE_C}

axes_a = []
for j, (key, ylab, fmt) in enumerate(metrics):
    ax = fig.add_subplot(gs_a[j])
    vals = {c: ing.loc[ing.clade3 == c, key].dropna().values for c in ORDER}
    vals["Leprosy"] = pd.to_numeric(prod[DECAY_KEY[key]], errors="coerce").dropna().values
    violin_box(ax, vals, ylab, groups=GROUPS_A, colours=COLOURS_A, seed=j)
    # The leprosy clade is part of SGM: its three genomes are highlighted on the
    # SGM violin as well as drawn separately, and the separate column is labelled
    # as a subset rather than as a fourth clade.
    ax.scatter([0] * len(vals["Leprosy"]), vals["Leprosy"], s=13, marker="D",
               color=LEPRAE_C, edgecolors=INK, linewidths=0.4, zorder=6)

    # Medians (and group n) as a row above the axes: no annotation can then collide with a
    # violin, a box or the ancestral-state marks.
    for i, c in enumerate(GROUPS_A):
        ax.text(i, 1.015, f"{fmt.format(np.median(vals[c]))}\nn={len(vals[c])}",
                fontsize=6.4, color=INK, va="bottom", ha="center", linespacing=1.15,
                transform=ax.get_xaxis_transform())
    ax.set_xticks(range(len(GROUPS_A)))
    ax.set_xticklabels([LAB_A[g] for g in GROUPS_A], fontsize=7)
    ax.set_xlim(-0.60, 5.05 if j == 0 else 3.60)
    lo_y, hi_y = ax.get_ylim()            # keep the lowest leprosy points off the axis
    ax.set_ylim(lo_y - 0.04 * (hi_y - lo_y), hi_y)
    axes_a.append((ax, vals))
    if j == 0:
        sizes = vals
    if j == 2:
        gc = vals

ax0 = axes_a[0][0]
ax0.axhline(INGROUP_ANC, ls="--", lw=0.8, color=RULE, zorder=1)
ax0.text(5.00, INGROUP_ANC + 0.16, f"ingroup ancestor {INGROUP_ANC} Mb", fontsize=6.4,
         color=INK, va="bottom", ha="right")
ax0.annotate("", xy=(2, float(np.median(sizes["FSM"]))), xytext=(2, FSM_ANC),
             arrowprops=dict(arrowstyle="-|>", color=FSM_C, lw=1.2, shrinkA=0,
                             shrinkB=2))
ax0.plot([1.55, 2.45], [FSM_ANC, FSM_ANC], color=FSM_C, lw=1.2, ls=(0, (3, 2)),
         zorder=8)
ax0.text(2.42, FSM_ANC - 0.18, f"node before FSM stem\n{FSM_ANC} Mb (95% CI 2.1–7.4)", fontsize=6.4,
         color=INK, va="top", ha="left")

# No row subtitle: the three y-axis labels name the metrics and the manuscript
# legend carries the finding.
panel_tag(ax0, "a", dx=-0.42, dy=1.10)

# ---------------- b: orthologue-guided core-gene decay screen --------------
# Replaces the Bakta pseudogene boxplots. One stacked bar per genome: the 890
# orthologue-guided core groups split into intact, disrupted (frameshift or
# in-frame stop), missing (no hit) and edge-excluded (hit running off a contig
# end, not assessable). Table S5i; per-genome values (inputs/tsv/S5i_...) are
# the per_genome_decay_screen_grouped.tsv written by the decay-screen analysis.
# disrupted % printed here = pct_fs_or_stop x assessed / 890: a share of ALL
# 890 groups, not of assessable loci (the Results quote the latter).
dg = pd.read_csv(TSV / "S5i_per_genome_decay_screen_grouped.tsv", sep="\t")
GB = ["FSM SPIRE MAG (>=95%)", "RGM SPIRE MAG (>=95%), matched",
      "Closed RGM, intact", "Closed RGM, shredded",
      "Leprosy clade, intact", "Leprosy clade, shredded"]
GB_LAB = ["FSM MAGs\n\u2265 95%", "matched RGM\nMAGs", "closed RGM\nintact",
          "closed RGM\nshredded", "leprosy\nintact", "leprosy\nshredded"]
dg = dg[dg.group.isin(GB)].copy()
N = float(dg.core_groups.iloc[0])
dg["f_edge"] = dg.hits_at_contig_end / N
dg["f_missing"] = (N - dg.hits) / N
dg["f_disrupted"] = dg.pct_fs_or_stop / 100.0 * dg.assessed / N
dg["f_intact"] = 1 - dg.f_edge - dg.f_missing - dg.f_disrupted
dg["gi"] = dg.group.map({g: i for i, g in enumerate(GB)})
dg = dg.sort_values(["gi", "f_intact"], ascending=[True, False]).reset_index(drop=True)

ax = fig.add_subplot(gs[1])
CAT = [("f_intact", "intact", "#BFD8C6"),
       ("f_disrupted", "disrupted (frameshift or in-frame stop)", "#C0392B"),
       ("f_missing", "missing (no hit)", "#7F8FA6"),
       ("f_edge", "edge-excluded (runs off a contig end)", "#E3E3E3")]
# Each group occupies the same horizontal span (bar width = span / group size),
# so the three-genome leprosy groups are as legible as the 57-genome MAG
# panels; bar width therefore carries no information and the group size is
# printed under each block.
SPAN, GAP = 10.0, 2.6
xs, ws = [], []
for i, g in enumerate(GB):
    k = int((dg.gi == i).sum())
    w = SPAN / k
    x0 = i * (SPAN + GAP)
    xs += [x0 + (j + 0.5) * w for j in range(k)]
    ws += [w] * k
x = np.array(xs); w = np.array(ws)
bottom = np.zeros(len(dg))
bhandles = []
for key, lab, col in CAT:
    v = dg[key].values * 100
    h = ax.bar(x, v, bottom=bottom, width=w, color=col, linewidth=0, label=lab)
    bottom += v
    bhandles.append(h)
ax.set_ylim(0, 100)
ax.set_xlim(-1.0, (len(GB) - 1) * (SPAN + GAP) + SPAN + 1.0)
ax.set_yticks([0, 25, 50, 75, 100])
ax.set_ylabel("Core groups (%)")
cent = [i * (SPAN + GAP) + SPAN / 2 for i in range(len(GB))]
ax.set_xticks(cent)
ax.set_xticklabels([f"{GB_LAB[i]}\n(n = {int((dg.gi == i).sum())})" for i in range(len(GB))],
                   fontsize=6.6)
for i in range(len(GB)):
    med = dg.loc[dg.gi == i, "f_disrupted"].median() * 100
    ax.text(cent[i], 102, f"disrupted {med:.1f}%", ha="center", va="bottom",
            fontsize=6.4, color=INK)
ax.spines[["top", "right"]].set_visible(False)
ax.legend(handles=[h[0] for h in bhandles], labels=[c[1] for c in CAT],
          loc="upper center", bbox_to_anchor=(0.5, -0.46), ncol=4, frameon=False,
          fontsize=6.3, handlelength=1.2, handletextpad=0.45, columnspacing=1.3)
panel_tag(ax, "b", dx=-0.075, dy=1.22)

# ---------------- c: two routes to the same genome size --------------------
DECAY = "leprae|lepromatosis|uberis"
s1["mode"] = np.where(s1.species.astype(str).str.contains(DECAY, case=False, na=False), "decay",
             np.where(s1.clade3 == "FSM", "FSM",
             np.where(s1.clade3.isin(["SGM", "RGM"]) & (s1.genome_size_Mb < 4), "small",
             np.where(s1.clade3.isin(["SGM", "RGM"]), "large", "out"))))
MODE_C = {"FSM": FSM_C, "decay": LEPRAE_C, "small": SMALL_C, "large": LARGE_C}
MODE_L = {"large": "large-genome relatives (≥ 4 Mb)",
          "small": "sporadic small non-FSM (< 4 Mb)",
          "FSM": "FSM clade"}

gs_c = gs[2].subgridspec(1, 2, wspace=0.30)
handles = []
for k, (col, ylab) in enumerate([("coding_density_pct", "Coding density (%)"),
                                 ("GC_pct", "GC content (%)")]):
    ax = fig.add_subplot(gs_c[k])
    for md in ["large", "small", "FSM"]:
        sub = s1[s1["mode"] == md]
        h = ax.scatter(sub.genome_size_Mb, sub[col], s=5.0, color=MODE_C[md],
                       alpha=0.60, marker="o", edgecolors="none", zorder=2,
                       label=MODE_L[md] if k == 0 else None)
        if k == 0:
            handles.append(h)

    # Method-consistent (Prodigal) values for the decayers. GC does not depend
    # on the annotation method. Only the coding-density panel carries the
    # paired PGAP marker; the gap between the two reflects, among other
    # differences between the callers, the pseudogenised fraction of the genome.
    pcol = "coding_density_prodigal" if col == "coding_density_pct" else "GC_pct"
    if k == 0:
        # Both decayer series come from Table S4b: filled = the uniform Prodigal
        # call, hollow = the deposited PGAP annotation of the same genome. The
        # PGAP value is read explicitly and never from the frame behind the
        # scatter, which carries the uniform call for every genome.
        pr = prod.sort_values("genome_size_Mb")
        for _, b1 in pr.iterrows():
            ax.plot([b1.genome_size_Mb] * 2,
                    [b1["coding_density_ncbi_pgap"], b1[pcol]], color=LEPRAE_C,
                    lw=0.7, ls=":", zorder=4)
        h1 = ax.scatter(pr.genome_size_Mb, pr["coding_density_ncbi_pgap"], s=16,
                        marker="D", facecolors="none", edgecolors=LEPRAE_C,
                        linewidths=0.8, zorder=5,
                        label="the same genomes, PGAP annotation")
    h2 = ax.scatter(prod.genome_size_Mb, prod[pcol], s=16, marker="D",
                    color=LEPRAE_C, edgecolors=INK, linewidths=0.5, zorder=6,
                    label="leprosy clade (Prodigal)" if k == 0 else None)
    if k == 0:
        handles += [h2, h1]

    ax.axvspan(2.9, 3.7, color="#000000", alpha=0.05, zorder=0)
    ax.set_xlabel("Genome size (Mb)")
    ax.set_ylabel(ylab)
    ax.set_xlim(1.8, 9.6)
    ax.spines[["top", "right"]].set_visible(False)
    if k == 0:
        panel_tag(ax, "c", "Genomic characteristics at similar genome size", dx=-0.185,
                  dy=1.06)

# Legend below the two subplots, outside both axes, which means nothing can collide
# with the decayer markers at the lower left.
fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
           fontsize=6.8, handletextpad=0.5, columnspacing=1.4,
           bbox_to_anchor=(0.53, 0.005), markerscale=1.6)

save(fig, "Figure2_reduction")
print("size medians  :", {c: round(float(np.median(sizes[c])), 2) for c in ORDER})
print("GC medians    :", {c: round(float(np.median(gc[c])), 1) for c in ORDER})
print("decay screen medians (% of 890 core groups):")
print(dg.groupby("group")[["f_intact","f_disrupted","f_missing","f_edge"]].median().mul(100).round(2).to_string())
print("matched-size medians:",
      s1[s1.genome_size_Mb < 4].groupby("mode")[
          ["genome_size_Mb", "coding_density_pct", "GC_pct"]].median().round(2).to_dict())
print("decayer Prodigal coding density median:", round(float(prod.coding_density_prodigal.median()), 2))
