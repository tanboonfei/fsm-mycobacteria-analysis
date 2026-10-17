#!/usr/bin/env python3
"""Figure 1: bac120 phylogeny of the Mycobacteriaceae with per-genome traits.

Radial tree (862 genomes) with, from the tips outward: genome size, urease,
ectoine synthase and microbial-opsin presence rings, the three ingroup clade
arcs (SGM, RGM in two segments, FSM), the Tsukamurella rooting-taxa arc, and the
Gupta et al. genus band. Branch colour = tip provenance / culture category.

Inputs : inputs/myco_bac120.rooted.treefile, inputs/S1_genome_metadata.tsv
         (the first 862 rows of Table S1, Additional file 2)
Output : figures/Figure1_phylogeny_traits.{png,pdf}
Usage  : python3 build_figure1_phylogeny.py

Two different "ingroup" sets are used and kept apart in the legend:
  * the 843-genome non-Tsukamurella crown (SGM, RGM, FSM and the 19 basal
    Mycobacteroides-grade genomes) defines the radial scale break (D_CROWN);
  * the 824-genome crown (SGM + RGM + FSM) is what the clade arcs cover.
"""
from pathlib import Path

import matplotlib as mpl
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import Phylo
from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.patches import Rectangle, Wedge

from figstyle import (ABSENT_C, BRANCH_ALIAS, BRANCH_C, CLADE_C, FSM_C,
                      GENUS_BAND_C, GENUS_NAMES, INK, INPUTS, INTERNAL_C,
                      OUTSIDE_C, PRESENT_C, RGM_C, SGM_C, WIDTH_FULL, save)

PRESENT = PRESENT_C

GAP_DEG = 9.0          # wedge left open at the top
START_DEG = 90 + GAP_DEG / 2
SPAN_DEG = 360 - GAP_DEG

R_IN, R_OUT = 0.055, 1.00     # tree occupies this radial band
# Rings are deliberately thin: the branching structure is the message, and
# five wide rings previously took more radius than the tree itself.
RINGS = [                      # (key, kind, width, colormap/None, label)
    ("genome_size_Mb",      "cont", 0.060, "GnBu",    "Genome size (Mb)"),
    ("urease",              "pres", 0.034, None,      "Urease"),
    ("ectoine_synthase",    "pres", 0.034, None,      "Ectoine synthase"),
    ("rhodopsin",           "pres", 0.040, None,      "Microbial rhodopsin"),
]
RING_GAP = 0.009

# ------------------------------------------------------------------ load
meta = pd.read_csv(INPUTS / "S1_genome_metadata.tsv", sep="\t").set_index("tip")
meta["clade3"] = meta["clade"].str.split(" ").str[0]
# "Outside ingroup (basal Mycobacteriaceae...)" -> "Outside"
meta.loc[meta["clade"].str.startswith("Outside"), "clade3"] = "Outside"
# Genus under the 2018 division, taken from position on the tree rather than
# from the deposited names: 15 rapid growers keep Mycobacterium names and one
# slow grower carries a Mycolicibacterium name (Table S1).
GENUS_OF = {"SGM": "Mycobacterium s.s.", "RGM": "Mycolicibacterium",
            "FSM": "Mycolicibacterium", "Outside": "Mycobacteroides"}
meta["genus2018"] = meta["clade3"].map(GENUS_OF)

tree = Phylo.read(INPUTS / "myco_bac120.rooted.treefile", "newick")
tree.ladderize(reverse=True)
tips = [t.name for t in tree.get_terminals()]
n = len(tips)

# Radial position: root-to-node distance, compressed so the long outgroup
# branch does not squash the ingroup.
depths = tree.depths(unit_branch_lengths=False)
raw = np.array([depths[t] for t in tree.get_terminals()])
dmax = raw.max()

# Radial scale break. Root-to-crown is a long bare stem -- 0.61 of the 1.39
# substitutions per site spanned by the tree -- carrying nothing but the
# Tsukamurella outgroup. Drawn to scale it emptied the middle of the plot and
# crammed every divergence into the outer third of the radius. The stem is
# therefore compressed into the inner quarter and the crown of the 843
# non-Tsukamurella genomes is drawn linearly beyond it (74% of the radius
# instead of 34%).
_clade_of = meta["clade"].astype(str).to_dict()
_crown843 = [t for t in tree.get_terminals()
            if not _clade_of.get(t.name, "").startswith("Outgroup")]
D_CROWN = depths[tree.common_ancestor(_crown843)]
R_BREAK = 0.26


def radius(d):
    d = max(d, 0.0)
    if d <= D_CROWN:
        return R_IN + (R_BREAK - R_IN) * (d / D_CROWN)
    return R_BREAK + (R_OUT - R_BREAK) * (d - D_CROWN) / (dmax - D_CROWN)


angle_of = {t: START_DEG - (i + 0.5) / n * SPAN_DEG for i, t in enumerate(tips)}
node_angle, node_r = {}, {}
for cl in tree.find_clades(order="postorder"):
    node_r[cl] = radius(depths[cl])
    if cl.is_terminal():
        node_angle[cl] = angle_of[cl.name]
    else:
        node_angle[cl] = float(np.mean([node_angle[c] for c in cl.clades]))


def xy(r, deg):
    a = np.deg2rad(deg)
    return r * np.cos(a), r * np.sin(a)


def arc_pts(r, a0, a1, step=0.35):
    k = max(2, int(abs(a1 - a0) / step) + 1)
    a = np.linspace(a0, a1, k)
    return np.column_stack(xy(r, a))


# Branch colour: tips take their category, internal branches stay neutral
# unless every descendant is the outgroup.
def branch_colour(cl):
    if cl.is_terminal():
        row = meta.loc[cl.name] if cl.name in meta.index else None
        if row is None:
            return INTERNAL_C
        if row["clade3"] == "Outside":
            return OUTSIDE_C
        status = row["culture_status"]
        status = BRANCH_ALIAS.get(status, status)
        return BRANCH_C.get(status, INTERNAL_C)
    names = [t.name for t in cl.get_terminals()]
    sub = meta.reindex(names)
    if (sub["clade3"] == "Outgroup").all():
        return BRANCH_C["Outgroup"]
    if (sub["clade3"] == "Outside").all():
        return OUTSIDE_C
    return INTERNAL_C


FIG_H = 8.55
LEG_H = 0.250                  # legend band, as a fraction of figure height
fig = plt.figure(figsize=(WIDTH_FULL, FIG_H))
ax = fig.add_axes([0.0, LEG_H, 1.0, 1.0 - LEG_H])
ax.set_aspect("equal")
ax.axis("off")

segs, cols, lws = [], [], []
for cl in tree.find_clades(order="preorder"):
    c = branch_colour(cl)
    for ch in cl.clades:
        segs.append([xy(node_r[cl], node_angle[ch]), xy(node_r[ch], node_angle[ch])])
        cols.append(branch_colour(ch))
        lws.append(0.42 if ch.is_terminal() else 1.05)
    if len(cl.clades) > 1:
        a = sorted(node_angle[ch] for ch in cl.clades)
        pts = arc_pts(node_r[cl], a[0], a[-1])
        segs.extend([[pts[i], pts[i + 1]] for i in range(len(pts) - 1)])
        cols.extend([c] * (len(pts) - 1))
        lws.extend([1.05] * (len(pts) - 1))
ax.add_collection(LineCollection(segs, colors=cols, linewidths=lws, zorder=2))

# Ultrafast-bootstrap support on internal nodes.
sx, sy, sc = [], [], []
for cl in tree.get_nonterminals():
    v = cl.confidence
    if v is None:
        continue
    if v >= 95:
        col = "#111111"
    elif v >= 90:
        col = "#9A9A9A"
    else:
        continue
    x, y = xy(node_r[cl], node_angle[cl])
    sx.append(x), sy.append(y), sc.append(col)
ax.scatter(sx, sy, s=1.5, c=sc, lw=0, zorder=3)

# ------------------------------------------------------------------ rings
dtheta = SPAN_DEG / n
r0 = R_OUT + 0.026
ring_bounds = {}
for key, kind, width, cmap_name, _lab in RINGS:
    vals = meta.reindex(tips)[key]
    if kind == "cont":
        v = vals.astype(float)
        # clip at the 2nd and 98th percentiles (2.859 and 7.385 Mb), rounded to
        # 0.1 Mb so the colour-bar ends read as the clipping values
        lo = round(float(np.nanpercentile(v, 2)), 1)
        hi = round(float(np.nanpercentile(v, 98)), 1)
        norm = mpl.colors.Normalize(lo, hi)
        cmap = mpl.colors.LinearSegmentedColormap.from_list(
            cmap_name + "_trunc",
            plt.get_cmap(cmap_name)(np.linspace(0.18, 0.74, 256)))
        colours = [cmap(norm(x)) if np.isfinite(x) else "#FFFFFF" for x in v]
        ring_bounds[key] = (lo, hi, cmap_name)
    else:
        colours = [PRESENT[key] if s == "present" else ABSENT_C for s in vals]

    patches = []
    for i, t in enumerate(tips):
        a1 = START_DEG - i * dtheta
        patches.append(Wedge((0, 0), r0 + width, a1 - dtheta, a1, width=width))
    ax.add_collection(PatchCollection(patches, facecolors=colours,
                                      edgecolors="none", zorder=2))
    r0 += width + RING_GAP

# ------------------------------------------------------------------ clade arcs
r_arc = r0 + 0.018
labels_seq = [meta.at[t, "clade3"] for t in tips]

# Contiguous runs, which means the paraphyletic RGM grade is drawn as the two
# segments that flank the nested FSM clade instead of one arc spanning it.
runs = []
start = 0
for i in range(1, len(labels_seq) + 1):
    if i == len(labels_seq) or labels_seq[i] != labels_seq[start]:
        runs.append((labels_seq[start], start, i - 1))
        start = i

ARC_C = {**CLADE_C, "Outgroup": BRANCH_C["Outgroup"]}   # Tsukamurella arc in the rooting-taxa black
for cl_name, i0, i1 in runs:
    if cl_name not in ARC_C:
        continue
    a_hi = START_DEG - i0 * dtheta
    a_lo = START_DEG - (i1 + 1) * dtheta
    pts = arc_pts(r_arc, a_lo, a_hi, step=0.2)
    ax.plot(pts[:, 0], pts[:, 1], color=ARC_C[cl_name], lw=3.6,
            solid_capstyle="butt", zorder=4)

# Label each clade once, on its longest run.
for cl_name in ["SGM", "RGM", "FSM"]:
    cand = [r for r in runs if r[0] == cl_name]
    if not cand:
        continue
    _, i0, i1 = max(cand, key=lambda r: r[2] - r[1])
    amid = START_DEG - (i0 + i1 + 1) / 2 * dtheta
    lx, ly = xy(r_arc + 0.042, amid)
    # Tangential, then flipped back if it would read upside down. The old
    # rule keyed off the raw angle and inverted the FSM label.
    rot = (amid - 90) % 360
    if 90 < rot < 270:
        rot -= 180
    ax.text(lx, ly, cl_name, color=CLADE_C[cl_name], fontsize=10,
            fontweight="bold", ha="center", va="center", rotation=rot,
            rotation_mode="anchor", zorder=6,
            path_effects=[pe.withStroke(linewidth=2.6, foreground="white")])

# ------------------------------------------------- genus band (2018 division)
r_gen = r_arc + 0.072
LIM = r_gen + 0.32
IN_PER_UNIT = min(WIDTH_FULL, FIG_H * (1.0 - LEG_H)) / (2 * LIM)


def fits_tangentially(text, fontsize, r, span_deg):
    """True if the label is shorter than the arc of the run it belongs to."""
    width_in = len(text) * 0.56 * fontsize / 72.0
    arc_in = np.deg2rad(span_deg) * r * IN_PER_UNIT
    return arc_in >= 1.08 * width_in


def leader_label(text, amid, colour, fontsize, italic=False, dy=0.0):
    """Horizontal label on a two-segment leader, for runs too short to letter."""
    side = 1 if np.cos(np.deg2rad(amid)) >= 0 else -1
    x0, y0 = xy(r_gen + 0.030, amid)
    xe, ye = xy(r_gen + 0.255, amid)
    x_lab = side * max(abs(xe) + 0.03, 0.30)
    ax.plot([x0, xe, x_lab], [y0, ye, ye + dy], color=colour, lw=0.6,
            solid_joinstyle="round", zorder=5, alpha=0.8)
    ax.text(x_lab + side * 0.022, ye + dy, text, color=colour, fontsize=fontsize,
            fontweight="bold", style="italic" if italic else "normal",
            ha="left" if side > 0 else "right", va="center", zorder=6,
            path_effects=[pe.withStroke(linewidth=2.4, foreground="white")])
gen_seq = [meta.at[t, "genus2018"] for t in tips]
gruns, gstart = [], 0
for i in range(1, len(gen_seq) + 1):
    if i == len(gen_seq) or gen_seq[i] != gen_seq[gstart]:
        gruns.append((gen_seq[gstart], gstart, i - 1))
        gstart = i

for g_name, i0, i1 in gruns:
    if g_name not in GENUS_NAMES:
        continue
    a_hi = START_DEG - i0 * dtheta
    a_lo = START_DEG - (i1 + 1) * dtheta
    pts = arc_pts(r_gen, a_lo, a_hi, step=0.2)
    ax.plot(pts[:, 0], pts[:, 1], color=GENUS_BAND_C, lw=4.2,
            solid_capstyle="butt", zorder=4)

for g_name in GENUS_NAMES:
    cand = [r for r in gruns if r[0] == g_name]
    if not cand:
        continue
    _, i0, i1 = max(cand, key=lambda r: r[2] - r[1])
    amid = START_DEG - (i0 + i1 + 1) / 2 * dtheta
    span = (i1 - i0 + 1) * dtheta
    if not fits_tangentially(g_name, 9.0, r_gen + 0.046, span):
        # Mycobacteroides spans ~8 deg; set tangentially it ran across its
        # neighbours' bands and collided with the Mycobacterium s.s. label.
        leader_label(g_name, amid, GENUS_BAND_C, 8.5, italic=True)
        continue
    lx, ly = xy(r_gen + 0.046, amid)
    rot = (amid - 90) % 360
    if 90 < rot < 270:
        rot -= 180
    ax.text(lx, ly, g_name, color=GENUS_BAND_C, fontsize=9.0,
            fontweight="bold", style="italic", ha="center", va="center",
            rotation=rot, rotation_mode="anchor", zorder=6,
            path_effects=[pe.withStroke(linewidth=2.6, foreground="white")])

# The 19 Tsukamurella tips (rooting taxa) have a clade arc but no genus-band
# segment; label the arc through the same two-segment leader.
for cl_name, i0, i1 in runs:
    if cl_name == "Outgroup":
        amid = START_DEG - (i0 + i1 + 1) / 2 * dtheta
        leader_label("Tsukamurella (rooting taxa)", amid, BRANCH_C["Outgroup"], 8.0,
                     italic=True, dy=0.07)

# Orientation marks for well-known taxa. Each label sits on the middle tip of
# all tips carrying that species name, so it marks the species clade and not
# necessarily the reference isolate (M. smegmatis: 1 reference + 8 SPIRE MAGs,
# the leader goes to a MAG; M. fallax: 1 reference + 4 MAGs; M. tuberculosis,
# M. leprae, M. avium: one tip each).
# M. abscessus carries no label because it sits outside the 824-genome crown
# (in the basal Mycobacteroides-grade lineage), so a label would sit beside the
# wrong arc. GTDB-Tk R226 returns g__Mycobacterium for these genomes; the
# exclusion is topological, not taxonomic.
KEY = ["Mycobacterium tuberculosis", "Mycobacterium leprae",
       "Mycobacterium avium", "Mycobacterium smegmatis",
       "Mycobacterium fallax"]
species = meta.reindex(tips)["species"].astype(str)
marks = []
for sp in KEY:
    hit = [i for i, s in enumerate(species) if s.strip() == sp]
    if not hit:
        hit = [i for i, s in enumerate(species) if s.startswith(sp)]
    if not hit:
        continue
    i = hit[len(hit) // 2]
    a = START_DEG - (i + 0.5) * dtheta
    marks.append((sp.replace("Mycobacterium", "M."), a))

# Labels are set horizontally with two-segment leaders. Rotated diagonal text
# collided at the SGM/FSM boundary, where M. leprae, M. tuberculosis and
# M. fallax sit within a few degrees of one another.
R_ELBOW = r_gen + 0.105
R_LABEL = r_gen + 0.160
MIN_GAP = 0.105          # axis units; ~1.4 line heights at 7 pt

placed = []
for side in (+1, -1):
    grp = [(sp, a) for sp, a in marks if np.cos(np.deg2rad(a)) * side >= 0]
    grp = [(sp, a, xy(R_LABEL, a)[1]) for sp, a in grp]
    grp.sort(key=lambda t: t[2])
    ys = [t[2] for t in grp]
    # Single downward-then-upward relaxation pass: enough for five labels.
    for k in range(1, len(ys)):
        if ys[k] - ys[k - 1] < MIN_GAP:
            ys[k] = ys[k - 1] + MIN_GAP
    for k in range(len(ys) - 2, -1, -1):
        if ys[k + 1] - ys[k] < MIN_GAP:
            ys[k] = ys[k + 1] - MIN_GAP
    for (sp, a, _), y_lab in zip(grp, ys):
        placed.append((sp, a, side, y_lab))

for sp, a, side, y_lab in placed:
    x0, y0 = xy(r_gen + 0.012, a)
    xe, ye = xy(R_ELBOW, a)
    x_lab = xy(R_LABEL, a)[0]
    x_lab = side * max(abs(x_lab), abs(xe) + 0.02)
    ax.plot([x0, xe, x_lab], [y0, ye, y_lab], color="#6E6E6E", lw=0.5,
            solid_joinstyle="round", zorder=5)
    ax.text(x_lab + side * 0.022, y_lab, sp, fontsize=7, style="italic",
            color=INK, zorder=6, va="center",
            ha="left" if side > 0 else "right")

# Scale bar for the linear part of the radial axis, set in the empty corner
# below the circle. Placed inside the wedge at the top it landed on the RGM
# branches.
BAR_SUBS = 0.1
bar_len = (R_OUT - R_BREAK) * BAR_SUBS / (dmax - D_CROWN)
bx0, by0 = -LIM * 0.90, -LIM * 0.80
ax.plot([bx0, bx0 + bar_len], [by0, by0], color=INK, lw=1.1,
        solid_capstyle="butt", zorder=6)
for tx in (bx0, bx0 + bar_len):
    ax.plot([tx, tx], [by0 - 0.022, by0 + 0.022], color=INK, lw=1.1, zorder=6)
ax.text(bx0 + bar_len / 2, by0 - 0.058, f"{BAR_SUBS:g} substitutions per site",
        fontsize=6.4, color=INK, ha="center", va="top", zorder=6)

lim = LIM
ax.set_xlim(-lim, lim)
ax.set_ylim(-lim, lim)

# ------------------------------------------------------------------ legend
# Three columns under the tree. The legend used to be a tall right-hand
# column, which limited the circle to the height left over beside it.
lax = fig.add_axes([0.035, 0.012, 0.945, LEG_H - 0.030])
lax.axis("off")
lax.set_xlim(0, 1)
lax.set_ylim(0, 1)

COLX = [0.000, 0.345, 0.680]
LH, HH, SECGAP = 0.083, 0.098, 0.030
ytop = 0.985
ycur = {i: ytop for i in range(3)}


def head(col, text):
    lax.text(COLX[col], ycur[col], text, fontsize=8, fontweight="bold",
             va="top", color=INK)
    ycur[col] -= HH


def swatch(col, colour, label):
    y = ycur[col]
    lax.add_patch(Rectangle((COLX[col] + 0.004, y - LH * 0.86), 0.030, LH * 0.62,
                            facecolor=colour,
                            edgecolor="#9A9A9A" if colour == ABSENT_C else "none",
                            linewidth=0.4,
                            transform=lax.transAxes, clip_on=False))
    lax.text(COLX[col] + 0.045, y - LH * 0.55, label, fontsize=7,
             va="center", color=INK)
    ycur[col] -= LH


def dot(col, colour, label):
    y = ycur[col]
    lax.plot([COLX[col] + 0.019], [y - LH * 0.55], marker="o", ms=2.8,
             color=colour, transform=lax.transAxes, clip_on=False)
    lax.text(COLX[col] + 0.045, y - LH * 0.55, label, fontsize=7,
             va="center", color=INK)
    ycur[col] -= LH


def colourbar(col, key):
    y = ycur[col]
    lo, hi, cmap_name = ring_bounds[key]
    grad = np.linspace(0, 1, 256).reshape(1, -1)
    sub = lax.inset_axes([COLX[col] + 0.004, y - LH * 0.95, 0.190, LH * 0.62])
    sub.imshow(grad, aspect="auto",
               cmap=mpl.colors.LinearSegmentedColormap.from_list(
                   cmap_name + "_trunc",
                   plt.get_cmap(cmap_name)(np.linspace(0.18, 0.74, 256))))
    sub.set_xticks([]), sub.set_yticks([])
    for sp_ in sub.spines.values():
        sp_.set_linewidth(0.5)
    lax.text(COLX[col] + 0.004, y - LH * 1.15, f"\u2264 {lo:.1f}", fontsize=6.6,
             ha="left", va="top", color=INK)
    lax.text(COLX[col] + 0.194, y - LH * 1.15, f"\u2265 {hi:.1f}", fontsize=6.6,
             ha="right", va="top", color=INK)
    ycur[col] -= LH * 2.1


def note(col, label):
    """Key line without a marker (for the 'no dot' category)."""
    y = ycur[col]
    lax.text(COLX[col] + 0.045, y - LH * 0.55, label, fontsize=7, va="center",
             color=INK)
    ycur[col] -= LH


n_mag = int((meta["culture_status"] == "MAG / uncultured").sum())
n_mag_spire = int(((meta["culture_status"] == "MAG / uncultured") &
                   (meta["origin"] == "SPIRE MAG")).sum())
n_mag_ref = n_mag - n_mag_spire
n_basal = int((meta["clade3"] == "Outside").sum())
print(f"MAG-coloured tips {n_mag} = {n_mag_spire} SPIRE + {n_mag_ref} GTDB reference; "
      f"basal Mycobacteroides-grade tips {n_basal}")

head(0, "Clade (arcs)")
for k, c in [("SGM  slow-growing clade", SGM_C),
             ("RGM  large-genome grade (paraphyletic)", RGM_C),
             ("FSM  small genome, mainly freshwater", FSM_C)]:
    swatch(0, c, k)
ycur[0] -= SECGAP
head(0, "Node support (UFBoot)")
dot(0, "#111111", "\u2265 95")
dot(0, "#9A9A9A", "90\u201394")
note(0, "< 90: no dot")

head(1, "Branch colour (tip category)")
for k, lab in [("Cultured isolate", "Cultured isolate"),
               ("NTM (isolate)", "NTM isolate (cultured)"),
               ("MAG / uncultured", f"MAG: {n_mag_spire} SPIRE + {n_mag_ref} GTDB refs"),
               ("Tissue-derived (host-adapted pathogen)",
                "Tissue-derived (never cultured)")]:
    swatch(1, BRANCH_C[k], lab)
swatch(1, BRANCH_C["Outgroup"], "Rooting taxa (Tsukamurella)")
swatch(1, INTERNAL_C, f"Internal branches; {n_basal} basal tips*")
ycur[1] -= SECGAP * 0.4
lax.text(COLX[1], ycur[1],
         "* Mycobacteroides-grade genomes outside the\n  824-genome crown; always grey.",
         fontsize=6.4, va="top", color=INK)

head(2, "Ring 1  Genome size (Mb)")
colourbar(2, "genome_size_Mb")
ycur[2] -= SECGAP * 0.5
head(2, "Rings 2\u20134  Gene presence")
swatch(2, PRESENT["urease"], "Urease")
swatch(2, PRESENT["ectoine_synthase"], "Ectoine synthase")
swatch(2, PRESENT["rhodopsin"], "Microbial rhodopsin (opsin gene)")
swatch(2, ABSENT_C, "Absent")
ycur[2] -= SECGAP * 0.6
lax.text(COLX[2], ycur[2],
         "862 genomes: 246 SPIRE MAGs + 616 GTDB\n"
         f"references ({n_mag_ref} of them MAGs, in the MAG colour).\n"
         "Radial scale: substitutions per site above the crown of\n"
         "the 843 non-Tsukamurella genomes (stem compressed).\n"
         "Clade arcs cover the 824-genome crown only.\n"
         "GC: Fig. 2a; proteome acidity: Tables S14, S1.",
         fontsize=6.4, va="top", color=INK, linespacing=1.25)

save(fig, "Figure1_phylogeny_traits")
for cl_name in ["SGM", "RGM", "FSM", "Outgroup"]:
    sub = meta[meta.clade3 == cl_name]
    print(f"{cl_name:9s} n={len(sub):4d}  rhodopsin+={100*(sub.rhodopsin=='present').mean():5.1f}%"
          f"  median size={sub.genome_size_Mb.median():.2f} Mb")
