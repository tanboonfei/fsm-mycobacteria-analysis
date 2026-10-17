#!/usr/bin/env python3
"""Supplementary Figure S2: gene-content clustering of the ingroup genomes.

UPGMA (average-linkage) dendrogram of the 822 ingroup genomes with gene
predictions (FSM 134, RGM 389, SGM 299; the 824-genome ingroup minus two RGM
assemblies not re-called under the uniform Prodigal pipeline) from gene-family
presence/absence (MMseqs2 clusters, Jaccard distance over ALL 147,225 families
with 2 <= size < 822), with the 2,477 clade-discriminating families drawn as a
presence heatmap. Row colour = the clade in which the family is most prevalent
(not the clade of the genome); coloured cell = present, white = absent.

Family selection: prevalence max - min over FSM, RGM, SGM >= 0.5. The families
are chosen using the clade labels, so the block-diagonal colour structure is
built in and is not independent evidence that the clustering recovers the
clades. The recovery that IS independent is read from the leaf order of the
dendrogram, printed at the end (FSM is recovered except for one genome; RGM and
SGM are only partly separated).

Inputs : inputs/S1_genome_metadata.tsv   (Table S1, clade calls)
         $MYCO_UNI_DIR/uni824_cluster.tsv and disc.emapper.annotations
         (default: publication/analysis_uniform_clustering, 307 MB, not copied)
Output : supplementary_figures/FigureS2_genecontent_full.{png,pdf}
         The 2,477-family list is distributed as scripts/inputs/tsv/S19_discriminating_families.tsv
         (old table number S19; this is not a table of Additional file 2).
Usage  : python3 build_figureS2_genecontent.py
"""
import csv
import os
from pathlib import Path
from collections import Counter, defaultdict

import matplotlib.colors as mc
import matplotlib.patches as mp
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist

from figstyle import CLADE_C, INK, INPUTS, WIDTH_FULL, save

UNI = Path(os.environ.get("MYCO_UNI_DIR", "/path/to/project/"
                                          "publication/analysis_uniform_clustering"))

MIN_PREV_DIFF = 0.5
ORDER_CLADES = ["FSM", "RGM", "SGM"]
COGNAME = {"C": "Energy production", "E": "Amino-acid metab", "P": "Ion transport",
           "Q": "Secondary metabolites", "G": "Carbohydrate metab", "I": "Lipid metab",
           "K": "Transcription", "L": "Replication/repair", "H": "Coenzyme metab",
           "J": "Translation", "M": "Cell wall/envelope", "T": "Signalling",
           "F": "Nucleotide metab", "V": "Defense", "O": "PTM/chaperones",
           "U": "Secretion", "N": "Motility"}

# eggNOG annotation of the uniform clustering's discriminating representatives;
# the descriptions printed on the figure are the block-level COG tallies of these.
annf = UNI / "disc.emapper.annotations"
ANN = {}
for line in open(annf):
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    ANN[f[0]] = (f[8], f[6])          # Preferred_name, COG_category


def annof(rep):
    return ANN.get(rep) or ANN.get(rep.replace("|", "__"))


# Clade calls come from Table S1, and the families from the uniform Prodigal
# clustering of the 822 ingroup genomes (analysis_uniform_clustering/README.txt).
clade = {}
for r in csv.DictReader(open(INPUTS / "S1_genome_metadata.tsv"), delimiter="\t"):
    c = r["clade"].split(" ")[0]
    if c in ("SGM", "RGM", "FSM"):
        clade[r["genome_id_or_accession"]] = c

fam = defaultdict(set)
for line in open(UNI / "uni824_cluster.tsv"):
    rep, mem = line.rstrip("\n").split("\t")
    fam[rep].add(mem.split("|")[0])

# Restrict to the genomes that actually entered the uniform clustering: the two
# ingroup assemblies not re-called under that pipeline carry no proteins here and
# would otherwise be counted as genomes missing every family.
in_clustering = {g for members in fam.values() for g in members}
clade = {g: c for g, c in clade.items() if g in in_clustering}

G = sorted(clade)
gi = {g: i for i, g in enumerate(G)}
nG = len(G)
keep = [f for f, s in fam.items() if 2 <= len(s) < nG]
M = np.zeros((nG, len(keep)), np.int8)
for j, f in enumerate(keep):
    for g in fam[f]:
        if g in gi:
            M[gi[g], j] = 1
clades = np.array([clade[g] for g in G])

Z = linkage(pdist(M, "jaccard"), "average")
gorder = dendrogram(Z, no_labels=True, no_plot=True)["leaves"]

cl_prev = np.vstack([M[clades == c].mean(0) for c in ORDER_CLADES])
sel = np.where(cl_prev.max(0) - cl_prev.min(0) >= MIN_PREV_DIFF)[0]
assign = np.array(ORDER_CLADES)[cl_prev[:, sel].argmax(0)]
maxp = cl_prev[:, sel].max(0)
rank = {c: i for i, c in enumerate(ORDER_CLADES)}
o = sorted(range(len(sel)), key=lambda k: (rank[assign[k]], -maxp[k]))
sel, assign = sel[o], assign[o]

rgb = np.ones((len(sel), nG, 3))
for ri, (fj, a) in enumerate(zip(sel, assign)):
    rgb[ri, M[gorder, fj].astype(bool)] = mc.to_rgb(CLADE_C[a])

groups = [(c, int((assign == c).sum())) for c in ORDER_CLADES]

blockfun = {}
y0 = 0
for c, n in groups:
    cogs, names, nann = Counter(), Counter(), 0
    for k in range(y0, y0 + n):
        a = annof(keep[sel[k]])
        if a:
            nann += 1
            if a[0] not in ("-", ""):
                names[a[0]] += 1
            for ch in (a[1] or ""):
                if ch in COGNAME:
                    cogs[ch] += 1
    blockfun[c] = (cogs.most_common(3), [g for g, _ in names.most_common(5)], nann)
    y0 += n

# ---------------------------------------------------------------- plot
fig = plt.figure(figsize=(WIDTH_FULL, 5.10))
gs = GridSpec(3, 1, height_ratios=[1.5, 0.16, 8], hspace=0.04,
              left=0.115, right=0.605, top=0.955, bottom=0.095)

axd = fig.add_subplot(gs[0])
dendrogram(Z, no_labels=True, ax=axd, color_threshold=0,
           above_threshold_color="#8A8A8A")
axd.set_xticks([])
axd.set_yticks([])
for s in axd.spines:
    axd.spines[s].set_visible(False)
axd.set_title("Gene-content clustering of the 822 ingroup genomes", loc="left",
              fontsize=8, color=INK, pad=4)

axc = fig.add_subplot(gs[1])
for i, l in enumerate(gorder):
    axc.add_patch(plt.Rectangle((i, 0), 1, 1, color=CLADE_C[clades[l]], lw=0))
axc.set_xlim(0, nG)
axc.set_ylim(0, 1)
axc.axis("off")

axh = fig.add_subplot(gs[2])
axh.imshow(rgb, aspect="auto", interpolation="nearest", rasterized=True)
axh.set_xticks([])
axh.set_yticks([])

y = 0
for c, n in groups:
    axh.text(-nG * 0.02, y + n / 2, f"{c}-enriched\ngene families\n({n})",
             ha="right", va="center", fontsize=7, fontweight="bold",
             color=CLADE_C[c], linespacing=1.4)
    cogs, genes, nann = blockfun[c]
    func = "Top functions: " + "; ".join(COGNAME[ch] for ch, _ in cogs)
    egtxt = ("e.g. " + ", ".join(genes[:5])) if genes else ""
    # Near-black at 7 pt in its own column: the shipped version set this in the
    # clade colour at a size that did not survive print.
    axh.text(nG * 1.035, y + n / 2, f"{func}\n{egtxt}", ha="left", va="center",
             fontsize=7, color=INK, linespacing=1.5)
    y += n
    if y < len(sel):
        axh.axhline(y - 0.5, color="white", lw=1.2)
axh.set_xlabel(f"{nG} genomes in dendrogram order (tree built on all {len(keep):,} families). Rows: {len(sel):,} families with\n"
               "clade prevalence gap \u2265 0.5, selected using the clade labels; row colour = clade of highest prevalence; white = absent.",
               fontsize=6.6, linespacing=1.3)

axd.legend(handles=[mp.Patch(color=CLADE_C[c],
                             label=f"{c} (n = {int((clades == c).sum())})")
                    for c in ORDER_CLADES],
           loc="upper left", bbox_to_anchor=(1.02, 1.06), frameon=False,
           fontsize=7, title="Clade", title_fontsize=7.5, handlelength=1.1,
           handletextpad=0.5, labelspacing=0.35)

save(fig, "FigureS2_genecontent_full", supp=True)
print(f"{nG} genomes, {len(sel)} discriminating families (prevalence diff >= "
      f"{MIN_PREV_DIFF}); groups: {groups}")
print("annotation source:", annf)
print(f"families used for the dendrogram: {len(keep):,}")
# independent check of recovery: runs of clade labels along the leaf order
labs = [clades[l] for l in gorder]
runs = []
for lab in labs:
    if runs and runs[-1][0] == lab:
        runs[-1][1] += 1
    else:
        runs.append([lab, 1])
print("leaf-order runs:", ", ".join(f"{c}{n}" for c, n in runs))
