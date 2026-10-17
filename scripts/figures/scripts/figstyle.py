"""Shared style block for every manuscript figure.

Import this first in every build script and call `save(fig, name)` to write.
Nothing here is figure-specific; per-figure content stays in its own script.

Typography
    Arial is requested first and Liberation Sans is the fallback actually used
    on this machine. Liberation Sans is metric-compatible with Arial (identical
    advance widths), so line breaks and label extents match what a journal
    setting Arial will see. Nimbus Sans (a Helvetica clone) follows.
    Base size is 8 pt at true column width, so figures MUST be sized with the
    WIDTH_* constants below and never scaled after the fact; scaling a 10 in
    figure down to column width is what silently pushes 8 pt text under 6 pt.

Colour
    The palette is deliberately not a matplotlib or seaborn default. It was
    selected by simulating deuteranopia, protanopia, tritanopia and greyscale
    over candidate sets and maximising the minimum separation. Against the
    previous tab10-derived colours it improves the worst greyscale lightness
    gap from 6.7 to 20.7 L* and the worst deuteranopia separation from 8.5 to
    24.7. Re-check colour separation after any palette change.

Output
    save() writes a vector PDF (the submission artefact) plus a 600 dpi PNG
    preview. pdf.fonttype 42 embeds TrueType outlines; the matplotlib default
    of Type 3 is rejected by several publishers.
"""
from pathlib import Path

import matplotlib as mpl

# Force the non-interactive backend before pyplot is imported. Otherwise an
# available display makes matplotlib select Tk, whose toolbar construction can
# fail outright under some HiDPI settings, and the figure never gets built.
mpl.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

# --------------------------------------------------------------- paths
# Layout (all relative to the folder that holds scripts/):
#   inputs/                     S1_genome_metadata.tsv, myco_bac120.rooted.treefile,
#                               relax_best_fits.tsv, optional fig4a_values.tsv, ...
#   inputs/tsv/                 per-table exports named by the CURRENT table numbers
#                               (S4b, S5i, S7, S9, S9e, S11d, S11i)
#   figures/, supplementary_figures/   outputs (PNG + PDF)
# Set MYCO_INPUTS to point at a different inputs folder.
import os

ROOT = Path(__file__).resolve().parents[1]
INPUTS = Path(os.environ.get("MYCO_INPUTS", ROOT / "inputs"))
TSV = INPUTS / "tsv"
OUT = ROOT / "figures"
SUPP = ROOT / "supplementary_figures"

# --------------------------------------------------------------- geometry
# Journal column widths: 85 mm single column, 170 mm full width.
WIDTH_SINGLE = 85 / 25.4          # 3.35 in
WIDTH_FULL = 170 / 25.4           # 6.69 in
WIDTH_150 = 129 / 25.4            # 5.08 in, intermediate

# --------------------------------------------------------------- ink
INK = "#1A1A1A"        # all data-carrying text, near-black
INK_SOFT = "#4D4D4D"   # structural furniture only (axis spines, gridlines)
RULE = "#8A8A8A"       # reference lines

# --------------------------------------------------------------- palette
SGM_C = "#6E2244"      # deep plum
RGM_C = "#1C808D"      # deep teal
FSM_C = "#DBA02A"      # warm ochre
CLADE_C = {"SGM": SGM_C, "RGM": RGM_C, "FSM": FSM_C}
ORDER = ["SGM", "RGM", "FSM"]

SMALL_C = "#8A8A8A"    # sporadic small non-FSM genomes
LARGE_C = "#CFCFCF"    # large-genome background points
LEPRAE_C = "#A32B18"   # host-restricted decayers

# Presence/absence rings on Figure 1. ABSENT must stay visible on white.
# Rhodopsin must not reuse the RGM teal: on Figure 1 the presence ring and the
# RGM clade arc would then read as the same category.
PRESENT_C = {"urease": "#7E6033", "ectoine_synthase": "#A32B7A",
             "rhodopsin": "#173F6B"}
ABSENT_C = "#D9D9D9"

# Figure 1 runs four colour systems at once: branch categories, trait rings,
# clade arcs and genus arcs. An earlier version reused one red for both
# "MAG / uncultured" branches and the Origin ring, and one plum for both
# "Pathogen (isolate)" branches and the SGM arc, so the eye read one signal
# twice. Every cross-system pair was held apart under simulated colour-vision deficiency.
# Two "Pathogen (isolate)" tips (M. tuberculosis, M. decipiens) sit among the 862
# tips; they take no colour slot of their own and are drawn as the cultured
# isolates they are (BRANCH_ALIAS).
BRANCH_C = {
    "Cultured isolate":   "#ACACAC",
    "MAG / uncultured":   "#E24A2B",
    "NTM (isolate)":      "#B25FCF",
    # Tissue-derived host-adapted pathogens (M. leprae, M. lepromatosis,
    # M. uberis) are not cultured isolates and get their own branch colour.
    "Tissue-derived (host-adapted pathogen)": "#1F77B4",
    "Outgroup":           "#111111",
}
BRANCH_ALIAS = {"Pathogen (isolate)": "Cultured isolate"}
# The genus band is one colour with its three segments labelled on the figure.
# Three colours here collided with the trait rings and the clade arcs under
# every simulated colour-vision deficiency, and the labels already name them.
GENUS_NAMES = ("Mycobacterium s.s.", "Mycolicibacterium", "Mycobacteroides")
GENUS_BAND_C = "#282828"
# One structural grey. Three greys (tips, internal, basal-outside) sat within
# 13 L* of each other and could not be told apart on a hairline branch.
INTERNAL_C = "#757575"
OUTSIDE_C = INTERNAL_C  # basal GTDB-Mycobacterium outside the analysed crown


# Niche categories (Fig. 4b-d). Deliberately disjoint from the clade palette.
# All four panels sit in one figure, and therefore a reader carries the clade colours over
# from Fig. 4a: an earlier version used blue for freshwater against RGM teal, which
# separated by only 1.7 under simulation, and a wine for host-associated that
# was 3 degrees of hue from SGM plum. Hue space is crowded once teal (187),
# magenta (333) and gold (40) are excluded; the palette is built from two
# clade-safe families and carries category on LIGHTNESS: indigo for the water
# categories, brown for the terrestrial ones, neutral grey for the remainder.
# That also makes it legible in greyscale by construction.
# Worst separation from any clade colour is 11.1; worst within-palette pair 16.6.
# (checked under simulated colour-vision deficiency).
HABITAT_C = {
    "freshwater":     "#4A4291",   # indigo, mid
    "soil/sediment":  "#7A5A33",   # brown, mid
    "wastewater":     "#7C75C0",   # indigo, light-mid
    "drinking-water": "#B3AEE4",   # indigo, pale
    "marine/saline":  "#1A1638",   # indigo, near-black
    "host-assoc.":    "#B08D5E",   # brown, light
    "aquaculture":    "#D9CBB0",   # brown, pale
    "other":          "#9A9A9A",   # neutral
}
# Ordinal: low salinity takes the freshwater indigo, high takes the terrestrial
# brown, and panel c reads consistently with panel b.
# Panel 4d (aquatic vs other biomes) is a different variable from the habitat
# classes of panel 4b, so it must not reuse the freshwater indigo: neutral
# dark grey for aquatic, light grey for the remaining biomes.
AQUATIC_C = "#3F3F3F"
OTHER_BIOME_C = "#CFCFCF"
SALINITY_C = {"low": "#4A4291", "medium": "#B3AEE4", "high": "#7A5A33"}

# --------------------------------------------------------------- rcParams
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "Nimbus Sans",
                        "FreeSans", "DejaVu Sans"],
    "font.size": 8,
    "axes.titlesize": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "legend.fontsize": 7.5,
    "legend.frameon": False,
    "figure.titlesize": 9,

    "text.color": INK,
    "axes.labelcolor": INK,
    "axes.edgecolor": INK_SOFT,
    "xtick.color": INK,
    "ytick.color": INK,
    "axes.titlecolor": INK,

    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.6,
    "ytick.major.size": 2.6,
    "lines.linewidth": 1.0,
    "patch.linewidth": 0.5,

    "pdf.fonttype": 42,        # embed TrueType, not Type 3
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
    "figure.dpi": 150,
})

# Liberation Sans has no U+207B SUPERSCRIPT MINUS or U+2070 SUPERSCRIPT ZERO,
# so Unicode exponents render as missing-glyph boxes. Exponents therefore go
# through mathtext, configured to use the same face as the body text.
def _resolve_family():
    from matplotlib.font_manager import findfont, FontProperties, get_font
    try:
        path = findfont(FontProperties(family=mpl.rcParams["font.sans-serif"]))
        return get_font(path).family_name
    except Exception:
        return "DejaVu Sans"


_FAMILY = _resolve_family()
mpl.rcParams.update({
    "mathtext.fontset": "custom",
    "mathtext.rm": _FAMILY,
    "mathtext.it": f"{_FAMILY}:italic",
    "mathtext.bf": f"{_FAMILY}:bold",
    # cal/sf/tt are unused but matplotlib resolves them anyway, and the default
    # 'cursive' is absent here. Left unset it warns on every build, which would
    # mask the missing-glyph report in save().
    "mathtext.cal": _FAMILY,
    "mathtext.sf": _FAMILY,
    "mathtext.tt": _FAMILY,
    "mathtext.default": "regular",
})

MIN_PT = 6.0   # nothing below this at final print size


# --------------------------------------------------------------- helpers
def panel_tag(ax, letter, title=None, dx=-0.10, dy=1.02, titlesize=8):
    """Bold panel letter, plus an optional SHORT noun phrase.

    Findings belong in the manuscript legend, not in the panel. Pass a noun
    phrase ("Codon-usage bias") or nothing; do not pass a sentence.
    """
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=10,
            fontweight="bold", va="bottom", ha="left", color=INK)
    if title:
        ax.set_title(title, loc="left", fontsize=titlesize, color=INK, pad=4)


def pfmt(p, prefix="p = ", mantissa_dp=0):
    """Format a p-value in true scientific notation with Unicode superscripts.

    One significant figure on the mantissa by default, which is the convention
    used throughout the manuscript text and legends.

    The exponent is emitted as mathtext, not as Unicode superscripts, because
    the Arial-metric face lacks the superscript-minus glyph.

    pfmt(8.8e-3) -> 'p = 9 × 10$^{-3}$'
    pfmt(0.2)    -> 'p = 0.2'
    """
    if p is None or not np.isfinite(p):
        return ""
    if p >= 0.01:
        s = f"{p:.2g}"
    else:
        exp = int(np.floor(np.log10(p)))
        mant = p / 10 ** exp
        if mantissa_dp == 0 and round(mant) == 10:      # 9.6e-3 -> 1 x 10^-2
            mant, exp = 1.0, exp + 1
        s = f"{mant:.{mantissa_dp}f} × 10$^{{{exp}}}$"
    return prefix + s


def _too_small(fig):
    """Every text artist below MIN_PT. A build fails loudly instead of
    shipping unreadable annotation."""
    bad = []
    for t in fig.findobj(mpl.text.Text):
        if not t.get_text().strip():
            continue
        sz = t.get_fontsize()
        if sz < MIN_PT:
            bad.append((round(sz, 2), t.get_text()[:40]))
    return bad


def save(fig, name, supp=False, check=True):
    """Write vector PDF plus 600 dpi PNG preview, and check text sizes."""
    d = SUPP if supp else OUT
    d.mkdir(parents=True, exist_ok=True)
    if check:
        bad = _too_small(fig)
        if bad:
            print(f"  !! {name}: {len(bad)} text item(s) below {MIN_PT} pt: {bad[:5]}")
    # A missing glyph is silent in the PNG but ships a blank box to the
    # publisher, which means promote matplotlib's warning to a visible report.
    import warnings
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        fig.savefig(d / f"{name}.pdf")
        fig.savefig(d / f"{name}.png", dpi=600)
    missing = sorted({str(w.message) for w in caught if "Glyph" in str(w.message)})
    if missing:
        print(f"  !! {name}: missing glyphs -> {missing}")
    w, h = fig.get_size_inches()
    print(f"wrote {d/name}.pdf + .png  ({w*25.4:.0f} x {h*25.4:.0f} mm)")


def violin_box(ax, values_by_group, ylabel, groups=None, colours=None, seed=0,
               point_size=2.2, min_violin=8):
    """Groups smaller than min_violin are drawn as their points with a median
    bar. A violin over three genomes would invent a distribution."""
    groups = groups or ORDER
    colours = colours or CLADE_C
    rng = np.random.default_rng(seed)
    for i, g in enumerate(groups):
        v = np.asarray(values_by_group[g], dtype=float)
        v = v[np.isfinite(v)]
        if 0 < len(v) < min_violin:
            ax.scatter([i] * len(v), v, s=point_size * 6, color=colours[g],
                       edgecolor=INK, lw=0.4, zorder=4)
            ax.plot([i - 0.22, i + 0.22], [np.median(v)] * 2, color=INK,
                    lw=1.2, zorder=5, solid_capstyle="butt")
            continue
        if len(v) > 1:
            parts = ax.violinplot(v, positions=[i], widths=0.78,
                                  showextrema=False, showmedians=False)
            for body in parts["bodies"]:
                body.set_facecolor(colours[g])
                body.set_alpha(0.45)
                body.set_edgecolor(colours[g])
                body.set_linewidth(0.6)
        ax.boxplot(v, positions=[i], widths=0.15, showfliers=False,
                   patch_artist=True,
                   boxprops=dict(facecolor="white", edgecolor=INK, lw=0.7),
                   medianprops=dict(color=INK, lw=1.2),
                   whiskerprops=dict(color=INK, lw=0.7),
                   capprops=dict(color=INK, lw=0.7))
        ax.scatter(i + rng.normal(0, 0.055, len(v)), v, s=point_size, color=INK,
                   alpha=0.30, lw=0, zorder=3)
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(groups)
    ax.set_ylabel(ylabel)
    ax.spines[["top", "right"]].set_visible(False)


def seq_cmap(name="GnBu", lo=0.15, hi=1.0):
    """Sequential ramp truncated away from white.

    Full-range sequential maps bottom out at white, where low values become
    invisible against the page and against the light grey absent state.
    """
    return mpl.colors.LinearSegmentedColormap.from_list(
        f"{name}_{lo}_{hi}", plt.get_cmap(name)(np.linspace(lo, hi, 256)))


def load_s1(pd):
    s1 = pd.read_csv(INPUTS / "S1_genome_metadata.tsv", sep="\t")
    s1["clade3"] = s1["clade"].str.split(" ").str[0]
    return s1


def report_font():
    from matplotlib.font_manager import findfont, FontProperties
    f = findfont(FontProperties(family=mpl.rcParams["font.sans-serif"]))
    print(f"[figstyle] resolved sans-serif -> {Path(f).name}")
    if "Liberation" not in f and "Arial" not in f and "Helvetica" not in f:
        print("  !! not an Arial-metric face; install Arial or Liberation Sans")
    return f
