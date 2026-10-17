"""Shared paths and the Table S1 loader.

The figure code is not part of this package; this module carries the few
definitions the analysis scripts took from it.
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parent
# Everything the scripts read that is not in Additional file 2 is bundled here.
INPUTS = ROOT / "scripts" / "inputs"
# inputs/tsv mirrors the workbook and is written only by export_workbook_tsv.py;
# a rerun of any analysis script writes its results here instead
TSV = INPUTS / "tsv"
RUNOUT = ROOT / "scripts" / "run_output"
RUNOUT.mkdir(exist_ok=True)

ORDER = ["SGM", "RGM", "FSM"]


def load_s1(pd):
    s1 = pd.read_csv(TSV / "S1_genome_metadata.tsv", sep="\t")
    s1["clade3"] = s1["clade"].str.split(" ").str[0]
    return s1
