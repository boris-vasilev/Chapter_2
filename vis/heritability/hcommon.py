"""Shared helpers for the heritability figures (Fig. 4 and supplementary).

Puts chapter2/vis on the import path so ``common`` is shared with the other
chapter 2 figures, and defines the annotation order, labels and the loader
for the tidy table written by 01_collect_results.py.
"""

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from common import (BUILD_DIR, add_variant_args,  # noqa: E402,F401
                    load_diseases, resolve_variant)
# The design setting and every input / output root, shared with common.py
# (single source: set CH2_VARIANT to switch design).
from common import (FIG_DIR, LDSC_DIR, PERM_DIR, RESULTS, ROOT,  # noqa: E402,F401
                    VARIANT)

# Column order used on every heritability figure (brief: anchors, regulators,
# targets). Colours are the chapter's ANNOT_COLORS.
H_ANNOTS = ["H4_genes", "regulators_H4", "targets_H4"]
H_ANNOT_LABEL = {"H4_genes": "Anchors", "regulators_H4": "Regulators",
                 "targets_H4": "Targets"}
FDR = 0.05
LONG_CSV = "heritability_long.csv"

# Status of each trait x condition x annotation cell (01_collect_results.py).
FITTED = "fitted"
EMPTY = "empty set"
DROPPED = "dropped (< 1,000 SNPs)"
DROPPED_GENES = "dropped (< 3 genes)"
FAILED = "run failed"
NOT_RUN = "not run"


def dropped_status(combine_txt, annot):
    """Why combine_annots.py dropped ``annot``, from its [DROP] line in
    <cond>_combined.combine.txt ("[DROP] <annot>: <n> SNPs, <k> genes").
    A set with enough SNPs was dropped by the gene rule; anything else
    (too few SNPs, or an older combine.txt with no gene count) is DROPPED."""
    m = re.search(rf"\[DROP\] {annot}: (\d+) SNPs, (\d+) genes",
                  combine_txt.read_text()) if combine_txt.exists() else None
    min_snps = int(re.search(r"_snp(\d+)", VARIANT).group(1))
    return DROPPED_GENES if m and int(m.group(1)) >= min_snps else DROPPED


def add_args(ap):
    add_variant_args(ap)
    ap.add_argument("--long", default=None,
                    help=f"Tidy table from 01_collect_results.py (default: "
                         f"<out-dir>/{LONG_CSV}).")


def resolve(args):
    """Variant paths as in common.resolve_variant, but staging in a
    heritability/ subfolder of build/<VARIANT>/ unless --out-dir is given."""
    if args.out_dir is None:
        args.out_dir = str(BUILD_DIR / "heritability")
    args = resolve_variant(args)
    args.long = Path(args.long or args.out_dir / LONG_CSV)
    return args


def load_long(path):
    df = pd.read_csv(path)
    df["log2_enrichment"] = np.log2(df["enrichment"].where(df["enrichment"] > 0))
    return df


def flagged_name(diseases, study):
    """Display name, with an asterisk when the GWAS ancestry differs from the
    LD reference panel used by S-LDSC (multi-ancestry GWAS on the EUR panel)."""
    n = diseases.loc[study, "display_name"]
    return n + " *" if diseases.loc[study, "ld_mismatch"] else n
