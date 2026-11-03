from pathlib import Path
import polars as pl


# --- Path resolution ------------------------------------------------

WORKSPACE = Path(config["workspace"]).resolve()


def _ws(rel):
    return str(WORKSPACE / rel)


# Module selection — module-specific paths are read from the matching
# config block (config['coloc'], config['l2g'], ...).
MODULE = config["gene_set_module"]
mod_cfg = config[MODULE]

# Shared / module-agnostic.
DE_INPUT_H5AD = _ws(config["de_input_h5ad"])

# Per-module raw and munged sumstats locations.
STUDY_INFO   = _ws(mod_cfg["study_info"])
GWAS_DIR     = _ws(mod_cfg["gwas_dir"])
SUMSTATS_DIR = _ws(mod_cfg["sumstats_dir"])

# Coloc-specific input — read from the *active* module so coloc_strict
# (or any future coloc-variant) can point at a different coloc table.
COLOC_DAT = _ws(mod_cfg.get("coloc_dat", "data/coloc_dat.tsv"))

# L2G-specific inputs.
L2G_DIR = _ws(config["l2g"].get("data_dir", "data/opentargets"))

# Derived outputs — every module gets its own subtree so coloc and
# l2g (and any future module) never collide.
DE_RESULTS             = _ws(config["de_results"])
GENE_SETS_DIR          = _ws(f"{config['gene_sets_dir']}/{MODULE}")
ANNOTATIONS_DIR        = _ws(f"{config['annotations_dir']}/{MODULE}")
LD_SCORES_DIR          = _ws(f"{config['ld_scores_dir']}/{MODULE}")
LD_SCORES_FILTERED_DIR = _ws(f"{config['ld_scores_filtered_dir']}/{MODULE}")
LDSC_RESULTS_DIR       = _ws(f"{config['ldsc_results_dir']}/{MODULE}")


# --- LDSC tool paths (inside container) -----------------------------

LDSC           = f"{config['ldsc_dir']}/ldsc.py"
MUNGE_SUMSTATS = f"{config['ldsc_dir']}/munge_sumstats.py"
MAKE_ANNOT     = f"{config['ldsc_dir']}/make_annot.py"


# --- Study list -----------------------------------------------------

_study_df = pl.read_csv(STUDY_INFO)
STUDIES = _study_df["study"].to_list()
_STUDY_INFO = {row["study"]: row for row in _study_df.to_dicts()}

CHROMS = list(range(1, 23))
CONDITIONS = config["conditions"]


# --- Per-study reference selection ---------------------------------

def study_trait(study):
    return study.split("_")[0]


def _ref_key(study):
    info = _STUDY_INFO[study]
    ancestry = info["ancestry"]
    build = info["build"]
    if ancestry in ("EUR", "Multi"):
        prefix = "EUR"
    elif ancestry == "EAS":
        prefix = "EAS"
    else:
        raise ValueError(f"Unsupported ancestry for {study}: {ancestry}")
    key = f"{prefix}_{build}"
    return key


def study_ancestry(study):
    return "EAS" if _STUDY_INFO[study]["ancestry"] == "EAS" else "EUR"


def study_bim_prefix(study):
    return _ws(config["bim"][_ref_key(study)])


def study_gene_coord(study):
    return _ws(config["gene_coord"][_STUDY_INFO[study]["build"]])


def study_snp_list(study):
    return _ws(config["snp_list"][study_ancestry(study)])


def study_baseline_ld(study):
    return _ws(config["baseline_ld"][_ref_key(study)])


def study_weights(study):
    return _ws(config["weights"][_ref_key(study)])


def study_allele_frq(study):
    return _ws(config["allele_frq"][_ref_key(study)])


# --- Wildcard constraints ------------------------------------------

wildcard_constraints:
    study = r"[^/]+",
    cond  = "|".join(CONDITIONS),
    chrom = r"\d+",
