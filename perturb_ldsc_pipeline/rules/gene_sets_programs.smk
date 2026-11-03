# Gene-set module: Perturb-seq programmes (Zhu clusters), cell-type-specific
# (per-annotation) S-LDSC.
#
# This module is structurally different from coloc / l2g and therefore ships
# its OWN annot -> ld-score -> heritability backend instead of reusing the
# generic rules in annots.smk / ld_scores.smk / heritability.smk:
#
#   * Programme gene sets are STUDY-INDEPENDENT (a programme is the same gene
#     list for every GWAS trait). LD scores therefore only depend on
#     (programme, condition, reference-panel) and are computed ONCE per
#     reference panel and reused across all studies sharing that panel.
#   * Each programme is tested on its own (Finucane 2018 cell-type-specific
#     design): --ref-ld-chr baseline,<one programme>. The programme's
#     Coefficient_z-score is the enrichment signal.
#
# Outputs collected into {LDSC_RESULTS_DIR}/program_enrichment.csv.

from collections import defaultdict

import polars as pl

_prog_cfg = mod_cfg
PROGRAM_ASSIGNMENTS = _ws(_prog_cfg["program_assignments"])
MIN_PROGRAM_SIZE = int(_prog_cfg.get("min_program_size", 10))

# --- enumerate (condition, cluster) targets at parse time ---------------
_pa = (
    pl.read_csv(PROGRAM_ASSIGNMENTS)
    .filter(pl.col("condition").is_in(CONDITIONS))
    .unique(subset=["condition", "cluster", "perturb_ensg"])
)
_sizes = _pa.group_by(["condition", "cluster"]).len()
PROGRAMS = sorted(
    (r["condition"], int(r["cluster"]))
    for r in _sizes.filter(pl.col("len") >= MIN_PROGRAM_SIZE).iter_rows(named=True)
)
if not PROGRAMS:
    raise ValueError(
        f"No programmes with >= {MIN_PROGRAM_SIZE} genes in {PROGRAM_ASSIGNMENTS}."
    )

# Reference panels actually needed by the study list (e.g. EUR_GRCh37).
REF_KEYS = sorted({_ref_key(s) for s in STUDIES})


# --- reference-panel-keyed lookups (vs the study-keyed ones in common) ---

def _refkey_build(rk):
    return rk.split("_", 1)[1]


def _refkey_ancestry(rk):
    return "EAS" if rk.startswith("EAS") else "EUR"


def refkey_bim(rk):
    return _ws(config["bim"][rk])


def refkey_gene_coord(rk):
    return _ws(config["gene_coord"][_refkey_build(rk)])


def refkey_snp_list(rk):
    return _ws(config["snp_list"][_refkey_ancestry(rk)])


def refkey_baseline(rk):
    return _ws(config["baseline_ld"][rk])


def refkey_weights(rk):
    return _ws(config["weights"][rk])


def refkey_allele_frq(rk):
    return _ws(config["allele_frq"][rk])


PROGRAMS_INDEX = f"{GENE_SETS_DIR}/programs_index.csv"

# Final pipeline target for this module (referenced by Snakefile's rule all).
PROGRAM_FINAL = f"{LDSC_RESULTS_DIR}/program_enrichment.csv"


wildcard_constraints:
    cluster = r"\d+",
    ref_key = "|".join(REF_KEYS),


rule make_gene_sets:
    """All programme gene sets + index in one shot (study-independent)."""
    input:
        programs = PROGRAM_ASSIGNMENTS,
    output:
        gene_sets = [f"{GENE_SETS_DIR}/{c}_program{cl}.txt" for c, cl in PROGRAMS],
        index = PROGRAMS_INDEX,
    params:
        out_dir    = GENE_SETS_DIR,
        conditions = " ".join(CONDITIONS),
        min_size   = MIN_PROGRAM_SIZE,
    shell:
        """
        python scripts/make_gene_sets_programs.py \
            --programs {input.programs} \
            --conditions {params.conditions} \
            --min-program-size {params.min_size} \
            --out-dir {params.out_dir}
        """


rule program_make_annot:
    """Binary SNP annotation for one (ref_panel, cond, programme), chr 1..22.

    Keyed by reference panel, not study, so it is built once per panel.
    run_make_annot.py writes .annot.gz; we gunzip so --l2 reads .annot.
    """
    input:
        gene_set = f"{GENE_SETS_DIR}/{{cond}}_program{{cluster}}.txt",
    output:
        expand(
            f"{LD_SCORES_DIR}/{{{{ref_key}}}}/{{{{cond}}}}_program{{{{cluster}}}}.{{chrom}}.annot",
            chrom=CHROMS,
        ),
    params:
        bim_prefix = lambda w: refkey_bim(w.ref_key),
        gene_coord = lambda w: refkey_gene_coord(w.ref_key),
        out_prefix = lambda w: f"{LD_SCORES_DIR}/{w.ref_key}/{w.cond}_program{w.cluster}",
        window     = config["window_size"],
        make_annot = MAKE_ANNOT,
    threads: 22
    shell:
        """
        mkdir -p $(dirname {params.out_prefix})
        python scripts/run_make_annot.py \
            --make-annot {params.make_annot} \
            --gene-set-file {input.gene_set} \
            --gene-coord-file {params.gene_coord} \
            --windowsize {params.window} \
            --bim-prefix {params.bim_prefix} \
            --out-prefix {params.out_prefix}
        for c in $(seq 1 22); do
            gunzip -f {params.out_prefix}.$c.annot.gz
        done
        """


rule program_ld_scores:
    """LD scores for one (ref_panel, cond, programme), chr 1..22."""
    input:
        annots = expand(
            f"{LD_SCORES_DIR}/{{{{ref_key}}}}/{{{{cond}}}}_program{{{{cluster}}}}.{{chrom}}.annot",
            chrom=CHROMS,
        ),
    output:
        ldscore = expand(
            f"{LD_SCORES_DIR}/{{{{ref_key}}}}/{{{{cond}}}}_program{{{{cluster}}}}.{{chrom}}.l2.ldscore.gz",
            chrom=CHROMS,
        ),
        m = expand(
            f"{LD_SCORES_DIR}/{{{{ref_key}}}}/{{{{cond}}}}_program{{{{cluster}}}}.{{chrom}}.l2.M",
            chrom=CHROMS,
        ),
        m550 = expand(
            f"{LD_SCORES_DIR}/{{{{ref_key}}}}/{{{{cond}}}}_program{{{{cluster}}}}.{{chrom}}.l2.M_5_50",
            chrom=CHROMS,
        ),
    params:
        bim_prefix = lambda w: refkey_bim(w.ref_key),
        snp_list   = lambda w: refkey_snp_list(w.ref_key),
        prefix     = lambda w: f"{LD_SCORES_DIR}/{w.ref_key}/{w.cond}_program{w.cluster}",
        ldsc       = LDSC,
    threads: 22
    shell:
        """
        python scripts/run_ldsc_l2.py \
            --ldsc {params.ldsc} \
            --bim-prefix {params.bim_prefix} \
            --snp-list {params.snp_list} \
            --prefix {params.prefix}
        """


rule program_heritability:
    """Cell-type-specific S-LDSC: baseline + one programme, per (study, cond)."""
    input:
        sumstats = f"{SUMSTATS_DIR}/{{study}}.sumstats.gz",
        ldscores = lambda w: expand(
            f"{LD_SCORES_DIR}/{_ref_key(w.study)}/{w.cond}_program{w.cluster}.{{chrom}}.l2.ldscore.gz",
            chrom=CHROMS,
        ),
        m_files = lambda w: expand(
            f"{LD_SCORES_DIR}/{_ref_key(w.study)}/{w.cond}_program{w.cluster}.{{chrom}}.l2.M_5_50",
            chrom=CHROMS,
        ),
    output:
        f"{LDSC_RESULTS_DIR}/{{study}}/{{cond}}_program{{cluster}}.results",
    params:
        baseline   = lambda w: refkey_baseline(_ref_key(w.study)),
        weights    = lambda w: refkey_weights(_ref_key(w.study)),
        allele_frq = lambda w: refkey_allele_frq(_ref_key(w.study)),
        ld_prefix  = lambda w: f"{LD_SCORES_DIR}/{_ref_key(w.study)}/{w.cond}_program{w.cluster}.",
        out_prefix = lambda w: f"{LDSC_RESULTS_DIR}/{w.study}/{w.cond}_program{w.cluster}",
        ldsc       = LDSC,
    shell:
        """
        # Force UTF-8 so studies with non-ASCII names (e.g. SCZ_Pardiñas)
        # don't trip LDSC's path handling under a C/POSIX locale.
        export LC_ALL=C.UTF-8 LANG=C.UTF-8
        mkdir -p $(dirname {params.out_prefix})
        python {params.ldsc} \
            --h2 {input.sumstats} \
            --ref-ld-chr {params.baseline},{params.ld_prefix} \
            --w-ld-chr {params.weights} \
            --frqfile-chr {params.allele_frq} \
            --overlap-annot \
            --print-coefficients \
            --out {params.out_prefix}
        """


rule collect_program_results:
    """Gather every programme's coefficient into one tidy CSV."""
    input:
        results = [
            f"{LDSC_RESULTS_DIR}/{s}/{c}_program{cl}.results"
            for s in STUDIES
            for c, cl in PROGRAMS
        ],
        index = PROGRAMS_INDEX,
    output:
        PROGRAM_FINAL,
    params:
        results_dir = LDSC_RESULTS_DIR,
        study_info  = STUDY_INFO,
    shell:
        """
        python scripts/collect_program_results.py \
            --results-dir {params.results_dir} \
            --study-info {params.study_info} \
            --programs-index {input.index} \
            --out {output}
        """
