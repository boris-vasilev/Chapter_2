# Heritability rule reads from ld_scores/ when filtering is off,
# ld_scores_filtered/ when on (DAG branches via the include in
# Snakefile + this directory selection).

_LD_DIR = LD_SCORES_FILTERED_DIR if config["filter_low_snp_annots"] else LD_SCORES_DIR


rule calculate_heritability:
    """Stratified LD-score regression: --h2 with the combined annotation."""
    input:
        sumstats = f"{SUMSTATS_DIR}/{{study}}.sumstats.gz",
        ldscores = lambda w: expand(
            f"{_LD_DIR}/{w.study}/{w.cond}_combined.{{chrom}}.l2.ldscore.gz",
            chrom=CHROMS,
        ),
        m_files = lambda w: expand(
            f"{_LD_DIR}/{w.study}/{w.cond}_combined.{{chrom}}.l2.M_5_50",
            chrom=CHROMS,
        ),
    output:
        f"{LDSC_RESULTS_DIR}/{{study}}/{{cond}}_combined.results",
    params:
        baseline   = lambda w: study_baseline_ld(w.study),
        weights    = lambda w: study_weights(w.study),
        allele_frq = lambda w: study_allele_frq(w.study),
        ld_prefix  = lambda w: f"{_LD_DIR}/{w.study}/{w.cond}_combined.",
        out_prefix = lambda w: f"{LDSC_RESULTS_DIR}/{w.study}/{w.cond}_combined",
        ldsc       = LDSC,
    shell:
        """
        python {params.ldsc} \
            --h2 {input.sumstats} \
            --ref-ld-chr {params.baseline},{params.ld_prefix} \
            --w-ld-chr {params.weights} \
            --frqfile-chr {params.allele_frq} \
            --overlap-annot \
            --print-coefficients \
            --print-cov \
            --out {params.out_prefix}
        """
