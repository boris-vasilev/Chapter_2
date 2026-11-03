rule calculate_ld_scores:
    """LDSC --l2 over chr 1..22 in parallel, per (study, cond)."""
    input:
        annots = expand(
            f"{LD_SCORES_DIR}/{{{{study}}}}/{{{{cond}}}}_combined.{{chrom}}.annot",
            chrom=CHROMS,
        ),
    output:
        ldscore = expand(
            f"{LD_SCORES_DIR}/{{{{study}}}}/{{{{cond}}}}_combined.{{chrom}}.l2.ldscore.gz",
            chrom=CHROMS,
        ),
        m = expand(
            f"{LD_SCORES_DIR}/{{{{study}}}}/{{{{cond}}}}_combined.{{chrom}}.l2.M",
            chrom=CHROMS,
        ),
        m550 = expand(
            f"{LD_SCORES_DIR}/{{{{study}}}}/{{{{cond}}}}_combined.{{chrom}}.l2.M_5_50",
            chrom=CHROMS,
        ),
    params:
        bim_prefix = lambda w: study_bim_prefix(w.study),
        snp_list   = lambda w: study_snp_list(w.study),
        prefix     = lambda w: f"{LD_SCORES_DIR}/{w.study}/{w.cond}_combined",
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
