# Optional pass: drop annotation columns whose chr1 LD scores have
# variance == 0 or fewer than `min_snps_per_annot` total tagged SNPs,
# then recompute LD scores on the filtered set. Mirrors the
# 6_filter_annots.py + 7_calculate_heritability_filtered.sbatch flow.

rule filter_annot_files:
    """Filter low-SNP-count annot columns over chr 1..22, per (study, cond)."""
    input:
        ldscore_chr1 = f"{LD_SCORES_DIR}/{{study}}/{{cond}}_combined.1.l2.ldscore.gz",
        annots = expand(
            f"{LD_SCORES_DIR}/{{{{study}}}}/{{{{cond}}}}_combined.{{chrom}}.annot",
            chrom=CHROMS,
        ),
    output:
        annots = expand(
            f"{LD_SCORES_FILTERED_DIR}/{{{{study}}}}/{{{{cond}}}}_combined.{{chrom}}.annot",
            chrom=CHROMS,
        ),
    params:
        in_prefix  = lambda w: f"{LD_SCORES_DIR}/{w.study}/{w.cond}_combined",
        out_prefix = lambda w: f"{LD_SCORES_FILTERED_DIR}/{w.study}/{w.cond}_combined",
        min_snps   = config["min_snps_per_annot"],
    shell:
        """
        python scripts/filter_annots.py \
            --ldscore {input.ldscore_chr1} \
            --in-prefix {params.in_prefix} \
            --out-prefix {params.out_prefix} \
            --min-snps {params.min_snps}
        """


rule calculate_ld_scores_filtered:
    """LDSC --l2 over chr 1..22 in parallel, per (study, cond), filtered set."""
    input:
        annots = expand(
            f"{LD_SCORES_FILTERED_DIR}/{{{{study}}}}/{{{{cond}}}}_combined.{{chrom}}.annot",
            chrom=CHROMS,
        ),
    output:
        ldscore = expand(
            f"{LD_SCORES_FILTERED_DIR}/{{{{study}}}}/{{{{cond}}}}_combined.{{chrom}}.l2.ldscore.gz",
            chrom=CHROMS,
        ),
        m = expand(
            f"{LD_SCORES_FILTERED_DIR}/{{{{study}}}}/{{{{cond}}}}_combined.{{chrom}}.l2.M",
            chrom=CHROMS,
        ),
        m550 = expand(
            f"{LD_SCORES_FILTERED_DIR}/{{{{study}}}}/{{{{cond}}}}_combined.{{chrom}}.l2.M_5_50",
            chrom=CHROMS,
        ),
    params:
        bim_prefix = lambda w: study_bim_prefix(w.study),
        snp_list   = lambda w: study_snp_list(w.study),
        prefix     = lambda w: f"{LD_SCORES_FILTERED_DIR}/{w.study}/{w.cond}_combined",
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
