rule munge_sumstats:
    """Run LDSC munge_sumstats.py on a single GWAS file."""
    input:
        gwas       = lambda w: f"{GWAS_DIR}/{w.study}.tsv.gz",
        study_info = STUDY_INFO,
    output:
        f"{SUMSTATS_DIR}/{{study}}.sumstats.gz",
    params:
        out_prefix = lambda w: f"{SUMSTATS_DIR}/{w.study}",
        munge      = MUNGE_SUMSTATS,
    shell:
        """
        python scripts/run_munge_sumstats.py \
            --gwas {input.gwas} \
            --study {wildcards.study} \
            --study-info {input.study_info} \
            --out-prefix {params.out_prefix} \
            --munge {params.munge}
        """
