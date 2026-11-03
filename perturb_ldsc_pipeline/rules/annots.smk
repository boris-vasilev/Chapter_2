rule make_annot:
    """Build a binary SNP annotation file for one (study, cond, annot) — loops chr 1..22."""
    input:
        gene_set = lambda w: f"{GENE_SETS_DIR}/{w.study}_{w.cond}_{w.annot}.txt",
    output:
        expand(
            f"{ANNOTATIONS_DIR}/{{{{study}}}}/{{{{cond}}}}_{{{{annot}}}}.{{chrom}}.annot.gz",
            chrom=CHROMS,
        ),
    params:
        bim_prefix = lambda w: study_bim_prefix(w.study),
        gene_coord = lambda w: study_gene_coord(w.study),
        out_prefix = lambda w: f"{ANNOTATIONS_DIR}/{w.study}/{w.cond}_{w.annot}",
        window     = config["window_size"],
        make_annot = MAKE_ANNOT,
    wildcard_constraints:
        annot = "|".join(ANNOT_TYPES),
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
        """


rule combine_annots:
    """Combine per-annot files into one combined .annot per (study, cond) — loops chr 1..22."""
    input:
        annots = lambda w: [
            f"{ANNOTATIONS_DIR}/{w.study}/{w.cond}_{a}.{c}.annot.gz"
            for a in ANNOT_TYPES
            for c in CHROMS
        ],
    output:
        expand(
            f"{LD_SCORES_DIR}/{{{{study}}}}/{{{{cond}}}}_combined.{{chrom}}.annot",
            chrom=CHROMS,
        ),
    params:
        in_dir    = lambda w: f"{ANNOTATIONS_DIR}/{w.study}",
        out_dir   = lambda w: f"{LD_SCORES_DIR}/{w.study}",
        condition = lambda w: w.cond,
        names     = " ".join(ANNOT_TYPES),
        min_snps  = config.get("min_snps_per_annot", 1000),
    shell:
        """
        python scripts/combine_annots.py \
            --in-dir {params.in_dir} \
            --out-dir {params.out_dir} \
            --condition {params.condition} \
            --names {params.names} \
            --min-snps {params.min_snps}
        """
