# Gene-set module: l2g (OpenTargets Locus-to-Gene).
#
# Each row in GWAS_study_info.csv's `study` column is the OpenTargets
# studyId; no separate mapping is needed. The script joins l2g →
# credible_set on studyLocusId to get geneIds per study.

ANNOT_TYPES = [
    "regulators_l2g",
    "targets_l2g",
    "l2g_genes",
]


rule make_gene_sets:
    input:
        de = DE_RESULTS,
    output:
        expand(
            f"{GENE_SETS_DIR}/{{study}}_{{cond}}_{{annot}}.txt",
            study=STUDIES,
            cond=CONDITIONS,
            annot=ANNOT_TYPES,
        ),
    params:
        out_dir         = GENE_SETS_DIR,
        opentargets_dir = L2G_DIR,
        studies         = " ".join(STUDIES),
        conditions      = " ".join(CONDITIONS),
        score_threshold = config["l2g"]["score_threshold"],
    shell:
        """
        python scripts/make_gene_sets_l2g.py \
            --de-results {input.de} \
            --opentargets-dir {params.opentargets_dir} \
            --studies {params.studies} \
            --conditions {params.conditions} \
            --score-threshold {params.score_threshold} \
            --out-dir {params.out_dir}
        """
