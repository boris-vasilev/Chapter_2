# Gene-set module: coloc.
#
# A "module" defines:
#   * ANNOT_TYPES — annotation columns that downstream rules combine
#   * rule make_gene_sets — produces {GENE_SETS_DIR}/{study}_{cond}_{annot}.txt
#
# Keying on study identifier (= imd column in coloc_dat.tsv) so the
# coloc filter actually matches. The original 2_create_gene_sets.py
# behaved the same way — it pulled the trait list from H4['imd'].unique().
#
# To swap in another module (e.g. L2G from OpenTargets), create
# rules/gene_sets_l2g.smk providing the same outputs and set
# config['gene_set_module'] = 'l2g'.

# Read params from the *active* module's config block — this lets
# coloc and coloc_strict share rules/gene_sets_coloc.smk but use
# different filter thresholds (e.g. min_anchors_per_regulator).
_coloc_cfg = mod_cfg

if _coloc_cfg["h4_only"]:
    ANNOT_TYPES = ["regulators_H4", "targets_H4", "H4_genes"]
else:
    ANNOT_TYPES = [
        "regulators_H4",
        "regulators_H3",
        "targets_H4",
        "targets_H3",
        "H4_genes",
        "H3_genes",
    ]


rule make_gene_sets:
    input:
        de    = DE_RESULTS,
        coloc = COLOC_DAT,
    output:
        expand(
            f"{GENE_SETS_DIR}/{{study}}_{{cond}}_{{annot}}.txt",
            study=STUDIES,
            cond=CONDITIONS,
            annot=ANNOT_TYPES,
        ),
    params:
        out_dir    = GENE_SETS_DIR,
        studies    = " ".join(STUDIES),
        conditions = " ".join(CONDITIONS),
        cell_type  = _coloc_cfg["cell_type"],
        h3         = _coloc_cfg["h3_threshold"],
        h4         = _coloc_cfg["h4_threshold"],
        h4_only    = "--h4-only" if _coloc_cfg["h4_only"] else "",
        # Accept either short (min_reg_k / min_tgt_k) or long names in config.
        # Default = 1 (no filtering — every regulator/target counts). The
        # strict variant must opt in by setting min_reg_k / min_tgt_k > 1.
        min_reg_k  = _coloc_cfg.get("min_reg_k",
                                    _coloc_cfg.get("min_anchors_per_regulator", 1)),
        min_tgt_k  = _coloc_cfg.get("min_tgt_k",
                                    _coloc_cfg.get("min_perturbations_per_target", 1)),
    shell:
        """
        python scripts/make_gene_sets_coloc.py \
            --de-results {input.de} \
            --coloc {input.coloc} \
            --studies {params.studies} \
            --conditions {params.conditions} \
            --cell-type {params.cell_type} \
            --h3-threshold {params.h3} \
            --h4-threshold {params.h4} \
            --min-anchors-per-regulator {params.min_reg_k} \
            --min-perturbations-per-target {params.min_tgt_k} \
            {params.h4_only} \
            --out-dir {params.out_dir}
        """
