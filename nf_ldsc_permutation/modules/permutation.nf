// The permutation-specific steps. Everything LDSC is LDSC_CHAIN from the main pipeline.

include { ws; refs } from '../../nf_ldsc_pipeline/modules/refs'

// Per-condition screen baseMean, the expression axis of the matching.
// In: DE h5ad. Out: screen_expression.tsv. storeDir: an existing table is reused.
process BUILD_EXPRESSION {
    storeDir "${ws(params.expression_dir)}"

    input:
    path h5ad

    output:
    path 'screen_expression.tsv'

    script:
    """
    python ${moduleDir}/../bin/build_screen_expression.py --h5ad ${h5ad} --out screen_expression.tsv
    """
}

// Gene sets for real and perms 1..n_perm of one study x condition, one call per perm,
// so seeds are perm_rng(seed, perm, cond) as in Stage 2. Takes seconds per perm.
// In: row (study table row plus cond, outdir), DE table, expression, coloc table.
// Out: row, {study}/perm{P}/gene_sets/* (the three sets and _matching.tsv).
// Files, not directories: a published directory would replace the one another
// condition of the same study has already published.
process PERM_GENE_SETS {
    tag "${row.study} ${row.cond}"
    publishDir path: { row.outdir }, mode: 'copy'
    // A real set that differs from reference_gene_sets stops the whole run.
    errorStrategy 'terminate'

    input:
    tuple val(row), path(de), path(expression), path(coloc)

    output:
    tuple val(row), path("${row.study}/perm*/gene_sets/*")

    script:
    def ref = params.reference_gene_sets ? "--reference-dir ${ws(params.reference_gene_sets)}" : ''
    """
    gene_sets() {
        python ${moduleDir}/../bin/make_permutation_gene_sets.py \
            --de-results ${de} --coloc ${coloc} --expression ${expression} \
            --gene-coord ${refs(row).gene_coord} \
            --protein-coding-coord ${ws(params.protein_coding_coord[row.build])} \
            --study '${row.study}' --condition ${row.cond} \
            --cell-type ${params.cell_type} --h4-threshold ${params.h4_threshold} \
            --min-anchors-per-regulator ${params.min_reg_k} \
            --min-perturbations-per-target ${params.min_tgt_k} \
            --min-abs-logfc ${params.min_abs_logfc} \
            --n-degree-bins ${params.n_degree_bins} --n-expr-bins ${params.n_expr_bins} \
            --knn-k ${params.knn_k} --base-seed ${params.seed} ${ref} \
            --perm \$1 --out-dir '${row.study}'/perm\$1/gene_sets
    }
    gene_sets real
    for ((p = 1; p <= ${params.n_perm}; p++)); do gene_sets \$p; done
    """
}

// Per-run and pooled summaries plus null histograms.
// In: manifest of "relative path <tab> work path" lines, rebuilt here as a
// {study}/perm{P}/... tree (published copies may still be in flight); run list.
// Out: {study}_{cond}/ per run and the cross-run CSVs.
process COLLECT {
    publishDir path: { outdir }, mode: 'copy'

    input:
    path manifest
    val runs
    val outdir

    output:
    path '*.csv', optional: true
    path '*_*', type: 'dir', optional: true

    script:
    """
    while IFS=\$'\\t' read -r rel src; do
        mkdir -p "tree/\$(dirname "\$rel")"
        ln -s "\$src" "tree/\$rel"
    done < ${manifest}
    python ${moduleDir}/../bin/collect_permutation_ldsc.py --root tree --out-dir . --runs ${runs}
    """
}
