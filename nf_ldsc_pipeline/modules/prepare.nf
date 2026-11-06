// DE table, gene sets and munged sumstats: the inputs to LDSC_CHAIN.

include { ws; variant } from './refs'

// Filter the DE h5ad. In: h5ad. Out: DE_results.csv (one row per kept edge).
process EXTRACT_DE {
    publishDir "${params.outdir}/${variant()}", mode: 'copy'

    input:
    path h5ad

    output:
    path 'DE_results.csv'

    script:
    def on_target = params.on_target ? '--on-target' : '--no-on-target'
    """
    python ${moduleDir}/../bin/extract_de_results.py --h5ad ${h5ad} --out DE_results.csv \
        --padj-threshold ${params.padj} \
        --z-thresh ${params.min_abs_z} \
        --min-degs ${params.min_de_genes} \
        --max-degs-quantile ${params.max_degs_quantile} \
        --min-cells ${params.min_cells} ${on_target}
    """
}

// Regulator, target and anchor sets. In: DE table, coloc table, study and
// condition lists. Out: gene_sets/{study}_{cond}_{annot}.txt.
process MAKE_GENE_SETS {
    publishDir "${params.outdir}/${variant()}", mode: 'copy'

    input:
    path de
    path coloc
    val studies
    val conditions

    output:
    path 'gene_sets'

    script:
    """
    python ${moduleDir}/../bin/make_gene_sets_coloc.py \
        --de-results ${de} --coloc ${coloc} \
        --studies ${studies.collect { "'${it}'" }.join(' ')} \
        --conditions ${conditions.join(' ')} \
        --cell-type ${params.cell_type} --h4-threshold ${params.h4_threshold} \
        --min-anchors-per-regulator ${params.min_reg_k} \
        --min-perturbations-per-target ${params.min_tgt_k} \
        --min-abs-logfc ${params.min_abs_logfc} --h4-only --out-dir gene_sets
    """
}

// Munge one study's GWAS. In: study row, study table. Out: row, sumstats.
// storeDir reuses sumstats that already exist, so nothing is re-munged.
process MUNGE_SUMSTATS {
    tag "${row.study}"
    storeDir "${params.workspace}/${params.sumstats_dir}"

    input:
    val row
    path study_info

    output:
    tuple val(row), path("${row.study}.sumstats.gz")

    script:
    """
    python ${moduleDir}/../bin/run_munge_sumstats.py \
        --gwas '${ws(params.gwas_dir)}/${row.study}.tsv.gz' --study '${row.study}' \
        --study-info ${study_info} --out-prefix '${row.study}' \
        --munge ${params.ldsc_dir}/munge_sumstats.py
    """
}
