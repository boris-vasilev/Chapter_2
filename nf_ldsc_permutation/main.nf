// Matched-null S-LDSC permutation: for each study x condition in params.runs, the real
// coloc anchors and n_perm matched pseudo-anchor sets, each through the main
// pipeline's LDSC_CHAIN. Parameters are in nextflow.config; see README.md.

include { EXTRACT_DE } from '../nf_ldsc_pipeline/modules/prepare'
include { LDSC_CHAIN } from '../nf_ldsc_pipeline/modules/ldsc_chain'
include { ws; variant } from '../nf_ldsc_pipeline/modules/refs'
include { BUILD_EXPRESSION; PERM_GENE_SETS; COLLECT } from './modules/permutation'

// Last n components of a path, e.g. {study}/perm{P}/gene_sets/{file}.
def tail(p, n) { p.subpath(p.nameCount - n, p.nameCount).toString() }

workflow {
    // Edge filters (variant) plus matching options name the output folder.
    def tag = "${variant()}_bins${params.n_degree_bins}-${params.n_expr_bins}" +
              "_knn${params.knn_k}_seed${params.seed}"
    def root = "${ws(params.perm_root)}/${tag}".toString()
    if (root.contains('/permutations/heritability/'))
        error "Stage 2 outputs in data/permutations/heritability/ are read-only"

    // --studies and --conditions (main pipeline parameters) select from params.runs.
    def wanted = params.studies ? params.studies.tokenize(',') : params.runs.keySet()
    def conds = params.conditions.tokenize(',')
    def runs = params.runs.findAll { s, cs -> s in wanted }
        .collectEntries { s, cs -> [s, cs.findAll { it in conds }] }.findAll { s, cs -> cs }
    def rows = file(ws(params.study_info)).splitCsv(header: true).findAll { it.study in runs }
    def missing = (wanted - params.runs.keySet()) + (runs.keySet() - rows*.study)
    if (missing) error "Not in params.runs or ${params.study_info}: ${missing.join(', ')}"

    // params.de_results reuses an existing table and skips extraction (which
    // publishes to the main pipeline's outdir, where it can be reused later).
    de = params.de_results ? channel.value(file(ws(params.de_results)))
                           : EXTRACT_DE(file(ws(params.de_h5ad)))
    expression = BUILD_EXPRESSION(file(ws(params.de_h5ad)))

    jobs = channel.fromList(rows.collectMany { r ->
            runs[r.study].collect { c -> r + [cond: c, outdir: root] } })
        .combine(de).combine(expression)
        .map { row, d, e -> tuple(row, d, e, file(ws(params.coloc_dat))) }
    PERM_GENE_SETS(jobs)
    if (params.gene_sets_only) return

    // toList() waits for every gene-set job, so no chain starts before all real
    // sets have passed the reference check. Then one chain per perm, given the
    // perm's gene_sets directory in the task's work directory.
    chains = PERM_GENE_SETS.out.toList().flatMap { done ->
        done.collectMany { run ->
            def (row, files) = run
            [files].flatten()*.parent.unique().collect { d ->
                def perm = "${row.study}/${d.parent.name}"
                def meta = row + [res_dir: "${perm}/ldsc_results".toString(),
                                  ld_dir: "${perm}/ld_scores/${row.cond}".toString()]
                def sumstats = file(ws("${params.sumstats_dir}/${row.study}.sumstats.gz"), checkIfExists: true)
                tuple(meta, d, sumstats)
            }
        }
    }
    LDSC_CHAIN(chains)

    // Collect from work-directory paths laid out as the published tree.
    gene_set_files = PERM_GENE_SETS.out.map { it[1] }.flatten().map { f -> [tail(f, 4), f] }
    chain_files = LDSC_CHAIN.out[1].flatten().filter { it.name.endsWith('.results') }
        .mix(LDSC_CHAIN.out[2])
        .map { f -> [tail(f, 4), f] }
    manifest = gene_set_files.mix(chain_files)
        .map { rel, f -> "${rel}\t${f}\n" }
        .collectFile(name: 'manifest.tsv', sort: true)
    COLLECT(manifest, runs.collect { s, cs -> "'${s}:${cs.join(',')}'" }.join(' '), root)
}
