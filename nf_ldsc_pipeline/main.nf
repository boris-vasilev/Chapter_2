// Perturb-seq gene sets for GWAS traits -> S-LDSC with baselineLD v2.2.
// Parameters are in nextflow.config; see CSD3_RUNBOOK.md to run it.

include { EXTRACT_DE; MAKE_GENE_SETS; MUNGE_SUMSTATS } from './modules/prepare'
include { LDSC_CHAIN } from './modules/ldsc_chain'
include { ws; variant } from './modules/refs'

workflow {
    def conditions = params.conditions.tokenize(',')
    def wanted = params.studies ? params.studies.tokenize(',') : null
    def outdir = "${params.outdir}/${variant()}"
    def study_info = file(ws(params.study_info))

    studies = channel.fromPath(study_info)
        .splitCsv(header: true)
        .filter { row -> !wanted || row.study in wanted }

    // params.de_results reuses an existing table and skips extraction.
    de = params.de_results ? channel.value(file(ws(params.de_results)))
                           : EXTRACT_DE(file(ws(params.de_h5ad)))

    MAKE_GENE_SETS(de, file(ws(params.coloc_dat)), studies.map { row -> row.study }.collect(), conditions)
    MUNGE_SUMSTATS(studies, study_info)

    // One chain per study x condition. Paths in meta decide the published layout.
    chains = MUNGE_SUMSTATS.out
        .combine(MAKE_GENE_SETS.out)
        .combine(channel.fromList(conditions))
        .map { row, sumstats, gene_sets, cond ->
            def meta = row + [cond: cond, outdir: outdir,
                              res_dir: "ldsc_results/${row.study}".toString(),
                              ld_dir: "ld_scores/${row.study}".toString()]
            tuple(meta, gene_sets, sumstats)
        }
    LDSC_CHAIN(chains)
}
