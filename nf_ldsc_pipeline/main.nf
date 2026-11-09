// Perturb-seq gene sets for GWAS traits -> S-LDSC with baselineLD v2.2.
// Parameters are in nextflow.config; see CSD3_RUNBOOK.md to run it.

include { EXTRACT_DE; MAKE_GENE_SETS; MUNGE_SUMSTATS } from './modules/prepare'
include { LDSC_CHAIN } from './modules/ldsc_chain'
include { ws; variant; hub_quantiles } from './modules/refs'

workflow {
    def conditions = params.conditions.tokenize(',')
    def wanted = params.studies ? params.studies.tokenize(',') : null
    def outdir = "${params.outdir}/${variant()}"
    def study_info = file(ws(params.study_info))
    def hub_annot = params.hub_annot_quantile != null

    studies = channel.fromPath(study_info)
        .splitCsv(header: true)
        .filter { row -> !wanted || row.study in wanted }

    // params.de_results reuses an existing table and skips extraction. The hub
    // list needs extraction: hubs are ranked before any edge filtering.
    if (params.de_results && hub_annot)
        error "--hub_annot_quantile needs EXTRACT_DE; don't combine it with --de_results"
    if (hub_annot && (params.max_degs_quantile as BigDecimal) < 1)
        error "--hub_annot_quantile annotates hubs; don't combine it with hub exclusion (--max_degs_quantile < 1)"
    def hub_qs = hub_quantiles()
    if (hub_annot && !hub_qs.every { q -> q.isBigDecimal() && (q as BigDecimal) > 0 && (q as BigDecimal) < 1 })
        error "--hub_annot_quantile: every quantile must be a number in (0, 1), got '${params.hub_annot_quantile}'"
    if (hub_annot && (0..<hub_qs.size() - 1).any { i -> (hub_qs[i] as BigDecimal) >= (hub_qs[i + 1] as BigDecimal) })
        error "--hub_annot_quantile: quantiles must be strictly increasing, got '${params.hub_annot_quantile}'"
    if (params.de_results) {
        de = channel.value(file(ws(params.de_results)))
    } else {
        extract = EXTRACT_DE(file(ws(params.de_h5ad)))
        de = extract.de
    }
    // Hub list for the trait-agnostic annotation, or a placeholder when it is off.
    hubs = hub_annot ? extract.hubs.ifEmpty { error "EXTRACT_DE wrote no hubs.csv" }
                     : channel.value(file("${projectDir}/assets/NO_FILE"))

    MAKE_GENE_SETS(de, file(ws(params.coloc_dat)), studies.map { row -> row.study }.collect(),
                   conditions, hubs)
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