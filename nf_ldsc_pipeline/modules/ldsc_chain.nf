// One chain: make_annot x3 -> combine -> LD scores -> S-LDSC h2.
// One process because each task is one SLURM job (see PLAN.md, section 7).
// Knows nothing about real vs null runs: `meta` says what to run and where
// to publish, so the permutation wrapper reuses it unchanged.
//
// In:  meta (study table row plus cond, outdir, res_dir, ld_dir), gene-set
//      directory, sumstats.
// Out: {res_dir}/{cond}_combined.{log,results,cov}, combine.txt (sets kept),
//      and LD scores under {ld_dir} when params.keep_ld is set.

include { refs } from './refs'

process LDSC_CHAIN {
    tag "${meta.study} ${meta.cond}"
    label 'chain'
    publishDir path: { meta.outdir }, mode: 'copy'
    // Exit 3 = every set under min_snps: an expected non-result, no files published.
    errorStrategy { task.exitStatus == 3 ? 'ignore' : task.exitStatus in [137, 140, 143] ? 'retry' : 'finish' }
    maxRetries 2

    input:
    tuple val(meta), path(gene_sets), path(sumstats)

    output:
    path "${meta.res_dir}/${meta.cond}_combined.log"
    path "${meta.res_dir}/${meta.cond}_combined.{results,cov}", optional: true
    path "${meta.res_dir}/${meta.cond}_combined.combine.txt"
    path "${meta.ld_dir}/*", optional: true

    script:
    def r = refs(meta)
    def bin = "${moduleDir}/../bin"
    def annots = 'regulators_H4 targets_H4 H4_genes'
    def gene_files = annots.tokenize().collect { "'${gene_sets}/${meta.study}_${meta.cond}_${it}.txt'" }.join(' ')
    def ld = "${meta.ld_dir}/${meta.cond}_combined"
    def out = "${meta.res_dir}/${meta.cond}_combined"
    """
    mkdir -p annot '${meta.ld_dir}' '${meta.res_dir}'
    for a in ${annots}; do
        python ${bin}/run_make_annot.py --make-annot ${params.ldsc_dir}/make_annot.py \
            --gene-set-file '${gene_sets}/${meta.study}_${meta.cond}_'\$a.txt \
            --gene-coord-file ${r.gene_coord} --windowsize ${params.window_size} \
            --bim-prefix ${r.bim} --out-prefix annot/${meta.cond}_\$a \
            --max-parallel ${task.cpus}
    done

    python ${bin}/combine_annots.py --in-dir annot --out-dir '${meta.ld_dir}' \
        --condition ${meta.cond} --names ${annots} \
        --min-snps ${params.min_snps_per_annot} \
        --gene-sets ${gene_files} \
        --min-genes ${params.min_genes_per_annot} > '${out}.combine.txt' 2> combine.err || {
        cat combine.err >&2
        grep -q 'refusing to write' combine.err && exit 3
        exit 1
    }

    python ${bin}/run_ldsc_l2.py --ldsc ${params.ldsc_dir}/ldsc.py --bim-prefix ${r.bim} \
        --snp-list ${r.snp_list} --prefix '${ld}' --max-parallel ${task.cpus}

    # A failed fit keeps its .log, which chapter2/vis reports as "run failed".
    python ${params.ldsc_dir}/ldsc.py --h2 '${sumstats}' \
        --ref-ld-chr ${r.baseline},'${ld}.' --w-ld-chr ${r.weights} --frqfile-chr ${r.frq} \
        --overlap-annot --print-coefficients --print-cov --out '${out}' || true

    rm -rf annot ${params.keep_ld ? '' : "'${meta.ld_dir}'"}
    """
}
