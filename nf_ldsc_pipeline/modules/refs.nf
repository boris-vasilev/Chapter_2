// Helpers shared by all processes.

// Absolute path from a path relative to params.workspace.
def ws(rel) { "${params.workspace}/${rel}" }

// Reference files for one row of the study table: Multi uses the EUR panel,
// as in rules/common.smk.
def refs(row) {
    if (!(row.ancestry in ['EUR', 'EAS', 'Multi']))
        error "Unsupported ancestry for ${row.study}: ${row.ancestry}"
    def anc = row.ancestry == 'EAS' ? 'EAS' : 'EUR'
    def key = "${anc}_${row.build}".toString()
    [ bim:        ws(params.bim[key]),
      snp_list:   ws(params.snp_list[anc]),
      gene_coord: ws(params.gene_coord[row.build]),
      baseline:   ws(params.baseline_ld[key]),
      weights:    ws(params.weights[key]),
      frq:        ws(params.allele_frq[key]) ]
}

// Hub-annotation quantiles from the comma-separated hub_annot_quantile, kept as
// given (0.90 stays 0.90), and their annotation names: 'hubs' for one quantile,
// 'hubs_q<quantile>' each for several (as make_gene_sets_coloc.py writes them).
def hub_quantiles() {
    params.hub_annot_quantile == null ? [] :
        params.hub_annot_quantile.toString().tokenize(',').collect { q -> q.trim() }
}

def hub_annots() {
    def qs = hub_quantiles()
    qs.size() == 1 ? ['hubs'] : qs.collect { q -> "hubs_q${q}".toString() }
}

// Output folder name: every setting that changes the gene sets.
def variant() {
    def q   = params.max_degs_quantile
    def hub = (q as BigDecimal) < 1 ? "hub${q}_" : ''
    def ann = params.hub_annot_quantile != null ? "_hubannot${hub_quantiles().join('-')}" : ''
    "${hub}deg${params.min_de_genes}_cells${params.min_cells}_ot${params.on_target ? 1 : 0}" +
    "_p${params.padj}_z${params.min_abs_z}_lfc${params.min_abs_logfc}_h4${params.h4_threshold}" +
    "_k${params.min_reg_k}-${params.min_tgt_k}_w${params.window_size}_snp${params.min_snps_per_annot}${ann}"
}