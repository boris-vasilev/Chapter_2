rule extract_de_results:
    """Extract significant DE stats from the Perturb-seq h5ad."""
    input:
        h5ad = DE_INPUT_H5AD,
    output:
        de = DE_RESULTS,
    params:
        padj = config["padj_threshold"],
    shell:
        """
        python scripts/extract_de_results.py \
            --h5ad {input.h5ad} \
            --out {output.de} \
            --padj-threshold {params.padj}
        """
