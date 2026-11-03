#!/usr/bin/env python
"""Extract significant DE stats from a Perturb-seq DE h5ad.

Keeps perturbations passing the on-target, DEG-count and cell-count filters
and writes a long-format CSV of gene-pair DE stats passing the downstream
filters. Edges are selected on the matrices before any long table is built,
so memory scales with kept edges, not perturbations x genes.
"""

import argparse

import numpy as np
import polars as pl
import scanpy as sc


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--padj-threshold", type=float, default=0.05)
    parser.add_argument(
        "--min-degs",
        type=int,
        default=75,
        help="Minimum n_total_de_genes for strong perturbation "
        "(default 75, matching Zhu et al.).",
    )

    parser.add_argument(
        "--min-cells",
        type=int,
        default=50,
        help="Minimum n_cells_target (default 50, Zhu et al.).",
    )

    parser.add_argument(
        "--z-thresh",
        type=float,
        default=2.0,
        help="|zscore| threshold for downstream gene " "(default 2.0).",
    )

    parser.add_argument(
        "--on-target",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Require on-target KD and a significant on-target effect.",
    )

    args = parser.parse_args()

    adata = sc.read_h5ad(args.h5ad)
    mask = np.ones(adata.n_obs, dtype=bool)
    if args.on_target:
        mask &= (adata.obs["ontarget_effect_category"] == "on-target KD").to_numpy()
        print(f"  on-target KD obs : {mask.sum()} / {len(mask)}")
        mask &= adata.obs["ontarget_significant"].astype(bool).to_numpy()

    # strong perturbation filters
    mask &= (adata.obs["n_total_de_genes"] >= args.min_degs).to_numpy()
    mask &= (adata.obs["n_cells_target"] >= args.min_cells).to_numpy()

    print(f"  after strong-perturbation filter: {mask.sum()} obs")

    on_target = adata[mask].copy()
    del adata

    obs_df = pl.from_pandas(
        on_target.obs[
            [
                "target_contrast",
                "target_contrast_gene_name",
                "culture_condition",
                "ontarget_effect_size",
                "target_baseMean",
                "n_total_de_genes",
                "n_cells_target",
            ]
        ].reset_index(names="obs_id")
    ).with_columns(pl.col("target_contrast_gene_name").cast(pl.String))

    var_df = pl.from_pandas(
        on_target.var[["gene_ids", "gene_name"]].reset_index(names="var_id")
    )

    # Strong downstream effect filter: padj < threshold AND |zscore| >= z_thresh.
    # ~(< z) keeps NaN z-scores, as the previous polars filter did.
    print("Applying downstream effect filters …")
    keep = on_target.layers["adj_p_value"] < args.padj_threshold
    keep &= ~(np.abs(on_target.layers["zscore"]) < args.z_thresh)
    obs_idx, var_idx = np.nonzero(keep)  # row-major, the previous row order
    obs_idx, var_idx = pl.Series(obs_idx), pl.Series(var_idx)

    def layer(name):
        return on_target.layers[name][keep]

    # Exclude on-target row (cis_ensg == trans_ensg) from targets
    de_filtered = pl.DataFrame(
        {
            "padj": layer("adj_p_value"),
            "base_DESeq2": layer("baseMean"),
            "lfcSE": layer("lfcSE"),
            "logFC": layer("log_fc"),
            "zscore": layer("zscore"),
        }
    ).with_columns(
        [
            obs_df["culture_condition"].gather(obs_idx).alias("culture_condition"),
            obs_df["target_contrast"].gather(obs_idx).alias("cis_ensg"),
            obs_df["target_contrast_gene_name"].gather(obs_idx).alias("cis_gene"),
            obs_df["ontarget_effect_size"].gather(obs_idx).alias("perturb_eff"),
            obs_df["target_baseMean"].gather(obs_idx).alias("cis_baseMean"),
            obs_df["n_total_de_genes"].gather(obs_idx).alias("n_total_de_genes"),
            obs_df["n_cells_target"].gather(obs_idx).alias("n_cells_target"),
            var_df["gene_ids"].gather(var_idx).alias("trans_ensg"),
            var_df["gene_name"].gather(var_idx).alias("trans_gene"),
        ]
    ).filter(pl.col("cis_ensg") != pl.col("trans_ensg"))
    print(f"  rows after downstream filter: {len(de_filtered):,}")

    de_filtered.write_csv(args.out)


if __name__ == "__main__":
    main()
