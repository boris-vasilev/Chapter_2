#!/usr/bin/env python
"""Extract significant DE stats from a Perturb-seq DE h5ad.

Keeps gRNA observations whose ``ontarget_effect_category`` is
``"on-target KD"`` and writes a long-format CSV of gene-pair DE stats
with rows passing ``padj < --padj-threshold``.
"""

import argparse

import numpy as np
import polars as pl
import scanpy as sc


def layer_to_df(adata, layer_name, value_name):
    mat = adata.layers[layer_name]
    return pl.DataFrame(
        {
            "obs_id": np.repeat(adata.obs_names, mat.shape[1]),
            "var_id": np.tile(adata.var_names, mat.shape[0]),
            value_name: mat.ravel(order="C"),
        }
    )


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

    args = parser.parse_args()

    adata = sc.read_h5ad(args.h5ad)
    # on-target KD only
    mask = adata.obs["ontarget_effect_category"] == "on-target KD"
    print(f"  on-target KD obs : {mask.sum()} / {len(mask)}")

    # strong perturbation filters
    mask &= (
        adata.obs["n_total_de_genes"] >= args.min_degs
    )  # min DEGs of perturbation
    mask &= (
        adata.obs["n_cells_target"] >= args.min_cells
    )  # min cells with this perturbation
    mask &= adata.obs["ontarget_significant"].astype(
        bool
    )  ## significant on-target effect

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

    adj_p = layer_to_df(on_target, "adj_p_value", "padj")
    base_mean = layer_to_df(on_target, "baseMean", "base_DESeq2")
    lfc_se = layer_to_df(on_target, "lfcSE", "lfcSE")
    log_fc = layer_to_df(on_target, "log_fc", "logFC")
    zscore = layer_to_df(on_target, "zscore", "zscore")

    de_df = pl.concat(
        [
            adj_p,
            base_mean.select("base_DESeq2"),
            lfc_se.select("lfcSE"),
            log_fc.select("logFC"),
            zscore.select("zscore"),
        ],
        how="horizontal",
    )

    n_obs = obs_df.height
    n_var = var_df.height
    obs_idx = pl.Series("obs_idx", [i for i in range(n_obs) for _ in range(n_var)])
    var_idx = pl.Series("var_idx", list(range(n_var)) * n_obs)

    de_full = de_df.with_columns(
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
    ).drop("obs_id", "var_id")

    # Strong downstream effect filter
    # padj < threshold  AND  |zscore| >= z_thresh
    # Exclude on-target row (cis_ensg == trans_ensg) from targets
    print("Applying downstream effect filters …")
    de_filtered = (
        de_full.filter(pl.col("padj") < args.padj_threshold)
        .filter(pl.col("zscore").abs() >= args.z_thresh)
        .filter(pl.col("cis_ensg") != pl.col("trans_ensg"))  # exclude on-target
    )
    print(f"  rows after downstream filter: {len(de_filtered):,}")

    de_filtered.write_csv(args.out)


if __name__ == "__main__":
    main()
