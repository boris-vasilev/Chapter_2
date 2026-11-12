"""Hub definition shared by the pipeline (extract_de_results.py) and the
chapter figures (vis/). numpy only, so it runs in every environment."""

import numpy as np


def hub_mask(mask, n_degs, cond, genes, q, label):
    """Perturbations above the q-quantile of n_total_de_genes within their
    culture_condition, among those passing `mask`."""
    hub = np.zeros_like(mask)
    print(f"  {label} (q={q}, per culture_condition):")
    for c in np.unique(cond[mask]):
        in_c = mask & (cond == c)
        cutoff = np.quantile(n_degs[in_c], q)
        hub_c = in_c & (n_degs > cutoff)
        hub |= hub_c
        print(
            f"    {c}: cutoff {cutoff:.0f} DEGs, {hub_c.sum()} / {in_c.sum()} obs: {sorted(genes[hub_c])}"
        )
    return hub
