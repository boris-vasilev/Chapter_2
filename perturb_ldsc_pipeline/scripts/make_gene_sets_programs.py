#!/usr/bin/env python
"""Gene-set module: Perturb-seq programmes (Zhu clusters).

Unlike the coloc / l2g modules, programme gene sets are **study-independent**
— a programme is the same set of perturbed genes regardless of which GWAS
trait it is later tested against. They are, however, condition-specific
(cluster membership differs between Rest / Stim8hr / Stim48hr).

For every (condition, cluster) with at least ``--min-program-size`` member
genes this writes::

    {out_dir}/{condition}_program{cluster}.txt

one Ensembl gene ID per line, no header (the ``perturb_ensg`` of the
cluster's members). It also writes ``{out_dir}/programs_index.csv`` mapping
each (condition, cluster) to its programme number, manual annotation and
size, for labelling downstream S-LDSC results.
"""

import argparse
from pathlib import Path

import polars as pl


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--programs", required=True,
                   help="data/programs_zhu/program_assignments.csv")
    p.add_argument("--conditions", nargs="+", required=True)
    p.add_argument("--min-program-size", type=int, default=10)
    p.add_argument("--out-dir", required=True)
    args = p.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pl.read_csv(args.programs)
    needed = {"perturb_ensg", "condition", "cluster", "manual_annotation"}
    missing = needed - set(df.columns)
    if missing:
        raise SystemExit(f"program_assignments missing columns: {missing}")
    has_program = "program" in df.columns

    df = (
        df.filter(pl.col("condition").is_in(args.conditions))
        .unique(subset=["condition", "cluster", "perturb_ensg"])
    )

    index_rows = []
    for (cond, cluster), grp in df.group_by(["condition", "cluster"]):
        genes = grp["perturb_ensg"].drop_nulls().unique().to_list()
        if len(genes) < args.min_program_size:
            continue
        cluster = int(cluster)
        meta = grp.row(0, named=True)
        prog = int(meta["program"]) if has_program and meta["program"] is not None else cluster
        fname = out_dir / f"{cond}_program{cluster}.txt"
        pl.DataFrame({"GENE": genes}).write_csv(
            fname, separator="\t", include_header=False
        )
        index_rows.append({
            "condition": cond,
            "cluster": cluster,
            "program": prog,
            "manual_annotation": meta["manual_annotation"] or f"cluster {cluster}",
            "n_genes": len(genes),
        })

    if not index_rows:
        raise SystemExit(
            f"No programmes reached --min-program-size={args.min_program_size}."
        )

    index = pl.DataFrame(index_rows).sort(["condition", "cluster"])
    index.write_csv(out_dir / "programs_index.csv")
    print(f"Wrote {len(index_rows)} programme gene sets + programs_index.csv "
          f"to {out_dir}")
    for cond in args.conditions:
        n = index.filter(pl.col("condition") == cond).height
        print(f"  {cond}: {n} programmes (>= {args.min_program_size} genes)")


if __name__ == "__main__":
    main()
