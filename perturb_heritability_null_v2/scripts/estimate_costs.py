#!/usr/bin/env python
"""Cost estimates for v2 permutation runs (no jobs are run).

One slurm job = one perm x condition chain (perm_chain rule), plus one `real`
chain per study x condition. Per-job inputs are measured (see NOTES):
  --cpu-h-per-chain     CPU-hours used per chain (sacct TotalCPU)
  --wall-min-per-chain  wall minutes per chain job
  --alloc-cpus          cpus per chain job (cpus_per_task)
  --overhead-min        per-job scheduling + status-poll latency (minutes)
  --concurrency         concurrent slurm jobs (profile `jobs`)
  --intermediate-mb     annotations + LD scores per chain while it runs
  --kept-mb             gene sets + .results/.log/.cov kept per chain
  --pilot-cpu-h-per-chain  optional second CPU-h basis (pilot sacct)
  --actual-*            optional measured totals of a finished run, written as
                        an extra `stage2_actual` row (sacct; see NOTES)

  uv run python scripts/estimate_costs.py --out ../data/permutations/plan_costs.csv ...
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import polars as pl

ALL_TRAITS = 18
IMD_SUBSET = ["IBD", "RA (EUR)", "ASTCO", "T1D", "AITD", "SLE"]
CONTROLS = ["IS", "MYOI", "OA", "T2D", "ADFH", "SCZ"]

PLANS = {
    "stage2_pilot_rerun": ("AITD x {Rest,Stim8hr,Stim48hr} + ASTCO x {Rest,Stim48hr}", 5),
    "A_all_traits_3_conditions": (f"{ALL_TRAITS} traits x 3 conditions", ALL_TRAITS * 3),
    "B_imd_subset_plus_controls_stim8": (
        f"{'; '.join(IMD_SUBSET)} + controls {'; '.join(CONTROLS)} at Stim8hr",
        len(IMD_SUBSET) + len(CONTROLS)),
}


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", required=True)
    p.add_argument("--cpu-h-per-chain", type=float, required=True)
    p.add_argument("--wall-min-per-chain", type=float, required=True)
    p.add_argument("--alloc-cpus", type=int, default=8)
    p.add_argument("--overhead-min", type=float, default=4.0)
    p.add_argument("--concurrency", type=int, default=20)
    p.add_argument("--intermediate-mb", type=float, required=True)
    p.add_argument("--kept-mb", type=float, required=True)
    p.add_argument("--n-perm", type=int, nargs="+", default=[200, 500])
    p.add_argument("--pilot-cpu-h-per-chain", type=float, default=None)
    p.add_argument("--actual-jobs", type=int, default=None)
    p.add_argument("--actual-cpu-h", type=float, default=None)
    p.add_argument("--actual-alloc-cpu-h", type=float, default=None)
    p.add_argument("--actual-wall-h", type=float, default=None)
    args = p.parse_args()

    rows = []
    if args.actual_jobs is not None:
        rows.append({
            "plan": "stage2_actual", "description": PLANS["stage2_pilot_rerun"][0]
            + " (measured, sacct)", "n_runs": PLANS["stage2_pilot_rerun"][1],
            "n_perm": 200, "slurm_jobs": args.actual_jobs,
            "cpu_hours_used": args.actual_cpu_h,
            "cpu_hours_allocated": args.actual_alloc_cpu_h,
            "concurrency": args.concurrency,
            "expected_wall_hours": args.actual_wall_h})
    for plan, (desc, n_runs) in PLANS.items():
        for n in args.n_perm:
            if plan == "stage2_pilot_rerun" and n != 200:
                continue
            jobs = n_runs * (n + 1)
            waves = math.ceil(jobs / args.concurrency)
            wall_h = waves * (args.wall_min_per_chain + args.overhead_min) / 60
            rows.append({
                "plan": plan, "description": desc, "n_runs": n_runs, "n_perm": n,
                "slurm_jobs": jobs,
                "cpu_hours_used": round(jobs * args.cpu_h_per_chain, 1),
                "cpu_hours_allocated": round(
                    jobs * args.alloc_cpus * args.wall_min_per_chain / 60, 1),
                "concurrency": args.concurrency,
                "expected_wall_hours": round(wall_h, 1),
                "peak_intermediate_gb": round(
                    args.concurrency * args.intermediate_mb / 1024, 2),
                "kept_outputs_gb": round(jobs * args.kept_mb / 1024, 2),
                "peak_disk_gb": round(
                    (args.concurrency * args.intermediate_mb + jobs * args.kept_mb) / 1024, 2),
            })
            if args.pilot_cpu_h_per_chain is not None:
                rows[-1]["cpu_hours_used_pilot_basis"] = round(
                    jobs * args.pilot_cpu_h_per_chain, 1)
    df = pl.DataFrame(rows, infer_schema_length=None)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    df.write_csv(args.out)
    with pl.Config(tbl_cols=-1, tbl_width_chars=220):
        print(df.drop("description"))


if __name__ == "__main__":
    main()
