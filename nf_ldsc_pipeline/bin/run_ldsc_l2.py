#!/usr/bin/env python
"""Run LDSC ``--l2`` on chr 1..22 in parallel via subprocess.Popen.

Mirrors the original 6_calculate_ld_scores.py loop so a single
Snakemake job covers a whole (study, condition) instead of one job
per chromosome.
"""

import argparse
import subprocess
import time


def run_capped(cmds, max_parallel):
    """Run (label, cmd) pairs with at most `max_parallel` concurrent Popens.

    Polls for *any* finished child before launching the next, so an uneven
    chr1 (much slower than chr22) never blocks the rest from starting.
    Returns the list of (label, returncode) for non-zero exits.
    """
    max_parallel = max(1, max_parallel)
    running = []  # list of (label, Popen)
    failures = []
    for label, cmd in cmds:
        while len(running) >= max_parallel:
            for i, (lab, p) in enumerate(running):
                rc = p.poll()
                if rc is not None:
                    if rc != 0:
                        failures.append((lab, rc))
                    running.pop(i)
                    break
            else:
                time.sleep(0.5)
        print("Spawning:", " ".join(cmd), flush=True)
        running.append((label, subprocess.Popen(cmd)))
    for label, p in running:
        rc = p.wait()
        if rc != 0:
            failures.append((label, rc))
    return failures


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ldsc",       required=True, help="Path to ldsc.py")
    parser.add_argument("--bim-prefix", required=True,
                        help="bim file prefix; chromosome appended as .{chr}")
    parser.add_argument("--snp-list",   required=True)
    parser.add_argument("--prefix",     required=True,
                        help="Shared prefix for {prefix}.{chr}.annot inputs and outputs")
    parser.add_argument("--max-parallel", type=int, default=22,
                        help="Max concurrent chromosome processes (default 22 = "
                             "all at once, the original behaviour).")
    args = parser.parse_args()

    cmds = [
        (chrom, [
            "python", args.ldsc,
            "--l2",
            "--bfile",       f"{args.bim_prefix}.{chrom}",
            "--ld-wind-cm",  "1",
            "--annot",       f"{args.prefix}.{chrom}.annot",
            "--out",         f"{args.prefix}.{chrom}",
            "--print-snps",  args.snp_list,
            "--thin-annot",
        ])
        for chrom in range(1, 23)
    ]
    failures = run_capped(cmds, args.max_parallel)
    if failures:
        raise SystemExit(f"ldsc --l2 failed on chrom(s): {failures}")


if __name__ == "__main__":
    main()
