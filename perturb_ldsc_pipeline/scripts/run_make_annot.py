#!/usr/bin/env python
"""Run LDSC ``make_annot.py`` on chr 1..22 in parallel via subprocess.Popen.

Replaces the previous sequential bash loop so a single Snakemake job
fans out over chromosomes instead of running them serially.
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
    parser.add_argument("--make-annot",     required=True, help="Path to make_annot.py")
    parser.add_argument("--gene-set-file",  required=True)
    parser.add_argument("--gene-coord-file", required=True)
    parser.add_argument("--windowsize",     required=True)
    parser.add_argument("--bim-prefix",     required=True,
                        help="bim prefix; chromosome appended as .{chr}.bim")
    parser.add_argument("--out-prefix",     required=True,
                        help="output prefix; chromosome appended as .{chr}.annot.gz")
    parser.add_argument("--max-parallel", type=int, default=22,
                        help="Max concurrent chromosome processes (default 22 = "
                             "all at once, the original behaviour).")
    args = parser.parse_args()

    cmds = [
        (chrom, [
            "python", args.make_annot,
            "--gene-set-file",  args.gene_set_file,
            "--gene-coord-file", args.gene_coord_file,
            "--windowsize",     args.windowsize,
            "--bimfile",        f"{args.bim_prefix}.{chrom}.bim",
            "--annot-file",     f"{args.out_prefix}.{chrom}.annot.gz",
        ])
        for chrom in range(1, 23)
    ]
    failures = run_capped(cmds, args.max_parallel)
    if failures:
        raise SystemExit(f"make_annot failed on chrom(s): {failures}")


if __name__ == "__main__":
    main()
