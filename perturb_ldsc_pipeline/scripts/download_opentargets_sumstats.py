#!/usr/bin/env python
"""Download GWAS sumstats listed in GWAS_study_info.csv.

Reads --study-info, downloads each `summarystatsLocation` to
{out_dir}/{study}.tsv.gz, and skips files that already exist.

Handles two URL schemes:
  * gs://bucket/path  → rewritten to https://storage.googleapis.com/bucket/path
  * https://...       → passed straight through

Uses curl under the hood (works inside the apptainer container; no
extra Python deps needed).
"""

import argparse
import shutil
import sys
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import polars as pl

MIN_OK_BYTES = 100_000  # anything smaller is almost certainly an error page


def gs_to_https(url):
    if url.startswith("gs://"):
        url = "https://storage.googleapis.com/" + url[len("gs://") :]
    # FinnGen R12 public files live at summary_stats/release/<file> but
    # OpenTargets records just summary_stats/<file>. Inject the missing
    # path segment for those URLs.
    if (
        "/finngen-public-data-" in url
        and "/summary_stats/" in url
        and "/summary_stats/release/" not in url
    ):
        url = url.replace("/summary_stats/", "/summary_stats/release/", 1)
    return url


def download(url, out_path):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_path.with_suffix(out_path.suffix + ".part")

    req = urllib.request.Request(
        url,
        headers={"User-Agent": "perturb_ldsc_pipeline/1.0"},
    )

    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                with open(tmp, "wb") as f:
                    shutil.copyfileobj(resp, f, length=1024 * 1024)
            break
        except (urllib.error.URLError, TimeoutError):
            if attempt == 2:
                raise
            time.sleep(5)

    if tmp.stat().st_size < MIN_OK_BYTES:
        tmp.unlink(missing_ok=True)
        raise RuntimeError(
            f"downloaded file < {MIN_OK_BYTES} bytes — likely an error page"
        )
    tmp.rename(out_path)


def process(study, url, out_dir):
    out_path = out_dir / f"{study}.tsv.gz"
    if out_path.exists() and out_path.stat().st_size >= MIN_OK_BYTES:
        return study, "skip", out_path.stat().st_size, ""
    try:
        download(gs_to_https(url), out_path)
        return study, "ok", out_path.stat().st_size, ""
    except Exception as e:
        return study, "fail", 0, str(e)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--study-info", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Parallel downloads (default 4; do not crank too high — these come from a shared FTP)",
    )
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pl.read_csv(args.study_info)
    tasks = list(zip(df["study"].to_list(), df["summarystatsLocation"].to_list()))
    print(f"Downloading {len(tasks)} files into {out_dir} (threads={args.threads})")

    ok = skip = fail = 0
    failures = []
    with ThreadPoolExecutor(max_workers=args.threads) as ex:
        futs = [ex.submit(process, s, u, out_dir) for s, u in tasks]
        for fut in as_completed(futs):
            study, status, size, err = fut.result()
            if status == "ok":
                ok += 1
                print(f"  [OK]    {study:30s}  {size/1e6:7.1f} MB", flush=True)
            elif status == "skip":
                skip += 1
                print(
                    f"  [SKIP]  {study:30s}  {size/1e6:7.1f} MB  (already present)",
                    flush=True,
                )
            else:
                fail += 1
                failures.append((study, err))
                print(f"  [FAIL]  {study:30s}  {err}", file=sys.stderr, flush=True)

    print(f"\nDone: {ok} downloaded, {skip} skipped, {fail} failed.")
    if failures:
        print("\nFailures:", file=sys.stderr)
        for s, e in failures:
            print(f"  {s}: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
