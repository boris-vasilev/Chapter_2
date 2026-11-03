#!/usr/bin/env python
"""Draft a GWAS → OpenTargets studyId mapping by pubmedId.

Assumes GWAS study names follow ``TRAIT_FirstAuthor_PUBMEDID_INDEX``
(e.g. ``UC_DeLange_28067908_1``). Extracts the pubmedId from each
study name, queries OpenTargets ``study/*.parquet``, and writes
candidate matches.

The output is a *draft* — manual curation is essential because:
  * a single pubmedId often has many OT studyIds (different ancestries,
    sub-traits, replication phases, etc.);
  * some of your studies may not have an exact pubmedId match (e.g.
    DOI-based identifiers like AST_Sakaue_doi...). Those land in the
    "unmatched" warning list — you'll need to look them up manually.

Output CSV columns:
  study              your GWAS study identifier
  studyId            candidate OpenTargets studyId
  traitFromSource    OT trait label (for human review when curating)
  pubmedId           the pubmedId used to find this match
  studyType          GWAS / molQTL / etc.
"""

import argparse
import re
import sys
from pathlib import Path

import polars as pl


PUBMED_RX = re.compile(r"_(\d{7,9})_\d+$")


def extract_pubmed(study):
    m = PUBMED_RX.search(study)
    return m.group(1) if m else None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--study-info",      required=True,
                        help="GWAS_study_info.csv with a 'study' column")
    parser.add_argument("--opentargets-dir", required=True,
                        help="Root OpenTargets dir containing study/*.parquet")
    parser.add_argument("--out",             required=True,
                        help="Output draft mapping CSV")
    args = parser.parse_args()

    studies = pl.read_csv(args.study_info)["study"].to_list()
    studies_with_pubmed = [(s, extract_pubmed(s)) for s in studies]

    unmatched = [s for s, p in studies_with_pubmed if p is None]
    if unmatched:
        print(
            f"WARNING: could not extract pubmedId from {len(unmatched)} "
            f"study name(s) — curate manually: {unmatched}",
            file=sys.stderr,
        )

    ot_studies = (
        pl.scan_parquet(str(Path(args.opentargets_dir) / "study" / "*.parquet"))
        .select(["studyId", "pubmedId", "traitFromSource", "studyType"])
        .collect()
    )

    rows = []
    for study, pmid in studies_with_pubmed:
        if pmid is None:
            continue
        matches = ot_studies.filter(pl.col("pubmedId") == pmid)
        if matches.height == 0:
            print(
                f"WARNING: no OpenTargets studies for pubmed {pmid} ({study})",
                file=sys.stderr,
            )
            continue
        for row in matches.to_dicts():
            rows.append({
                "study":           study,
                "studyId":         row["studyId"],
                "traitFromSource": row.get("traitFromSource"),
                "pubmedId":        pmid,
                "studyType":       row.get("studyType"),
            })

    out = pl.DataFrame(rows) if rows else pl.DataFrame(
        schema={"study": pl.Utf8, "studyId": pl.Utf8,
                "traitFromSource": pl.Utf8, "pubmedId": pl.Utf8,
                "studyType": pl.Utf8}
    )
    out.write_csv(args.out)

    print(f"Wrote {out.height} candidate mappings to {args.out}", file=sys.stderr)
    print(
        "Curate manually — drop irrelevant matches, keep one or more "
        "studyIds per study. The final mapping needs columns 'study' and "
        "'studyId' (extras are ignored).",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
