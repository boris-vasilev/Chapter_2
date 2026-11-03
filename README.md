# Chapter 2

`perturb_ldsc_pipeline/` is a copy of the S-LDSC Snakemake pipeline. LDSC itself
runs in the cluster container. `vis/` holds the chapter 2 figure scripts.

## Setup

```bash
cd chapter2
uv sync                       # creates .venv from pyproject.toml / uv.lock
```

## Figures

Everything (figures, tables, data, index, contact sheet) is rebuilt from scratch with:

```bash
cd vis
uv run 00_build_all.py            # optional: --schematic path/to/schematic.png
```

Outputs land in `figures/{main,supp,data}/` with `INDEX.md` (one row per file) and
`contact_sheet.pdf`. All outputs use the `coloc_thr0.5` gene sets (|logFC| ≥ 0.5 edges),
the scripts' default. Individual scripts (`01_…` to `07_…`) can be run on their own and
write to the gitignored staging folder `build/coloc_thr0.5/`.
Disease display names and the IMD / non-IMD split live in `vis/diseases.csv`.
Notes (definitions, results, open questions): `notes/NOTES_overview.md`.
# Chapter_2
