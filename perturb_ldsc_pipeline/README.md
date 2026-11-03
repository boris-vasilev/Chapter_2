# Perturb-seq + GWAS S-LDSC pipeline

Snakemake pipeline that uses CD4+ T-cell Perturb-seq differential-expression
results to build SNP annotations and tests them against GWAS sumstats with
stratified LD score regression (S-LDSC). Two interchangeable gene-set modules
are supported:

| Module | Anchor genes | Source |
|---|---|---|
| `coloc` | Coloc-anchored: H3 (independent) and H4 (colocalised) genes per disease | `data/coloc_dat.tsv` |
| `l2g`   | OpenTargets Locus-to-Gene scored genes per disease | `data/opentargets/*.parquet` |

For each (study × condition × anchor set) the pipeline produces three annotation
columns:

- `<anchor>_genes` – the anchor genes themselves (e.g. `H4_genes`, `l2g_genes`)
- `regulators_<anchor>` – Perturb-seq cis genes whose perturbation changes
  expression of an anchor gene
- `targets_<anchor>` – DEGs (trans genes) when an anchor gene is itself perturbed

These three are combined into one .annot file per chromosome and run through
LDSC `--l2` and `--h2`.

---

## Repository layout

```
perturb_ldsc_pipeline/
├── Snakefile                       Pipeline entry point
├── config.yaml                     All tunable knobs
├── README.md                       This file
├── rules/
│   ├── common.smk                  Path/constant helpers, STUDIES list
│   ├── extract_de.smk              h5ad -> DE_results_sig.csv
│   ├── gene_sets_coloc.smk         coloc anchor sets
│   ├── gene_sets_l2g.smk           L2G anchor sets
│   ├── sumstats.smk                LDSC munge
│   ├── annots.smk                  make_annot + combine_annots
│   ├── ld_scores.smk               LDSC --l2
│   ├── filter_annots.smk           optional post-LD-score column filter
│   └── heritability.smk            LDSC --h2
├── scripts/                        Pipeline-internal scripts called by Snakemake rules
│   ├── extract_de_results.py
│   ├── make_gene_sets_coloc.py
│   ├── make_gene_sets_l2g.py
│   ├── run_make_annot.py           wraps LDSC make_annot, 22-chr Popen fan-out
│   ├── combine_annots.py           horizontal concat + low-SNP filter
│   ├── run_ldsc_l2.py              wraps ldsc --l2, 22-chr Popen fan-out
│   ├── run_munge_sumstats.py       FinnGen / GWC harmonised / GWC raw autodetect
│   ├── filter_annots.py            post-LD-score column filter (rarely used)
│   ├── select_imds_from_opentargets.py    L2G step 1
│   ├── compile_studies_from_opentargets.py L2G step 2
│   ├── download_opentargets_sumstats.py    L2G step 3
│   └── verify_sumstats_build.py            L2G optional QC
├── vis/                            Post-hoc analysis & plotting (run manually)
│   ├── parse_ldsc_results.py       results -> tidy CSV
│   └── plot_ldsc_results.py        heatmap + forest plot
└── slurm/
    ├── run.sbatch                  Single-job allocation runner (default driver)
    ├── submit_all.sh               Fan-out: one sbatch per study
    └── config.yaml                 Snakemake slurm-executor profile (alternative)
```

---

## How it works

### DAG (per study × condition)

```
                         extract_de_results
                                 │
                                 v
                          make_gene_sets    <-- coloc or l2g
                                 │
        ┌────────────────────────┴────────────────────────┐
        v                                                 v
  munge_sumstats                                    make_annot
        │                                                 │   (3 annot types ×
        │                                                 │    22 chromosomes
        │                                                 │    per study × cond,
        │                                                 │    fan-out inside)
        │                                                 v
        │                                          combine_annots
        │                                                 │   (drops sparse annots
        │                                                 │    with < min_snps total)
        │                                                 v
        │                                       calculate_ld_scores
        │                                                 │
        └──────────────────────┬──────────────────────────┘
                               v
                       calculate_heritability
                               │
                               v
            data/ldsc_results/<module>/{study}/{cond}_combined.results
```

### Path resolution

All paths in `config.yaml` are relative to `workspace`. Inside `common.smk`
we resolve via `Path.resolve()`, which follows the CSD3 symlinks (`/home/<user>/rds/...`
becomes `/rds/project/...`). The slurm runner binds both forms so the same
path works inside and outside the apptainer container.

### Module subtree

Derived outputs (gene sets, annotations, LD scores, results) live under a
`{module}/` subdirectory so coloc and L2G never collide. Flip
`gene_set_module: coloc | l2g` in `config.yaml` to switch.

Sumstats and the raw DE results stay outside the module subtree — they are
shared.

---

## One-time setup

### 1. Container

The pipeline runs nearly every step inside an apptainer container built from
`../ldsc.def`. The container ships LDSC + Python 3.9 + polars + scanpy +
matplotlib.

```bash
cd ..                            # at perturb_heritability_model/
mkdir -p containers
apptainer build --fakeroot containers/ldsc39_env.sif ldsc.def
```

If the build fails because of conda environment activation, ensure `ldsc.def`'s
`%environment` block points PATH at `/opt/conda/envs/ldsc/bin` (so `python` and
friends resolve to the right binary under `apptainer exec`, which doesn't run
`%runscript`).

### 2. Reference data

These are referenced by `config.yaml: bim / snp_list / baseline_ld / weights /
allele_frq / gene_coord`. Standard LDSC reference set:

```
data/
├── GRCh37/
│   ├── 1000G_Phase3_plink/1000G.EUR.QC.{1..22}.{bed,bim,fam}
│   ├── 1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.{1..22}.l2.ldscore.gz
│   └── baselineLD_v2.2/baselineLD.{1..22}.{annot.gz,l2.ldscore.gz,l2.M,l2.M_5_50}
├── GRCh37_EAS/                  # same structure for EAS
├── GRCh38/
│   ├── 1000G_Phase3_plink/1000G.EUR.hg38.{1..22}.{bed,bim,fam}
│   ├── 1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.{1..22}.l2.ldscore.gz
│   └── baselineLD_v2.2/baselineLD.{1..22}.*
├── hm3_no_MHC.list.txt
├── hm3_EAS_no_MHC.list.txt
├── gene_coord.txt               # GRCh37, from Ensembl 75 GTF
└── gene_coord_hg38.txt          # GRCh38, from Ensembl 115 GTF
```

Gene coord files are generated once by `../4_create_gene_coord_file.py` (not
part of the Snakemake DAG).

### 3. Perturb-seq DE input

Place at `data/GWCD4i.DE_stats.h5ad`. The `extract_de_results` rule keeps
gRNAs with `ontarget_effect_category == "on-target KD"` and rows with
`padj < padj_threshold` (default 0.05).

---

## Workflow A — coloc-anchored

### Inputs

| File | Description |
|---|---|
| `data/GWAS_study_info.csv` | Curated study list (`study,N,N_cases,N_controls,ancestry,build`) |
| `data/GWAS/{study}.tsv.gz` | Raw sumstats for each study |
| `data/coloc_dat.tsv` | Coloc results with `cell_type, imd, Gene_id, H3, H4` |
| `data/GWCD4i.DE_stats.h5ad` | Perturb-seq DE |

The `imd` column in `coloc_dat.tsv` must match `study` in
`GWAS_study_info.csv` exactly.

### Run

```bash
cd perturb_ldsc_pipeline

# Make sure config.yaml has:  gene_set_module: coloc

# One study end-to-end (replace UC_DeLange... with one of yours)
sbatch slurm/run.sbatch \
    /rds/project/rds-csoP2nj6Y6Y/biv22/perturb_heritability_model/data/ldsc_results/coloc/UC_DeLange_28067908_1/{Rest,Stim8hr,Stim48hr}_combined.results

# All studies (one sbatch per study)
MODULE=coloc bash slurm/submit_all.sh --dry-run   # confirm counts
MODULE=coloc bash slurm/submit_all.sh             # submit
```

---

## Workflow B — L2G-anchored (OpenTargets)

This module bootstraps its own GWAS study list from OpenTargets rather than
using a curated one. Four scripts run sequentially before the main pipeline.

### Step 1 — pick the IMDs

Drop the OpenTargets parquet export under `data/opentargets/` so the layout is:

```
data/opentargets/
├── disease/*.parquet
├── credible_set/*.parquet
├── l2g/*.parquet
├── study/*.parquet
└── target/*.parquet
```

Then:

```bash
cd ..      # at perturb_heritability_model/

# Inspect what therapeutic areas exist first
apptainer exec --fakeroot \
    -B /rds/project/rds-csoP2nj6Y6Y:/rds/project/rds-csoP2nj6Y6Y \
    -B /home/biv22/rds/rds-mrc-bsu-csoP2nj6Y6Y:/home/biv22/rds/rds-mrc-bsu-csoP2nj6Y6Y \
    containers/ldsc39_env.sif \
    python perturb_ldsc_pipeline/scripts/select_imds_from_opentargets.py \
        --opentargets-dir data/opentargets --inspect

# Filter to immune-system disease descendants (EFO_0000540)
apptainer exec ... \
    python perturb_ldsc_pipeline/scripts/select_imds_from_opentargets.py \
        --opentargets-dir data/opentargets \
        --therapeutic-area EFO_0000540 \
        --out data/opentargets/imd_list.csv

# Manually curate imd_list.csv if needed (keep only diseases of interest)
```

### Step 2 — compile GWAS study list

```bash
apptainer exec ... \
    python perturb_ldsc_pipeline/scripts/compile_studies_from_opentargets.py \
        --opentargets-dir data/opentargets \
        --imd-list data/opentargets/imd_list.csv \
        --out data/opentargets/GWAS_study_info.csv \
        --ancestry-filter EUR EAS
```

Filters: `studyType=gwas`, `hasSumstats=true`, at least one diseaseId in the
IMD list, at least one credible set, and (with `--ancestry-filter`) single
primary population mapped to EUR or EAS. gnomAD codes are mapped via:

| gnomAD | LDSC bucket |
|---|---|
| nfe | EUR |
| fin | EUR |
| eas | EAS |
| afr | AFR (not supported by the rest of the pipeline) |
| amr | AMR (not supported) |
| sas | SAS (not supported) |

Inspect the output: `head data/opentargets/GWAS_study_info.csv`.

*Note: GWAS studies from the GWAS catalog were marked as `hasSumstats=true` but\
had no `summarystatsLocation` specified. The sumstats destination was manually\
added for those from the GWAS catalog, taking the harmonised GRCh38 summary\
statistics.*

### Step 3 — download sumstats

```bash
salloc -A MRC-BSU2-SL2-CPU -p icelake --cpus-per-task=4 --time=4:00:00

apptainer exec ... \
    python perturb_ldsc_pipeline/scripts/download_opentargets_sumstats.py \
        --study-info data/opentargets/GWAS_study_info.csv \
        --out-dir data/opentargets/sumstats_raw \
        --threads 4
```

Handles both `gs://finngen-public-data-r12/...` (rewritten to
`https://storage.googleapis.com/.../summary_stats/release/<file>.gz`) and
`https://ftp.ebi.ac.uk/...` URLs. Resumable: rerun to fill in failures.

**Use no more than 4 threads** — EBI's FTP throttles aggressively and refuses
parallel connections beyond that.

### Step 4 — run the pipeline

```bash
cd perturb_ldsc_pipeline

# config.yaml:  gene_set_module: l2g

# Single study smoke test
sbatch slurm/run.sbatch \
    /rds/project/rds-csoP2nj6Y6Y/biv22/perturb_heritability_model/data/ldsc_results/l2g/GCST90476232/{Rest,Stim8hr,Stim48hr}_combined.results

# Fan out
MODULE=l2g bash slurm/submit_all.sh --dry-run
MODULE=l2g bash slurm/submit_all.sh
```

---

## Configuration (`config.yaml`)

Key knobs:

```yaml
workspace: "..."                        # project root (absolute)
container: "containers/ldsc39_env.sif"
gene_set_module: "coloc" | "l2g"

# DE filter
padj_threshold: 0.05

# Coloc-specific
coloc:
  study_info:    "data/GWAS_study_info.csv"
  gwas_dir:      "data/GWAS"
  sumstats_dir:  "data/sumstats"
  coloc_dat:     "data/coloc_dat.tsv"
  cell_type:     "CD4"
  h3_threshold:  0.5
  h4_threshold:  0.5
  h4_only:       false                  # set true to drop H3 annotations

# L2G-specific
l2g:
  study_info:      "data/opentargets/GWAS_study_info.csv"
  gwas_dir:        "data/opentargets/sumstats_raw"
  sumstats_dir:    "data/opentargets/sumstats"
  data_dir:        "data/opentargets"
  score_threshold: 0.5

# Annotation
window_size: 10000                      # bp added either side of each gene

# Conditions to model
conditions: [Rest, Stim8hr, Stim48hr]

# Reference data paths (per ancestry × build)
bim:        { EUR_GRCh37: "...", EAS_GRCh37: "...", EUR_GRCh38: "..." }
snp_list:   { EUR: "...", EAS: "..." }
baseline_ld:{ EUR_GRCh37: "...", EUR_GRCh38: "...", EAS_GRCh37: "..." }
weights:    { ... }
allele_frq: { ... }
gene_coord: { GRCh37: "...", GRCh38: "..." }

# Sparse-annotation filtering (applied during combine_annots)
min_snps_per_annot: 1000

# Post-LD-score filtering (rarely needed; off by default)
filter_low_snp_annots: false
```

---

## SLURM runner

`slurm/run.sbatch` is the canonical driver:

```
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=22 --time=1:30:00
```

It runs Snakemake **on the host** (your `env_nf` conda env must be active for
sbatch jobs) and Snakemake invokes the apptainer container per rule via
`software-deployment-method: apptainer`. No Snakemake slurm executor / no
recurring `squeue` polling — each job is a single allocation that completes
when it's done.

Per-disease wall time on 22 CPUs: ~20-30 minutes.

### CSD3 etiquette

CSD3 admins have explicitly asked users not to hammer slurmctld with rapid
queries. The pipeline avoids this by:

- Single-allocation mode (no per-rule sbatch).
- `submit_all.sh` paces `sbatch` submissions at 1/sec.
- No `watch squeue` anywhere — check status manually with `squeue -u $USER`
  no more often than once every couple of minutes, or use the filesystem:

  ```bash
  ls data/ldsc_results/<module>/*/Rest_combined.results | wc -l
  ```

The legacy `slurm/config.yaml` Snakemake-slurm profile is kept for reference
but **not the recommended runner** — it polls slurm continuously and was the
source of past warnings from CSD3 admins.

---

## Output structure

```
data/
├── DE_results_sig.csv                                 # padj-filtered DE
├── gene_sets/<module>/{study}_{cond}_{annot}.txt      # one rsID per line
├── annotations/<module>/{study}/{cond}_{annot}.{chr}.annot.gz
├── ld_scores/<module>/{study}/{cond}_combined.{chr}.{annot,l2.ldscore.gz,l2.M_5_50}
├── ldsc_results/<module>/{study}/{cond}_combined.results
└── (l2g only) opentargets/sumstats/{study}.sumstats.gz       # munged
└── (coloc only) sumstats/{study}.sumstats.gz                 # munged
```

The `.results` file has these columns (LDSC standard):
`Category, Prop._SNPs, Prop._h2, Prop._h2_std_error, Enrichment,
Enrichment_std_error, Enrichment_p, Coefficient, Coefficient_std_error,
Coefficient_z-score`.

Our annotation rows have names like `regulators_l2gL2_1`,
`targets_l2gL2_1`, `l2g_genesL2_1`, plus all the baseline categories.

---

## Plotting

```bash
# 1. Parse all .results into a single tidy CSV
apptainer exec ... \
    python perturb_ldsc_pipeline/vis/parse_ldsc_results.py \
        --results-dir data/ldsc_results/l2g \
        --study-info  data/opentargets/GWAS_study_info.csv \
        --annots regulators_l2g targets_l2g l2g_genes \
        --out         data/ldsc_results/l2g_summary.csv

# 2. Heatmap + forest plot
apptainer exec ... \
    python perturb_ldsc_pipeline/vis/plot_ldsc_results.py \
        --summary data/ldsc_results/l2g_summary.csv \
        --out-dir data/ldsc_results \
        --annots regulators_l2g targets_l2g l2g_genes \
        --prefix l2g
```

For coloc:

```bash
apptainer exec ... \
    python perturb_ldsc_pipeline/vis/parse_ldsc_results.py \
        --results-dir data/ldsc_results/coloc \
        --study-info  data/GWAS_study_info.csv \
        --annots regulators_H4 targets_H4 H4_genes \
                 regulators_H3 targets_H3 H3_genes \
        --out         data/ldsc_results/coloc_summary.csv

apptainer exec ... \
    python perturb_ldsc_pipeline/vis/plot_ldsc_results.py \
        --summary data/ldsc_results/coloc_summary.csv \
        --out-dir data/ldsc_results \
        --annots regulators_H4 targets_H4 H4_genes \
                 regulators_H3 targets_H3 H3_genes \
        --label-col trait_short \
        --prefix coloc
```

`plot_ldsc_results.py` picks row labels in order of preference:
`traitFromSource` (L2G) → `trait_short` (study prefix, e.g. `UC`) →
`study` (the raw identifier).

Outputs:

- `<prefix>_zscore_heatmap.pdf` — trait-mean coefficient z-score, rows
  sorted by signal strength
- `<prefix>_enrichment_forest.pdf` — per-annot forest plot of enrichment
  (log scale) ± 95% CI, with the three conditions offset vertically

### Gene-set overlap between studies / traits

`vis/study_overlap.py` measures how similar the per-study gene sets are to
each other (Jaccard index + overlap coefficient) for each (condition,
annotation) panel. Useful sanity check — typically the anchor sets
(`H4_genes`, `l2g_genes`) overlap little across diseases (different genes
prioritised per trait) but the perturbation-network sets (`regulators_*`,
`targets_*`) overlap substantially because many trans-effects converge on
shared regulatory hubs.

Coloc (18 studies, labelled by raw study name):

```bash
apptainer exec ... \
    python perturb_ldsc_pipeline/vis/study_overlap.py \
        --gene-sets-dir data/gene_sets/coloc \
        --annots H4_genes regulators_H4 targets_H4 \
        --out-prefix    plots/gene_set_analysis/coloc
```

L2G (~99 studies → ~17 traits by unioning gene sets across studies of the
same `traitFromSource`):

```bash
apptainer exec ... \
    python perturb_ldsc_pipeline/vis/study_overlap.py \
        --gene-sets-dir data/gene_sets/l2g \
        --study-info    data/opentargets/GWAS_study_info.csv \
        --annots l2g_genes regulators_l2g targets_l2g \
        --label-col     traitFromSource \
        --out-prefix    plots/gene_set_analysis/l2g
```

Outputs (per call):
- `<prefix>_jaccard.png` — pairwise Jaccard index, conditions on rows
  × annotations on cols
- `<prefix>_overlap.png` — pairwise overlap coefficient
  (intersection / size of the smaller set) — better than Jaccard when the
  two sets differ in size

Diagonal is masked. With ≤ 25 labels the cell values are printed inside the
heatmap; with more, only the colours are shown (keeps the L2G plot
readable).

`--label-col` controls how studies are grouped. Pass `study` (default) to
keep every study as a separate row; pass `traitFromSource` or `trait_short`
to collapse multiple studies of the same trait into one row by unioning
their gene sets. Useful when the per-study matrix would be too dense to
read.

---

## Common gotchas

**All zeros in the LD scores or NaN in h²**: an annotation column was empty
or too sparse. `combine_annots` now drops annotations with fewer than
`min_snps_per_annot` (default 1000) summed across chromosomes; rerun after
wiping `data/ld_scores/<module>/<study>` and `data/ldsc_results/<module>/<study>`.

**`(1, X) vs (1, Y)` shape mismatch in `--h2`**: column count mismatch between
`l2.M_5_50` and the annot matrix. Almost always caused by mixing pre- and
post-filter combined annots. Wipe the affected `data/ld_scores/...` dir and rerun.

**Snakemake wants to rerun everything**: probably the container hash changed.
Add `--rerun-triggers mtime` to ignore provenance (already on by default in
`slurm/run.sbatch`).

**EAS GRCh38 study fails at `make_annot`**: only `EUR_GRCh38` is in the
default reference set. Either download EAS GRCh38 1000G + baselineLD or
filter those studies out of `GWAS_study_info.csv`.

**FinnGen 404s on download**: OpenTargets records `gs://.../summary_stats/<f>.gz`
but the public bucket layout has them under
`gs://.../summary_stats/release/<f>.gz`. The downloader handles this; verify
by `curl -sI` on the rewritten URL.

**Sumstats munge runs out of memory**: each FinnGen file is ~20M rows. The
22 GB allocation in `run.sbatch` is comfortable; if a single bigger file blows
it up, bump `--mem` or `cpus-per-task` to get a bigger memory share.

---

## Acknowledgements / dependencies

- Snakemake (host, via `env_nf` conda env)
- Apptainer (container runtime)
- LDSC (Bulik-Sullivan & Finucane et al.)
- Polars + scanpy + matplotlib (in the container)

Pipeline runs against:

- 1000G Phase 3 reference panels
- baselineLD v2.2
- HapMap3 SNP list (MHC excluded)
- OpenTargets release (24.06 or later) — `l2g`, `credible_set`, `study`,
  `disease`, `target` parquet tables
