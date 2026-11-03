# Chapter 2 overview figures — notes

Source of this file: `chapter2/notes/NOTES_overview.md`. `00_build_all.py` copies it into
`figures/`. Edit the source, not the copy.

## Build

```bash
cd chapter2/vis
uv run 00_build_all.py                 # schematic: assets/ch2_schematic.png
```

The build empties `figures/` and `build/coloc_thr0.5/`, runs scripts 01–07 (each writes to
`build/`), copies the chapter outputs into `figures/` under their figure numbers, and writes
`INDEX.md` and `contact_sheet.pdf`. It fails if any file in `figures/` is missing from the
index. A clean run takes about 45 s.

```
figures/
  main/   Table 2.1, Fig. 2.1, Fig. 2.2, Fig. 2.6A
  supp/   Fig. S2A–D, S3, S4
  data/   every CSV behind a figure
  INDEX.md, contact_sheet.pdf, NOTES_overview.md
```

See `INDEX.md` for the file-by-file list (figure number, section, message, command).
PNG versions of every figure sit in `build/coloc_thr0.5/` (not copied to `figures/`).

---

## Variant used: `coloc_thr0.5` (only)

Every figure and table uses one gene-set variant (the scripts' default; `--variant` is kept
in the code, but the build produces nothing for other variants).

- **Anchors:** `coloc_dat.tsv` rows with `cell_type == "CD4"` and `H4 > 0.8`, Ensembl version stripped.
- **Edges:** `DE_results_strong.csv` (on-target KD, `ontarget_significant`,
  `n_total_de_genes ≥ 75`, `n_cells_target ≥ 50`, `padj < 0.05`, `|zscore| ≥ 2`, no
  self-edges) **plus `|logFC| ≥ 0.5`** (`make_gene_sets_coloc.py --min-abs-logfc 0.5`).
  This keeps 450,003 of 1,080,973 edges.
- **Regulators** of a disease = knockdowns with ≥1 edge into one of its anchors;
  **targets** = genes with ≥1 edge out of a knocked-down anchor (no convergence filter).
- **Check (measured):** rebuilding all 18 diseases × 3 conditions from this definition
  reproduces `data/gene_sets/coloc_thr0.5` exactly (162 / 162 sets, 0 genes differ).
  H4 > 0.5 does *not* reproduce it (152 / 162 sets differ), so "thr0.5" is the |logFC|
  threshold.
- **Universes** (for the expected Jaccard and the regulator reference lines): knockdowns
  with ≥1 edge (894 / 970 / 934 at Rest / 8 h / 48 h) and genes with ≥1 incoming edge
  (9,544 / 9,775 / 9,665). No list of all measured or perturbed genes exists locally (the
  DE h5ad is on the cluster), so these edge-based universes are a fallback.

Step 0 findings from the earlier round still hold: 05's edges equal the gene-set edges, and
3,341 is not reproduced by any definition tried (nearest 3,336, a coincidence).

---

## Display names (please check)

Used on every figure and in Table 2.1 (`display_name` in `vis/diseases.csv`). No
abbreviations appear on figures or in tables.

| Display name | Class | Study id |
|---|---|---|
| Autoimmune thyroid disease | IMD | AITD_Saevarsdottir_32581359_1 |
| Adult-onset asthma | IMD | ASTAO_Ferreira_30929738_1 |
| Childhood-onset asthma | IMD | ASTCO_Ferreira_30929738_1 |
| Asthma (EAS) | IMD | AST_Sakaue_doi1011012020102320213652_1 |
| Atopic dermatitis | IMD | ATD_BuduAggrey_37794016_1 |
| Crohn's disease | IMD | CD_DeLange_28067908_1 |
| Inflammatory bowel disease | IMD | IBD_DeLange_28067908_1 |
| Ulcerative colitis | IMD | UC_DeLange_28067908_1 |
| Rheumatoid arthritis (EAS) | IMD | RA_Ishigaki_36333501_1 |
| Rheumatoid arthritis (EUR) | IMD | RA_Ishigaki_36333501_2 |
| Systemic lupus erythematosus | IMD | SLE_Bentham_26502338_1 |
| Type 1 diabetes | IMD | T1D_Robertson_34127860_1 |
| Alzheimer's disease (proxy) | non-IMD | ADFH_Jansen_30617256_1 |
| Ischaemic stroke | non-IMD | IS_Sakaue_34594039_1 |
| Myocardial infarction | non-IMD | MYOI_Sakaue_34594039_1 |
| Osteoarthritis | non-IMD | OA_McDonald_36411363_3 |
| Schizophrenia | non-IMD | SCZ_Pardiñas_29483656_1 |
| Type 2 diabetes | non-IMD | T2D_Xue_30054458_1 |

Longest name: "Systemic lupus erythematosus" (28 characters). All disease labels are drawn
at 8–9 pt (≥ 7 pt everywhere). Figure widths and margins were increased so that no label
clips or overlaps.

## Diseases with empty target sets (measured)

| Disease | Conditions with zero targets |
|---|---|
| Asthma (EAS) | Rest, 8 h, 48 h |
| Ulcerative colitis | Rest, 8 h, 48 h |
| Osteoarthritis | Rest, 8 h, 48 h |
| Ischaemic stroke | Rest, 8 h, 48 h |
| Alzheimer's disease (proxy) | Rest, 8 h |
| Type 2 diabetes | Rest |

No disease has an empty regulator set. Fig. 2.2D shows empty target sets in a separate "0"
strip left of an axis break; Fig. 2.2E omits them; Fig. S3 excludes these six diseases from
the target panels.

---

## Measured results (coloc_thr0.5; Rest / 8 h / 48 h)

From `figures/data/`:

- **Degrees (overview_numbers.csv).** Out-degree: median 67 / 78 / 59, 99th percentile
  1,067 / 1,857 / 1,046, max 1,728 / 4,206 / 2,056. Top hubs: LCP2, ZAP70, CD3E.
  In-degree: median 7 / 12 / 7, 99th percentile 69 / 101 / 86, max 185 / 162 / 177.
  Top responders: TXNIP, SESN1, FDXR.
- **Set sizes (set_sizes.csv).** Median regulator set 142 / 217 / 182 genes; median target
  set 79 / 202 / 115 (non-empty and empty sets included). The median regulator set is
  16% / 22% / 19% of all knockdowns with an effect; the largest is 42% / 51% / 42%.
- **Genome coverage, median % of reference SNPs (genome_coverage.csv).** Anchors 0.05;
  regulators 0.38 / 0.59 / 0.40; targets 0.55 / 1.14 / 0.45 (Rest / 8 h / 48 h). Read from
  heritability_long.csv, so the single-set runs (osteoarthritis regulators, Alzheimer's
  proxy 48 h targets) are included. There are no S-LDSC results for schizophrenia (all
  conditions) or the Alzheimer's proxy (Rest, 8 h); sets under 1,000 SNPs were not fitted.
- **Anchor coverage (anchor_coverage.csv), out of 344 anchors per condition.**
  Both edge types: 15 / 22 / 24. Incoming only: 225 / 215 / 217. Outgoing only: 1 / 1 / 0.
  No edge: 103 / 106 / 103.
- **Largest-anchor share (largest_anchor_share.csv).** Targets: 25 of the 39 non-empty sets
  come entirely from one anchor (minimum share 0.49). Regulators: median share 0.33
  (minimum 0.21). Recurrent dominant target anchors are RASGRP1, CD247, UBE2L3, TUFM and
  TRAF3.
- **Excess Jaccard, regulators (observed − expected).** At 8 h, every pair is positive
  (0 / 153 below zero; range 0.005–0.58). At Rest 4 / 153 pairs and at 48 h 1 / 153 pairs
  are marginally negative (≥ −0.003). Median for pairs of two IMDs vs all other pairs:
  0.20 vs 0.07 (Rest), 0.25 vs 0.11 (8 h), 0.27 vs 0.07 (48 h).
- **Excess Jaccard, targets.** Mostly near zero: median for IMD pairs 0.03, other pairs
  0.004. A few pairs reach about 1: the two rheumatoid arthritis sets and AITD (RASGRP1);
  IBD and T1D (TUFM); SLE and CD (UBE2L3).
- **Set size by condition (set_size_by_condition_tests.csv; tests are not on the figure).**
  Regulators: n = 18, median log2 FC vs Rest +0.63 (8 h), +0.21 (48 h). Targets: n = 12,
  +0.94 (8 h), +0.20 (48 h).

## Robustness: |logFC| filter vs no filter (for one sentence)

The same pipeline without the |logFC| filter (`data/gene_sets/coloc`, edges
`|logFC| ≥ 0`) gives the numbers below. No figure or table is produced for that variant.

| Quantity (Rest / 8 h / 48 h) | No |logFC| filter | |logFC| ≥ 0.5 (used) |
|---|---|---|
| Edges | 1,080,973 | 450,003 |
| Regulator universe | 895 / 970 / 939 | 894 / 970 / 934 |
| Target universe | 10,125 / 10,210 / 10,163 | 9,544 / 9,775 / 9,665 |
| Median regulator set, % of universe | 23 / 31 / 28 | 16 / 22 / 19 |
| Out-degree median | 145 / 177 / 167 | 67 / 78 / 59 |
| In-degree median | 23 / 38 / 26 | 7 / 12 / 7 |
| Median % SNPs, regulators | 0.67 / 0.92 / 0.69 | 0.41 / 0.69 / 0.41 |
| Median % SNPs, targets | 2.37 / 4.89 / 1.69 | 0.55 / 1.14 / 0.46 |
| Anchors with both edge types | 16 / 23 / 24 | 15 / 22 / 24 |
| Target sets from one anchor | 25 / 39 | 25 / 39 |
| Excess Jaccard at 8 h, regulators: IMD pairs vs other pairs | 0.21 vs 0.15 | 0.25 vs 0.11 |
| Median log2 FC vs Rest, regulators (8 h, 48 h) | +0.49, +0.35 | +0.63, +0.21 |
| Median log2 FC vs Rest, targets (8 h, 48 h) | +1.30, +0.26 | +0.94, +0.20 |
| Diseases with no targets in any condition | 4 (same four) | 4 |

Suggested sentence: *Requiring |logFC| ≥ 0.5 roughly halved node degrees and reduced target
coverage about four-fold, but left the anchors that dominate each disease's expanded sets
unchanged and widened the gap in regulator overlap between pairs of immune-mediated
diseases and other pairs.*

## Interpretation (mine, to check)

- **Target sets are mostly one knockdown's response.** At most 24 of 344 anchors are
  strong knockdowns, and most diseases' targets come from one of them. So high
  target-set overlap means two diseases share that anchor (RASGRP1, TUFM, UBE2L3, CD247).
  It does not show convergent biology.
- **Regulator overlap above random is partly a degree effect.** The expectation assumes
  uniform sampling, but hub knockdowns hit many anchors. The more informative contrast is
  the higher excess for pairs of two IMDs. It is not tested here, and diseases share
  anchors and duplicate traits (RA ×2, asthma ×3, IBD/CD/UC).
- **Stimulation enlarges sets most at 8 h**, matching the 8 h peak in out-degree (A).
- **The largest knockdowns are TCR-signalling genes** (LCP2, ZAP70, CD3E) once weak
  edges are removed.

---

## Citation keys for Table 2.1 (match to your .bib)

First-author-year keys (`cite_key` in `vis/diseases.csv`). Author and year were checked
against PubMed for each PMID.

| Key | Paper | Traits |
|---|---|---|
| `saevarsdottir2020` | Saevarsdottir et al., Nature 2020, doi:10.1038/s41586-020-2436-0 | Autoimmune thyroid disease |
| `ferreira2019` | Ferreira et al., Am J Hum Genet 2019, doi:10.1016/j.ajhg.2019.02.022 | Childhood- and adult-onset asthma |
| `buduaggrey2023` | Budu-Aggrey et al., Nat Commun 2023, doi:10.1038/s41467-023-41180-2 | Atopic dermatitis |
| `delange2017` | de Lange et al., Nat Genet 2017, doi:10.1038/ng.3760 | IBD, Crohn's disease, ulcerative colitis |
| `ishigaki2022` | Ishigaki et al., Nat Genet 2022, doi:10.1038/s41588-022-01213-w | Rheumatoid arthritis (EUR, EAS) |
| `bentham2015` | Bentham et al., Nat Genet 2015, doi:10.1038/ng.3434 | Systemic lupus erythematosus |
| `robertson2021` | Robertson et al., Nat Genet 2021, doi:10.1038/s41588-021-00880-5 | Type 1 diabetes |
| `jansen2019` | Jansen et al., Nat Genet 2019, doi:10.1038/s41588-018-0311-9 | Alzheimer's disease (proxy) |
| `sakaue2021` | Sakaue et al., Nat Genet 2021, doi:10.1038/s41588-021-00931-x | Ischaemic stroke, myocardial infarction, **Asthma (EAS)** |
| `mcdonald2022` | McDonald et al., Nat Genet 2022, doi:10.1038/s41588-022-01221-w | Osteoarthritis |
| `pardinas2018` | Pardiñas et al., Nat Genet 2018, doi:10.1038/s41588-018-0059-2 | Schizophrenia |
| `xue2018` | Xue et al., Nat Commun 2018, doi:10.1038/s41467-018-04951-w | Type 2 diabetes |

The Asthma (EAS) study id points to the medRxiv preprint doi:10.1101/2020.10.23.20213652.
I keyed it to `sakaue2021`, assuming that preprint is the Sakaue cross-population atlas
later published in Nature Genetics. Check this, and use a separate key if you cite the preprint.

**Shared samples** (`sample_group` in `diseases.csv`; Table 2.1 footnote, Fig. 2.6A
outlines): de Lange 2017 (IBD / CD / UC), Ferreira 2019 (childhood- / adult-onset
asthma), Ishigaki 2022 (RA EUR / EAS). This is the list from the brief. **Probably
also overlapping, but not flagged:** ischaemic stroke and myocardial infarction, both from
Sakaue 2021 with the same biobank cohorts, and possibly Asthma (EAS) with RA (EAS),
since both include Biobank Japan. Add them to `sample_group` if you agree.

## Print size (main figures)

Table 2.1, Fig. 2.1, Fig. 2.2 and Fig. 2.6A are drawn at final size: 6.3 in wide, text
7–9 pt (ticks 7, labels 7.5, titles 8, panel letters 9), lines ≥ 0.5 pt. Each PDF is
exactly 6.3 in wide: 2.1 is 6.3 × 8.95 in, 2.2 is 6.3 × 7.1 in, 2.6A is 6.3 × 3.6 in. Each
was checked by rasterising the PDF and viewing it at about 100–150% zoom.

To keep 7 pt at this width:
- **Heatmaps (2.1C, 2.6A):** rows show "Name  k" and columns show the number k (rotated);
  there are no rotated disease names on the x axis.
- **Fig. 2.1C** no longer prints the Jaccard value in each cell, since ".37" does not fit a
  0.125 in cell at 7 pt. Values are in `data/jaccard_anchors_observed.csv`.
- **Fig. 2.2:** one condition legend at the top for all panels. The n per condition moved
  from the CCDF legends into the caption. Hub labels carry their condition, e.g.
  "CD3E (8 h)".
- **Fig. 2.1B:** the class labels at the right edge are gone (the names are already
  coloured by class, and the legend in C explains the colours). The asterisk note is in
  the caption only.
- **Fig. 2.1B colour:** "Both" is now yellow `#eda100`, used nowhere else. The dataviz
  validator gives CVD ΔE 9.1 against its stack neighbour (green). Magenta failed
  (ΔE 6.1), and violet is already the IMD colour.
- **Fig. 2.6A:** the excess panel uses a sequential red scale from 0 to the maximum. No
  pair is below 0 at 8 h, so nothing is clipped; the colour bar would show an arrow if one
  were. The diagonal is masked (light grey) in both panels.
- **Fig. 2.2C** uses violins instead of boxplots. The density is estimated on log10(% SNPs) and drawn on a log10 axis labelled with plain numbers (`common.plain_log10_ticks`); each violin rests on 12–17 diseases.
- **Log-axis ticks:** each decade is labelled; 2 and 5 ticks are added only when a decade
  has ≥ 0.9 in of axis, so ticks never crowd.
- **Schematic:** at 6.3 in, its smallest text ("R1", "G1") is about 7 pt. Export the new
  PNG for the same width, at ≥ 300 dpi (≥ 1,900 px wide).

## What changed in earlier rounds

- `display_name` in `diseases.csv` is used everywhere; there are no abbreviations on
  figures or in Table 2.1.
- Log axes show plain numbers through `common.plain_log_axis`, and symlog is gone. Fig.
  2.2D uses a log axis with a "0" strip and an axis break.
- Fig. 2.2 layout: A, B, C on the top row; D and E on the second row, sharing the disease
  labels. Panels were re-lettered to keep reading order.
- `07_fig_2_1.py`, `04 --pair`, `00_build_all.py`, and the `main/ supp/ data/` layout with
  `INDEX.md` and `contact_sheet.pdf`.

## Could not do / decisions for you

1. **Supplementary figures are not at print size.** The brief covered Table 2.1, Fig. 2.1,
   2.2 and 2.6A, so S2–S4 keep their earlier large-format layout and fonts shrink below
   7 pt when scaled to 6.3 in. S2 (three 18 × 18 heatmaps side by side) would need a
   different layout at print size, e.g. stacked panels with numbered columns as in 2.6A.
2. **Fig. 2.1 is 8.95 in tall**, nearly a full page. The caption will probably have to go
   on the facing page, or the schematic must be shorter (its height follows the PNG's
   aspect ratio).
3. **Supplementary numbering starts at S2**, assuming S1 is used elsewhere. Shift it in
   `MANIFEST` (`00_build_all.py`) if not.
4. **S2 is four files (S2A–D)**: regulators / targets × observed / excess.
5. **Chapter section names in INDEX.md are descriptive placeholders.**
6. **Schematic (Fig. 2.1A)** comes from `chapter2/assets/ch2_schematic.png`. Replace that
   file with your new PNG and rerun the build; the current one has red spell-check
   underlines baked in.
7. **Linear axes** keep default tick labels (no thousands separators); none of the four
   main figures has a linear axis that reaches 1,000.
8. **PDF only in `figures/`.** PNGs stay in `build/`.

## TODO for Boris

1. **Table 2.1** (`figures/main/table_2_1_traits.tex`) is a complete `table` float with
   `\caption` and `\label{tab:traits}`: put `\input{figures/main/table_2_1_traits}` where
   the table belongs and refer to it as `Table~\ref{tab:traits}`. Preamble needs
   `\usepackage{booktabs,tabularx}`. The tabularx fits `\linewidth` (tested clean at 345,
   441 and 495 pt, numeric and author-year citations). One eQTL column only: **CD4⁺ eQTL
   studies** = number of distinct CD4 eQTL studies (blood_cell suffix; several CD4 subsets
   of one study count once) with ≥ 1 colocalisation at PP.H4 > 0.8. `coloc_dat.tsv` is
   filtered to `cell_type == "CD4"`, so other cell types never count. The former
   "eQTL contexts" column (CD4 sub-datasets) is gone.
2. **Match the cite keys above to your .bib**, and confirm that Asthma (EAS) = `sakaue2021`.
2. **Shared samples:** confirm the three groups and decide on IS / MI and the Biobank Japan traits.
3. **eQTL study mapping is unverified.** `blood_cell` suffix → study, in
   `01_coloc_summary_table.py::EQTL_STUDY`: OneK1K → Yazar 2022, Perez → Perez 2022,
   Randolph → Randolph 2021, BP → BLUEPRINT, Sch → DICE (Schmiedel 2018), UT → ImmuNexUT
   (Ota 2021), Bo → Bossini-Castillo 2022.
4. **GWAS N** in `diseases.csv` was copied from the pipeline's `data/GWAS_study_info.csv`.
   Check against the papers.
5. **New schematic PNG** → `chapter2/assets/ch2_schematic.png`.
6. **3,341**: where does it appear? It is probably from the h5ad on the cluster.
7. **SCZ S-LDSC results** are missing locally for this variant (and ALZ Rest / 8 h).
