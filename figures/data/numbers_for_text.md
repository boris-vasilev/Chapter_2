# Numbers for the text — Chapter 2

## Heritability

S-LDSC, variant coloc_thr0.5 (baselineLD v2.2 + anchors, regulators and targets fitted jointly per trait x condition). 129 fitted sets of 162 trait x condition x annotation cells. Significant = BH FDR < 0.05, corrected within expanded sets (regulators + targets, 84 tests) and within anchors (45 tests) separately, for enrichment (two-sided) and for τ* (one-sided, H1: τ* > 0).

Flagged, never counted as significant (negative Prop._h2): Systemic lupus erythematosus Rest regulators (Prop._h2 = -0.0047, enrichment = -11.8).

### Sets significant, by annotation and condition (significant / fitted)

| Annotation | Condition | Enrichment, IMD | τ*, IMD | Enrichment, control | τ*, control |
|---|---|---|---|---|---|
| Anchors | Rest | 7/12 | 7/12 | 2/3 | 2/3 |
| Anchors | 8 h | 7/12 | 7/12 | 2/3 | 2/3 |
| Anchors | 48 h | 7/12 | 7/12 | 2/3 | 2/3 |
| Regulators | Rest | 7/12 | 0/12 | 0/4 | 0/4 |
| Regulators | 8 h | 8/12 | 0/12 | 1/4 | 0/4 |
| Regulators | 48 h | 10/12 | 0/12 | 0/4 | 0/4 |
| Targets | Rest | 6/10 | 0/10 | 0/1 | 0/1 |
| Targets | 8 h | 8/10 | 0/10 | 1/2 | 0/2 |
| Targets | 48 h | 8/10 | 0/10 | 0/3 | 0/3 |

### Median log2 enrichment (all conditions pooled)

| Annotation | IMD | Control |
|---|---|---|
| Anchors | 6.28 (n = 36; median enrichment 77.56) | 6.08 (n = 9; median enrichment 67.84) |
| Regulators | 2.39 (n = 35; median enrichment 5.21) | 0.45 (n = 12; median enrichment 1.37) |
| Targets | 1.71 (n = 30; median enrichment 3.28) | 1.56 (n = 6; median enrichment 2.97) |

log2 medians exclude sets with enrichment ≤ 0 (Systemic lupus erythematosus Rest regulators = -11.8); the plain median includes them.

### Marginal vs conditional

| Annotation | Fitted | Both | Enrichment only | τ* only | Neither |
|---|---|---|---|---|---|
| Anchors | 45 | 27 | 0 | 0 | 18 |
| Regulators | 48 | 0 | 26 | 0 | 22 |
| Targets | 36 | 0 | 23 | 0 | 13 |

- Expanded sets (regulators + targets), immune-mediated: 47/66 enriched, 0/66 τ*-significant, 47 enriched but not τ*.
- Expanded sets (regulators + targets), non-IMD control: 2/18 enriched, 0/18 τ*-significant, 2 enriched but not τ*.

### Anchors per trait (reference point)

Anchor estimates differ slightly between conditions only because the joint model's other annotations change; Rest shown, range over conditions in brackets.

| Trait | Class | Anchors | Enrichment ± SE | p | FDR < 0.05 | τ* ± SE | p | FDR < 0.05 |
|---|---|---|---|---|---|---|---|---|
| Autoimmune thyroid disease | IMD | 63 | 27.5 ± 6.8 [27.2–27.5] | 1.3e-04 | yes | 1.27 ± 0.34 [1.24–1.27] | 1.1e-04 | yes |
| Adult-onset asthma | IMD | 23 | 101.0 ± 18.4 [101.0–101.6] | 1.2e-08 | yes | 2.29 ± 0.39 [2.29–2.30] | 2.9e-09 | yes |
| Childhood-onset asthma | IMD | 27 | 75.7 ± 17.2 [75.7–75.8] | 5.0e-05 | yes | 2.11 ± 0.52 [2.09–2.11] | 2.5e-05 | yes |
| Asthma (EAS) | IMD | 9 | 95.5 ± 57.2 [95.5–96.1] | 9.9e-02 | no | 1.29 ± 0.79 [1.29–1.30] | 5.1e-02 | no |
| Atopic dermatitis | IMD | 18 | 122.9 ± 36.9 [122.8–123.2] | 1.0e-03 | yes | 2.47 ± 0.76 [2.47–2.48] | 6.0e-04 | yes |
| Crohn's disease | IMD | 28 | 54.4 ± 32.9 [54.3–54.5] | 1.2e-01 | no | 1.57 ± 1.03 [1.57–1.58] | 6.3e-02 | no |
| Inflammatory bowel disease | IMD | 24 | 70.4 ± 20.8 [70.2–70.5] | 8.1e-04 | yes | 1.68 ± 0.52 [1.67–1.68] | 5.7e-04 | yes |
| Ulcerative colitis | IMD | 9 | 123.7 ± 126.2 [123.7–123.8] | 3.4e-01 | no | 1.14 ± 1.24 [1.14–1.14] | 1.8e-01 | no |
| Rheumatoid arthritis (EAS) | IMD | 16 | 79.3 ± 52.1 [79.3–79.6] | 1.5e-01 | no | 1.68 ± 1.19 [1.68–1.68] | 8.0e-02 | no |
| Rheumatoid arthritis (EUR) | IMD | 25 | 45.7 ± 14.5 [45.6–45.7] | 1.9e-03 | yes | 1.08 ± 0.36 [1.08–1.08] | 1.3e-03 | yes |
| Systemic lupus erythematosus | IMD | 7 | 264.4 ± 286.2 [263.5–264.4] | 3.8e-01 | no | 3.24 ± 3.77 [3.23–3.24] | 1.9e-01 | no |
| Type 1 diabetes | IMD | 26 | 21.9 ± 9.1 [21.6–21.9] | 6.7e-03 | yes | 0.57 ± 0.22 [0.56–0.57] | 4.6e-03 | yes |
| Ischaemic stroke | control | 3 | 173.6 ± 338.9 [173.5–174.0] | 6.2e-01 | no | 1.83 ± 3.85 [1.83–1.84] | 3.2e-01 | no |
| Myocardial infarction | control | 20 | 67.9 ± 26.7 [67.8–67.9] | 1.2e-02 | yes | 1.45 ± 0.57 [1.44–1.45] | 6.0e-03 | yes |
| Type 2 diabetes | control | 21 | 29.6 ± 9.7 [29.6–29.7] | 3.1e-03 | yes | 0.69 ± 0.24 [0.69–0.69] | 1.8e-03 | yes |

### Missing runs and empty sets

- **dropped (< 1,000 SNPs)** (5): Alzheimer's disease (proxy) 48 h anchors (n = 3); Alzheimer's disease (proxy) 48 h regulators (n = 2); Osteoarthritis Rest anchors (n = 5); Osteoarthritis 8 h anchors (n = 5); Osteoarthritis 48 h anchors (n = 5)
- **empty set** (15): Asthma (EAS) Rest targets (n = 0); Asthma (EAS) 8 h targets (n = 0); Asthma (EAS) 48 h targets (n = 0); Ulcerative colitis Rest targets (n = 0); Ulcerative colitis 8 h targets (n = 0); Ulcerative colitis 48 h targets (n = 0); Alzheimer's disease (proxy) Rest targets (n = 0); Alzheimer's disease (proxy) 8 h targets (n = 0); Ischaemic stroke Rest targets (n = 0); Ischaemic stroke 8 h targets (n = 0); Ischaemic stroke 48 h targets (n = 0); Osteoarthritis Rest targets (n = 0); Osteoarthritis 8 h targets (n = 0); Osteoarthritis 48 h targets (n = 0); Type 2 diabetes Rest targets (n = 0)
- **not run** (4): Alzheimer's disease (proxy) Rest anchors (n = 3); Alzheimer's disease (proxy) Rest regulators (n = 3); Alzheimer's disease (proxy) 8 h anchors (n = 3); Alzheimer's disease (proxy) 8 h regulators (n = 5)
- **run failed** (9): Schizophrenia Rest anchors (n = 17); Schizophrenia Rest regulators (n = 81); Schizophrenia Rest targets (n = 25); Schizophrenia 8 h anchors (n = 17); Schizophrenia 8 h regulators (n = 109); Schizophrenia 8 h targets (n = 32); Schizophrenia 48 h anchors (n = 17); Schizophrenia 48 h regulators (n = 90); Schizophrenia 48 h targets (n = 145)
