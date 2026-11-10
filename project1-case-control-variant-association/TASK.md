# Task: HWE quality control changes which variants are significant

**Domain:** statistical genetics / GWAS quality control
**Skills tested:** exact Hardy-Weinberg test, allelic association, odds ratios, multiple-testing correction
**Type:** multi-step calculation with a single verifiable answer

## Problem

A case-control study genotyped 8 biallelic variants in 1,000 cases and 1,000 controls.
Each entry is the genotype count `(hom-ref, het, hom-alt)`.

| Variant | Cases (RR, RA, AA) | Controls (RR, RA, AA) |
|---------|--------------------|-----------------------|
| V1 | 430, 450, 120 | 490, 420, 90 |
| V2 | 440, 440, 120 | 500, 460, 40 |
| V3 | 600, 340, 60 | 640, 320, 40 |
| V4 | 355, 485, 160 | 360, 480, 160 |
| V5 | 780, 200, 20 | 810, 180, 10 |
| V6 | 200, 490, 310 | 250, 500, 250 |
| V7 | 700, 270, 30 | 700, 270, 30 |
| V8 | 140, 470, 390 | 160, 480, 360 |

Apply this pipeline:

1. **QC:** run an *exact* Hardy-Weinberg equilibrium test (Wigginton et al., 2005) in the **controls only**. Exclude every variant with HWE p < 0.001.
2. **Association:** for each remaining variant, run a **two-sided Fisher's exact test** on the 2x2 table of allele counts (alternate allele vs. reference allele, cases vs. controls). Compute the allelic odds ratio (alt allele, cases vs. controls) as `(a*d)/(b*c)`, with a Woolf 95% confidence interval.
3. **Multiple testing:** apply the **Benjamini-Hochberg** procedure **only across the variants that passed QC**. Call a variant significant if its BH-adjusted q < 0.05.

**Questions**

a. Which variant(s) fail QC, and what is the exact HWE p-value of the worst one?
b. How many variants are significant after step 3, and which?
c. For the most significant variant, give the allelic odds ratio and its 95% CI (3 decimals).
d. Without QC (BH across all 8 variants), would V3 be called significant? What is its q-value in each analysis?

## Answer

a. **V2** is the only failure (heterozygote excess in controls); exact HWE p = **6.03 x 10^-8**. All other variants have p >= 0.49.
b. **2 variants** are significant: **V6** (q = 0.00388) and **V1** (q = 0.00910).
c. Most significant is **V6**: Fisher p = 5.55 x 10^-4, OR = **1.247** (95% CI **1.101-1.412**).
d. **Yes without QC, no with QC.** V3's q is 0.0463 across all 8 variants but 0.0540 across the 7 that pass QC. The artifactual V2 (p = 1.8 x 10^-6) inflates the number of small p-values, which loosens the BH thresholds for its neighbours. Removing it restores the honest result.

## Solution outline

```python
from variant_assoc import hwe_exact_pvalue, allelic_test, benjamini_hochberg
```

Allele counts for V6 cases: alt = 2*310 + 490 = 1110, ref = 2*200 + 490 = 890.
Controls: alt = 2*250 + 500 = 1000, ref = 1000. OR = (1110 x 1000) / (890 x 1000) = 1.2472.
SE(log OR) = sqrt(1/1110 + 1/890 + 1/1000 + 1/1000) = 0.0634, so the 95% CI is exp(ln 1.2472 ± 1.96 x 0.0634) = 1.101-1.412.

The script `verify_task.py` in this folder reproduces every number above.
