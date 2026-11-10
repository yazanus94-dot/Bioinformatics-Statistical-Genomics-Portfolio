# Task: normalise a count matrix, then test gene-set enrichment

**Domain:** RNA-seq analysis
**Skills tested:** median-of-ratios normalisation, fold changes, hypergeometric over-representation test
**Type:** multi-step calculation with verifiable numeric answers

## Problem

A small RNA-seq experiment has 2 control and 2 treated samples and 6 genes.
Raw counts:

| Gene | ctrl_1 | ctrl_2 | treat_1 | treat_2 |
|------|-------:|-------:|--------:|--------:|
| A | 120 | 100 | 400 | 360 |
| B | 800 | 900 | 1500 | 1750 |
| C | 50 | 60 | 55 | 70 |
| D | 0 | 3 | 40 | 50 |
| E | 300 | 250 | 100 | 90 |
| F | 1000 | 1200 | 2100 | 2300 |

1. **Size factors (median of ratios).** For each gene, compute the geometric mean of its counts across the four samples. Genes with a zero count in any sample have no usable geometric mean and are **excluded** from this step. For each sample, the size factor is the median over usable genes of `count / geometric mean`.
2. **Normalise:** divide each sample's counts by its size factor.
3. **Log2 fold change** (treated vs. control) per gene, using group means of the normalised counts and a pseudocount of 1:
   `log2( (mean_treated + 1) / (mean_control + 1) )`. This step uses **all six genes**, including any excluded in step 1.
4. **Enrichment.** The differentially expressed list contains 400 genes out of 12,000 testable genes. A pathway has 150 genes, of which 18 are in the list. Compute the one-sided hypergeometric p-value, P(X >= 18).

**Questions**

a. What is the size factor of `treat_2` (4 decimals)?
b. Which gene has the largest absolute log2 fold change, and what is it (3 decimals)?
c. What is the hypergeometric p-value (3 significant figures), and the fold enrichment?

## Answer

a. Gene D is excluded (it has a zero). Size factors are 0.6823, 0.7734, 1.2793, **1.4823** for `ctrl_1, ctrl_2, treat_1, treat_2`.
b. **Gene D**, log2FC = **3.510**. It is the largest even though it was excluded from normalisation: after normalisation its group means are 1.94 (control) and 32.5 (treated). The next largest in magnitude is gene E at -2.441. The other log2FCs are A 0.860, B 0.010, C -0.729, F 0.082.
c. Expected overlap = 400 x 150 / 12,000 = 5, so the fold enrichment is 18 / 5 = **3.6**, and
   P(X >= 18) = sum over x = 18..150 of C(150, x) C(11850, 400 - x) / C(12000, 400) = **2.48 x 10^-6**.

## Verify

`python verify_task.py` reproduces these numbers using both this repo's functions and an independent plain-NumPy/`math.comb` calculation.
