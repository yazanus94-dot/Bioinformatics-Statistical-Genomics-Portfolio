"""Statistical building blocks for a case-control variant association study.

Genotype counts are always given as ``(n_hom_ref, n_het, n_hom_alt)``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

import numpy as np
from scipy import stats


@dataclass(frozen=True)
class AlleleTestResult:
    """Result of an allelic (2x2 allele-count) association test."""

    odds_ratio: float
    ci_low: float
    ci_high: float
    p_value: float


def hwe_exact_pvalue(n_hom_ref: int, n_het: int, n_hom_alt: int) -> float:
    """Exact test for Hardy-Weinberg equilibrium.

    Implements the O(n) algorithm of Wigginton, Cutler & Abecasis (2005,
    Am. J. Hum. Genet. 76:887-893). The p-value is the total probability of
    all heterozygote counts that are no more likely than the observed one,
    conditional on the observed allele counts.
    """
    if min(n_hom_ref, n_het, n_hom_alt) < 0:
        raise ValueError("genotype counts must be non-negative")
    n = n_hom_ref + n_het + n_hom_alt
    if n == 0:
        return 1.0

    hom_common = max(n_hom_ref, n_hom_alt)
    hom_rare = min(n_hom_ref, n_hom_alt)
    rare_alleles = 2 * hom_rare + n_het

    probs = np.zeros(rare_alleles + 1)

    # Start at the most likely heterozygote count (same parity as rare_alleles).
    mid = rare_alleles * (2 * n - rare_alleles) // (2 * n)
    if (mid % 2) != (rare_alleles % 2):
        mid += 1

    cur_hets = mid
    cur_rare = (rare_alleles - mid) // 2
    cur_common = n - cur_hets - cur_rare
    probs[mid] = 1.0
    total = 1.0

    # Walk down in steps of two heterozygotes.
    while cur_hets > 1:
        probs[cur_hets - 2] = (
            probs[cur_hets] * cur_hets * (cur_hets - 1.0) / (4.0 * (cur_rare + 1.0) * (cur_common + 1.0))
        )
        total += probs[cur_hets - 2]
        cur_hets -= 2
        cur_rare += 1
        cur_common += 1

    # Walk up in steps of two heterozygotes.
    cur_hets = mid
    cur_rare = (rare_alleles - mid) // 2
    cur_common = n - cur_hets - cur_rare
    while cur_hets <= rare_alleles - 2:
        probs[cur_hets + 2] = (
            probs[cur_hets] * 4.0 * cur_rare * cur_common / ((cur_hets + 2.0) * (cur_hets + 1.0))
        )
        total += probs[cur_hets + 2]
        cur_hets += 2
        cur_rare -= 1
        cur_common -= 1

    probs /= total
    tol = probs[n_het] * (1.0 + 1e-7)
    return float(min(1.0, probs[probs <= tol].sum()))


def allele_counts(genotypes: Sequence[int]) -> tuple[int, int]:
    """Return ``(ref_alleles, alt_alleles)`` from ``(hom_ref, het, hom_alt)`` counts."""
    hom_ref, het, hom_alt = genotypes
    return 2 * hom_ref + het, 2 * hom_alt + het


def minor_allele_frequency(genotypes: Sequence[int]) -> float:
    """Minor allele frequency from genotype counts."""
    ref, alt = allele_counts(genotypes)
    total = ref + alt
    return min(ref, alt) / total if total else 0.0


def allelic_test(
    case_genotypes: Sequence[int], control_genotypes: Sequence[int], alpha: float = 0.05
) -> AlleleTestResult:
    """Allelic association test (alt allele, cases vs controls).

    The p-value is from Fisher's exact test on the 2x2 allele-count table.
    The odds ratio is the sample odds ratio with a Woolf (log-scale) Wald
    confidence interval; a Haldane-Anscombe 0.5 correction is applied only if
    a cell is zero. Treating the two alleles of an individual as independent
    assumes Hardy-Weinberg equilibrium, which is why HWE QC precedes this test.
    """
    case_ref, case_alt = allele_counts(case_genotypes)
    ctrl_ref, ctrl_alt = allele_counts(control_genotypes)

    _, p_value = stats.fisher_exact([[case_alt, case_ref], [ctrl_alt, ctrl_ref]])

    cells = np.array([case_alt, case_ref, ctrl_alt, ctrl_ref], dtype=float)
    if (cells == 0).any():
        cells = cells + 0.5
    a, b, c, d = cells
    odds_ratio = (a * d) / (b * c)
    se = math.sqrt(1.0 / a + 1.0 / b + 1.0 / c + 1.0 / d)
    z = stats.norm.ppf(1.0 - alpha / 2.0)
    return AlleleTestResult(
        odds_ratio=float(odds_ratio),
        ci_low=float(math.exp(math.log(odds_ratio) - z * se)),
        ci_high=float(math.exp(math.log(odds_ratio) + z * se)),
        p_value=float(p_value),
    )


def bonferroni(p_values: Sequence[float]) -> np.ndarray:
    """Bonferroni-adjusted p-values (capped at 1)."""
    p = np.asarray(p_values, dtype=float)
    return np.minimum(p * p.size, 1.0)


def benjamini_hochberg(p_values: Sequence[float]) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (q-values controlling the FDR)."""
    p = np.asarray(p_values, dtype=float)
    n = p.size
    if n == 0:
        return p.copy()
    order = np.argsort(p)
    scaled = p[order] * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(scaled[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(adjusted, 1.0)
    return out


def genomic_inflation_factor(p_values: Sequence[float]) -> float:
    """Genomic inflation factor lambda_GC = median(chi2) / expected median.

    Values near 1 indicate a well-calibrated test; values well above 1 suggest
    population stratification, batch effects or systematic genotyping error.
    """
    p = np.clip(np.asarray(p_values, dtype=float), 1e-300, 1.0)
    chisq = stats.chi2.isf(p, df=1)
    return float(np.median(chisq) / stats.chi2.ppf(0.5, df=1))
