"""Library-size normalisation for RNA-seq count matrices (genes x samples)."""

from __future__ import annotations

import numpy as np


def size_factors_median_of_ratios(counts: np.ndarray) -> np.ndarray:
    """Median-of-ratios size factors (Anders & Huber 2010; used by DESeq2).

    1. For each gene, take the geometric mean of its counts across samples.
    2. For each sample, divide each gene's count by that geometric mean.
    3. The sample's size factor is the median of those ratios.

    Genes with a zero count in *any* sample have an undefined (zero) geometric
    mean and are excluded from the estimate.
    """
    counts = np.asarray(counts, dtype=float)
    if counts.ndim != 2:
        raise ValueError("counts must be a 2-D genes x samples matrix")
    if (counts < 0).any():
        raise ValueError("counts must be non-negative")

    usable = (counts > 0).all(axis=1)
    if not usable.any():
        raise ValueError("no gene has non-zero counts in every sample")

    log_counts = np.log(counts[usable])
    log_geo_mean = log_counts.mean(axis=1, keepdims=True)
    return np.exp(np.median(log_counts - log_geo_mean, axis=0))


def normalize_counts(counts: np.ndarray, size_factors: np.ndarray) -> np.ndarray:
    """Divide each sample (column) by its size factor."""
    counts = np.asarray(counts, dtype=float)
    size_factors = np.asarray(size_factors, dtype=float)
    if counts.shape[1] != size_factors.size:
        raise ValueError("need one size factor per sample (column)")
    return counts / size_factors
