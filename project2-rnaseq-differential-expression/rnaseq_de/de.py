"""Differential expression testing on normalised counts."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from .normalization import normalize_counts, size_factors_median_of_ratios


def log2_fold_change(
    normalized: np.ndarray, is_treated: np.ndarray, pseudocount: float = 1.0
) -> np.ndarray:
    """log2((mean treated + pc) / (mean control + pc)) on *normalised* counts."""
    normalized = np.asarray(normalized, dtype=float)
    treated = np.asarray(is_treated, dtype=bool)
    mean_t = normalized[:, treated].mean(axis=1)
    mean_c = normalized[:, ~treated].mean(axis=1)
    return np.log2((mean_t + pseudocount) / (mean_c + pseudocount))


def welch_de_test(
    counts: pd.DataFrame,
    is_treated: np.ndarray,
    min_base_mean: float = 5.0,
) -> pd.DataFrame:
    """Differential expression by Welch's t-test on log2(normalised count + 1).

    Genes whose mean normalised count is below ``min_base_mean`` are not tested
    (independent filtering); Benjamini-Hochberg is applied to the tested genes only.

    This is a deliberately simple, transparent model. For production analyses
    with few replicates, prefer a negative-binomial model (DESeq2, edgeR) that
    shares dispersion information across genes.
    """
    is_treated = np.asarray(is_treated, dtype=bool)
    if is_treated.sum() < 2 or (~is_treated).sum() < 2:
        raise ValueError("need at least two samples per group")

    size_factors = size_factors_median_of_ratios(counts.to_numpy())
    norm = normalize_counts(counts.to_numpy(), size_factors)
    log_norm = np.log2(norm + 1.0)

    base_mean = norm.mean(axis=1)
    lfc = log2_fold_change(norm, is_treated)
    t_stat, p_value = stats.ttest_ind(
        log_norm[:, is_treated], log_norm[:, ~is_treated], axis=1, equal_var=False
    )

    result = pd.DataFrame(
        {"base_mean": base_mean, "log2_fc": lfc, "t_stat": t_stat, "p_value": p_value},
        index=counts.index,
    )
    tested = (base_mean >= min_base_mean) & np.isfinite(p_value)
    result.loc[~tested, ["t_stat", "p_value"]] = np.nan
    result["q_value"] = np.nan
    result.loc[tested, "q_value"] = stats.false_discovery_control(
        result.loc[tested, "p_value"].to_numpy(), method="bh"
    )
    result.attrs["size_factors"] = size_factors
    return result
