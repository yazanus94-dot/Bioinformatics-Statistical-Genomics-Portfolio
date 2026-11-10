#!/usr/bin/env python3
"""Reproduce the answers in TASK.md, with an independent cross-check."""

import math

import numpy as np
from scipy import stats

from rnaseq_de import log2_fold_change, normalize_counts, size_factors_median_of_ratios

GENES = list("ABCDEF")
COUNTS = np.array(
    [
        [120, 100, 400, 360],
        [800, 900, 1500, 1750],
        [50, 60, 55, 70],
        [0, 3, 40, 50],
        [300, 250, 100, 90],
        [1000, 1200, 2100, 2300],
    ],
    dtype=float,
)
TREATED = np.array([False, False, True, True])

# --- library functions ---
sf = size_factors_median_of_ratios(COUNTS)
lfc = log2_fold_change(normalize_counts(COUNTS, sf), TREATED)

# --- independent re-implementation ---
usable = COUNTS[(COUNTS > 0).all(axis=1)]
geo = np.exp(np.log(usable).mean(axis=1))
sf_check = np.median(usable / geo[:, None], axis=0)
norm_check = COUNTS / sf_check
lfc_check = np.log2((norm_check[:, TREATED].mean(1) + 1) / (norm_check[:, ~TREATED].mean(1) + 1))
assert np.allclose(sf, sf_check) and np.allclose(lfc, lfc_check)

print("a) size factors:", np.round(sf, 4), "-> treat_2 =", f"{sf[3]:.4f}")
top = int(np.argmax(np.abs(lfc)))
print("b) log2FC:", {g: round(float(x), 3) for g, x in zip(GENES, lfc)}, f"-> largest: gene {GENES[top]} = {lfc[top]:.3f}")

N, K, n, k = 12000, 150, 400, 18
p = stats.hypergeom.sf(k - 1, N, K, n)
p_exact = sum(math.comb(K, x) * math.comb(N - K, n - x) for x in range(k, K + 1)) / math.comb(N, n)
assert math.isclose(p, p_exact, rel_tol=1e-9)
print(f"c) p = {p:.3g}, fold enrichment = {k / (n * K / N):.1f}")
