"""Simulate an RNA-seq experiment with known differentially expressed genes.

Counts follow a negative binomial distribution with a mean-dependent dispersion
and sample-specific library sizes. Two gene sets are planted as enriched: most
of ``set_00`` is up-regulated and most of ``set_07`` is down-regulated, on top of
background DE genes scattered across the genome.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass
class SimulatedExperiment:
    counts: pd.DataFrame  # genes x samples
    is_treated: np.ndarray  # bool per sample
    gene_sets: dict[str, set[str]]
    truth: pd.DataFrame  # per gene: true_log2_fc, is_de
    true_size_factors: np.ndarray


def simulate_experiment(
    n_genes: int = 3000,
    n_per_group: int = 6,
    n_sets: int = 30,
    set_size: int = 60,
    seed: int = 1,
) -> SimulatedExperiment:
    rng = np.random.default_rng(seed)
    genes = [f"GENE{i:04d}" for i in range(n_genes)]
    n_samples = 2 * n_per_group
    is_treated = np.array([False] * n_per_group + [True] * n_per_group)

    gene_sets = {
        f"set_{k:02d}": set(genes[k * set_size : (k + 1) * set_size]) for k in range(n_sets)
    }

    true_lfc = np.zeros(n_genes)
    up = rng.choice(np.arange(0, set_size), size=40, replace=False)  # inside set_00
    down = 7 * set_size + rng.choice(np.arange(0, set_size), size=30, replace=False)  # set_07
    background_pool = np.arange(n_sets * set_size, n_genes)
    background = rng.choice(background_pool, size=50, replace=False)
    true_lfc[up] = rng.uniform(1.0, 2.5, size=up.size)
    true_lfc[down] = -rng.uniform(1.0, 2.5, size=down.size)
    true_lfc[background] = rng.choice([-1, 1], size=background.size) * rng.uniform(1.0, 2.5, size=background.size)

    base_mean = np.exp(rng.normal(np.log(150), 1.3, size=n_genes))
    dispersion = 0.1 + 1.0 / base_mean
    size_factors = np.exp(rng.normal(0.0, 0.2, size=n_samples))

    fold = 2.0 ** (np.outer(true_lfc, is_treated.astype(float)))
    mu = base_mean[:, None] * size_factors[None, :] * fold
    phi = dispersion[:, None]
    counts = rng.negative_binomial(n=1.0 / phi, p=1.0 / (1.0 + phi * mu))

    sample_names = [f"ctrl_{i + 1}" for i in range(n_per_group)] + [
        f"treat_{i + 1}" for i in range(n_per_group)
    ]
    return SimulatedExperiment(
        counts=pd.DataFrame(counts, index=genes, columns=sample_names),
        is_treated=is_treated,
        gene_sets=gene_sets,
        truth=pd.DataFrame({"true_log2_fc": true_lfc, "is_de": true_lfc != 0}, index=genes),
        true_size_factors=size_factors,
    )
