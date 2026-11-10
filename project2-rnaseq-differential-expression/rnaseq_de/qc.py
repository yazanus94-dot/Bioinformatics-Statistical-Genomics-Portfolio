"""Sample-level quality control."""

from __future__ import annotations

import numpy as np


def pca(log_expression: np.ndarray, n_components: int = 2) -> tuple[np.ndarray, np.ndarray]:
    """PCA of samples (columns) from a genes x samples log-expression matrix.

    Returns ``(scores, explained_variance_ratio)`` where scores has shape
    ``(n_samples, n_components)``. Genes are mean-centred before the SVD.
    """
    x = np.asarray(log_expression, dtype=float)
    centred = (x - x.mean(axis=1, keepdims=True)).T
    u, s, _ = np.linalg.svd(centred, full_matrices=False)
    variance = s**2 / (centred.shape[0] - 1)
    return (u * s)[:, :n_components], (variance / variance.sum())[:n_components]
