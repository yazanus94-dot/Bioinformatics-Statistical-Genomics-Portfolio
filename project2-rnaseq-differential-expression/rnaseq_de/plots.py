"""Figures for the RNA-seq analysis."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

BLUE = "#2a6f97"
RED = "#c8553d"
GREY = "#b8c0c8"


def make_figure(
    pca_scores: np.ndarray,
    pca_var: np.ndarray,
    is_treated: np.ndarray,
    de: pd.DataFrame,
    enrichment: pd.DataFrame,
    path: str | Path,
    fdr: float = 0.05,
    lfc_cut: float = 1.0,
    top_sets: int = 6,
) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4), constrained_layout=True)

    ax = axes[0]
    ax.scatter(*pca_scores[~is_treated].T, color=GREY, edgecolor="black", lw=0.5, s=55, label="Control")
    ax.scatter(*pca_scores[is_treated].T, color=BLUE, edgecolor="black", lw=0.5, s=55, label="Treated")
    ax.set_xlabel(f"PC1 ({pca_var[0]:.0%} variance)")
    ax.set_ylabel(f"PC2 ({pca_var[1]:.0%} variance)")
    ax.set_title("Samples (PCA of log2 normalised counts)")
    ax.legend(frameon=False)

    ax = axes[1]
    tested = de.dropna(subset=["q_value"])
    hit = (tested["q_value"] < fdr) & (tested["log2_fc"].abs() >= lfc_cut)
    y = -np.log10(np.clip(tested["p_value"], 1e-300, 1))
    ax.scatter(tested.loc[~hit, "log2_fc"], y[~hit], s=6, color=GREY, linewidths=0)
    ax.scatter(tested.loc[hit, "log2_fc"], y[hit], s=10, color=RED, linewidths=0, label=f"{int(hit.sum())} DE genes")
    ax.axvline(lfc_cut, color="black", lw=0.6, ls="--")
    ax.axvline(-lfc_cut, color="black", lw=0.6, ls="--")
    ax.set_xlabel("log2 fold change (treated / control)")
    ax.set_ylabel("−log10(p)")
    ax.set_title(f"Volcano plot (FDR < {fdr}, |log2FC| ≥ {lfc_cut:g})")
    ax.legend(frameon=False, loc="upper center")

    ax = axes[2]
    top = enrichment.head(top_sets).iloc[::-1]
    ax.barh(top["gene_set"], -np.log10(top["q_value"]), color=BLUE)
    ax.axvline(-np.log10(fdr), color="black", lw=0.8, ls="--")
    ax.set_xlabel("−log10(q)")
    ax.set_title("Top enriched gene sets (hypergeometric)")

    for a in axes:
        a.spines[["top", "right"]].set_visible(False)

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)
