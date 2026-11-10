"""Diagnostic figures for the association study."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

BLUE = "#2a6f97"
GREY = "#9aa5b1"
RED = "#c8553d"


def _qq(ax, p_values: np.ndarray, color: str, label: str) -> None:
    p = np.sort(np.clip(p_values, 1e-300, 1.0))
    expected = -np.log10((np.arange(1, p.size + 1) - 0.5) / p.size)
    ax.scatter(expected, -np.log10(p), s=8, color=color, alpha=0.7, label=label, linewidths=0)


def make_figure(results: pd.DataFrame, summary: dict, path: str | Path) -> None:
    """Save a QQ plot (before/after QC) and a Manhattan-style scatter."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), constrained_layout=True)

    ax = axes[0]
    _qq(ax, results["p_value"].to_numpy(), GREY, f"No QC (λ={summary['lambda_gc_naive']:.2f})")
    passed = results[results["pass_qc"]]
    _qq(ax, passed["p_value"].to_numpy(), BLUE, f"After QC (λ={summary['lambda_gc_qc']:.2f})")
    max_expected = -np.log10(0.5 / len(results))
    ax.plot([0, max_expected], [0, max_expected], color="black", lw=0.8, ls="--")
    ax.set_xlabel("Expected −log10(p)")
    ax.set_ylabel("Observed −log10(p)")
    ax.set_title("QQ plot of allelic association p-values")
    ax.legend(frameon=False)

    ax = axes[1]
    x = np.arange(len(results))
    y = -np.log10(np.clip(results["p_value"].to_numpy(), 1e-300, 1.0))
    ok = results["pass_qc"].to_numpy()
    ax.scatter(x[ok], y[ok], s=7, color=BLUE, alpha=0.6, linewidths=0, label="Pass QC")
    ax.scatter(x[~ok], y[~ok], s=22, color=RED, marker="x", label="Fail QC (HWE/MAF)")
    n_tests = max(int(ok.sum()), 1)
    ax.axhline(-np.log10(0.05 / n_tests), color="black", lw=0.8, ls="--")
    ax.text(len(results) * 0.30, -np.log10(0.05 / n_tests) + 0.25, "Bonferroni threshold (QC'd set)", ha="left", fontsize=8)
    ax.set_xlabel("Variant index")
    ax.set_ylabel("−log10(p)")
    ax.set_title("Association signal; QC failures flagged")
    ax.legend(frameon=False, loc="upper left")

    for a in axes:
        a.spines[["top", "right"]].set_visible(False)

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)
