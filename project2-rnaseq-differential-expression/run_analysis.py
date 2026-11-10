#!/usr/bin/env python3
"""Simulate an RNA-seq experiment, find DE genes and test gene-set enrichment."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from rnaseq_de import hypergeometric_enrichment, normalize_counts, welch_de_test
from rnaseq_de.plots import make_figure
from rnaseq_de.qc import pca
from rnaseq_de.simulate import simulate_experiment


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--n-genes", type=int, default=3000)
    parser.add_argument("--n-per-group", type=int, default=6)
    parser.add_argument("--fdr", type=float, default=0.05)
    parser.add_argument("--lfc", type=float, default=1.0, help="minimum |log2 fold change|")
    parser.add_argument("--out", type=Path, default=Path("results"))
    args = parser.parse_args()

    exp = simulate_experiment(n_genes=args.n_genes, n_per_group=args.n_per_group, seed=args.seed)
    de = welch_de_test(exp.counts, exp.is_treated)
    size_factors = de.attrs["size_factors"]

    called = (de["q_value"] < args.fdr) & (de["log2_fc"].abs() >= args.lfc)
    tested_genes = set(de.index[de["q_value"].notna()])
    enrichment = hypergeometric_enrichment(
        selected=set(de.index[called]), gene_sets=exp.gene_sets, universe=tested_genes
    )

    log_norm = np.log2(normalize_counts(exp.counts.to_numpy(), size_factors) + 1)
    scores, var = pca(log_norm)

    truth = exp.truth.loc[de.index]
    true_de = truth["is_de"]
    summary = {
        "n_genes": int(len(de)),
        "n_tested": int(de["q_value"].notna().sum()),
        "n_called_de": int(called.sum()),
        "n_true_de": int(true_de.sum()),
        "true_positives": int((called & true_de).sum()),
        "false_positives": int((called & ~true_de).sum()),
        "sensitivity_among_all_true_de": float((called & true_de).sum() / true_de.sum()),
        "empirical_fdr": float((called & ~true_de).sum() / max(called.sum(), 1)),
        "size_factor_correlation_with_truth": float(
            np.corrcoef(np.log(size_factors), np.log(exp.true_size_factors))[0, 1]
        ),
        "top_enriched_sets": enrichment.head(3)[["gene_set", "overlap", "set_size", "q_value"]]
        .assign(q_value=lambda d: d["q_value"].map(float))
        .to_dict(orient="records"),
    }

    args.out.mkdir(parents=True, exist_ok=True)
    de.join(truth).to_csv(args.out / "de_results.csv", index_label="gene")
    enrichment.to_csv(args.out / "enrichment.csv", index=False)
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    make_figure(scores, var, exp.is_treated, de, enrichment, args.out / "rnaseq_overview.png", args.fdr, args.lfc)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
