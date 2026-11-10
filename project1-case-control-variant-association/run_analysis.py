#!/usr/bin/env python3
"""Run the simulated case-control association study and write results to ./results."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from variant_assoc.pipeline import run_pipeline
from variant_assoc.plots import make_figure
from variant_assoc.simulate import simulate_cohort


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--n-variants", type=int, default=2000)
    parser.add_argument("--n-cases", type=int, default=1500)
    parser.add_argument("--n-controls", type=int, default=1500)
    parser.add_argument("--fdr", type=float, default=0.05)
    parser.add_argument("--hwe-threshold", type=float, default=1e-6)
    parser.add_argument("--maf-threshold", type=float, default=0.01)
    parser.add_argument("--out", type=Path, default=Path("results"))
    args = parser.parse_args()

    cohort = simulate_cohort(
        n_variants=args.n_variants,
        n_cases=args.n_cases,
        n_controls=args.n_controls,
        seed=args.seed,
    )
    results, summary = run_pipeline(
        cohort, hwe_threshold=args.hwe_threshold, maf_threshold=args.maf_threshold, fdr=args.fdr
    )

    args.out.mkdir(parents=True, exist_ok=True)
    results.to_csv(args.out / "variant_results.csv", index=False)
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    make_figure(results, summary, args.out / "association_plots.png")

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
