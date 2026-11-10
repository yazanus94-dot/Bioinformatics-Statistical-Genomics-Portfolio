"""Simulate genotype counts for a case-control cohort.

Three kinds of variants are generated:

* **null** variants: same allele frequency in cases and controls;
* **causal** variants: the alternate allele has odds ratio > 1 in cases;
* **artifact** variants: null in truth, but most control hom-alt individuals are
  miscalled as heterozygous (a differential genotyping error). This lowers the
  control alt-allele frequency, creating a spurious association *and* a
  heterozygote excess that Hardy-Weinberg QC in controls can detect.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def _genotype_counts(rng: np.random.Generator, n: int, alt_freq: float) -> np.ndarray:
    probs = [(1 - alt_freq) ** 2, 2 * alt_freq * (1 - alt_freq), alt_freq**2]
    return rng.multinomial(n, probs)


def simulate_cohort(
    n_variants: int = 2000,
    n_cases: int = 1500,
    n_controls: int = 1500,
    n_causal: int = 10,
    n_artifacts: int = 15,
    miscall_rate: float = 0.9,
    seed: int = 42,
) -> pd.DataFrame:
    """Return one row per variant with case/control genotype counts and truth labels."""
    rng = np.random.default_rng(seed)

    # Allele frequencies: skewed towards rare, bounded to keep every variant testable.
    alt_freq = np.clip(rng.beta(0.9, 3.0, size=n_variants), 0.03, 0.5)

    causal_idx = rng.choice(n_variants, size=n_causal, replace=False)
    remaining = np.setdiff1d(np.arange(n_variants), causal_idx)
    artifact_idx = rng.choice(remaining, size=n_artifacts, replace=False)

    causal_or = np.zeros(n_variants)
    causal_or[causal_idx] = rng.uniform(1.4, 1.9, size=n_causal)
    # Artifacts need enough hom-alt individuals to miscall: put them at higher frequency.
    alt_freq[artifact_idx] = rng.uniform(0.30, 0.50, size=n_artifacts)

    rows = []
    for i in range(n_variants):
        p0 = alt_freq[i]
        if causal_or[i] > 0:
            odds = causal_or[i] * p0 / (1 - p0)
            p1 = odds / (1 + odds)
        else:
            p1 = p0

        ctrl = _genotype_counts(rng, n_controls, p0)
        case = _genotype_counts(rng, n_cases, p1)

        is_artifact = i in set(artifact_idx.tolist())
        if is_artifact:
            miscalled = rng.binomial(ctrl[2], miscall_rate)
            ctrl = ctrl + np.array([0, miscalled, -miscalled])

        rows.append(
            {
                "variant_id": f"rs{100000 + i}",
                "case_hom_ref": int(case[0]),
                "case_het": int(case[1]),
                "case_hom_alt": int(case[2]),
                "ctrl_hom_ref": int(ctrl[0]),
                "ctrl_het": int(ctrl[1]),
                "ctrl_hom_alt": int(ctrl[2]),
                "is_causal": bool(causal_or[i] > 0),
                "is_artifact": bool(is_artifact),
            }
        )
    return pd.DataFrame(rows)
