import math

import numpy as np
import pytest
from scipy import stats
from scipy.special import gammaln

from variant_assoc import (
    allelic_test,
    benjamini_hochberg,
    bonferroni,
    genomic_inflation_factor,
    hwe_exact_pvalue,
)
from variant_assoc.pipeline import run_pipeline
from variant_assoc.simulate import simulate_cohort


def _hwe_brute_force(n_aa, n_ab, n_bb):
    """Enumerate every heterozygote count with the same allele totals (Haldane 1954)."""
    n = n_aa + n_ab + n_bb
    n_a = 2 * n_aa + n_ab
    n_b = 2 * n_bb + n_ab

    def log_prob(het):
        hom_a = (n_a - het) // 2
        hom_b = (n_b - het) // 2
        return (
            gammaln(n + 1) + gammaln(n_a + 1) + gammaln(n_b + 1) + het * math.log(2)
            - gammaln(hom_a + 1) - gammaln(het + 1) - gammaln(hom_b + 1) - gammaln(2 * n + 1)
        )

    hets = [h for h in range(min(n_a, n_b) + 1) if (n_a - h) % 2 == 0]
    logp = np.array([log_prob(h) for h in hets])
    p = np.exp(logp)
    obs = p[hets.index(n_ab)]
    return float(p[p <= obs * (1 + 1e-7)].sum())


@pytest.mark.parametrize(
    "counts",
    [(50, 40, 10), (25, 50, 25), (60, 10, 30), (100, 0, 100), (5, 2, 1), (300, 150, 20), (0, 10, 0)],
)
def test_hwe_matches_brute_force(counts):
    assert hwe_exact_pvalue(*counts) == pytest.approx(_hwe_brute_force(*counts), rel=1e-6)


def test_hwe_perfect_equilibrium_has_high_pvalue():
    assert hwe_exact_pvalue(250, 500, 250) > 0.9


def test_hwe_is_symmetric_in_homozygotes():
    assert hwe_exact_pvalue(70, 20, 10) == pytest.approx(hwe_exact_pvalue(10, 20, 70))


def test_hwe_rejects_negative_counts():
    with pytest.raises(ValueError):
        hwe_exact_pvalue(-1, 5, 5)


def test_allelic_test_matches_manual_fisher():
    case, ctrl = (400, 800, 300), (500, 750, 250)
    res = allelic_test(case, ctrl)
    case_alt, case_ref = 2 * 300 + 800, 2 * 400 + 800
    ctrl_alt, ctrl_ref = 2 * 250 + 750, 2 * 500 + 750
    _, p = stats.fisher_exact([[case_alt, case_ref], [ctrl_alt, ctrl_ref]])
    assert res.p_value == pytest.approx(p)
    assert res.odds_ratio == pytest.approx((case_alt * ctrl_ref) / (case_ref * ctrl_alt))
    assert res.ci_low < res.odds_ratio < res.ci_high


def test_allelic_test_handles_zero_cell():
    res = allelic_test((100, 0, 0), (90, 10, 0))
    assert math.isfinite(res.odds_ratio) and res.ci_low > 0


def test_benjamini_hochberg_matches_scipy():
    rng = np.random.default_rng(0)
    p = np.concatenate([rng.uniform(size=200), rng.uniform(0, 1e-4, size=20)])
    assert benjamini_hochberg(p) == pytest.approx(stats.false_discovery_control(p, method="bh"))


def test_bonferroni_caps_at_one():
    assert bonferroni([0.5, 0.001, 0.2]).tolist() == pytest.approx([1.0, 0.003, 0.6])


def test_genomic_inflation_is_one_for_uniform_pvalues():
    p = (np.arange(1, 10001) - 0.5) / 10000
    assert genomic_inflation_factor(p) == pytest.approx(1.0, abs=0.01)


def test_pipeline_qc_removes_artifacts_and_keeps_signal():
    cohort = simulate_cohort(n_variants=800, n_causal=8, n_artifacts=10, seed=7)
    results, summary = run_pipeline(cohort)
    # Every planted artifact is caught by HWE QC in controls.
    assert summary["artifacts_removed_by_qc"] == summary["n_artifacts_truth"]
    # Naive analysis is polluted by artifacts; QC'd analysis is not.
    assert summary["naive_artifact_hits"] > 0
    assert summary["qc_false_positives"] <= 1
    assert summary["qc_true_positives"] >= 5
    assert summary["lambda_gc_qc"] < summary["lambda_gc_naive"]
