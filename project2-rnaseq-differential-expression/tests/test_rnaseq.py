import math

import numpy as np
import pandas as pd
import pytest

from rnaseq_de import (
    hypergeometric_enrichment,
    log2_fold_change,
    normalize_counts,
    size_factors_median_of_ratios,
    welch_de_test,
)
from rnaseq_de.simulate import simulate_experiment


def test_size_factors_hand_computed_example():
    # Gene 4 has a zero and must be ignored. Sample 2 is exactly 2x sample 1.
    counts = np.array([[10, 20], [30, 60], [100, 200], [0, 5]])
    sf = size_factors_median_of_ratios(counts)
    assert sf[1] / sf[0] == pytest.approx(2.0)
    # geometric mean of sf is 1 by construction for proportional samples
    assert math.sqrt(sf[0] * sf[1]) == pytest.approx(1.0)


def test_size_factors_scale_with_library_size():
    rng = np.random.default_rng(0)
    base = rng.integers(20, 500, size=(300, 4)).astype(float)
    scaled = base.copy()
    scaled[:, 2] *= 3
    ratio = size_factors_median_of_ratios(scaled) / size_factors_median_of_ratios(base)
    assert ratio[2] / ratio[0] == pytest.approx(3.0, rel=1e-9)


def test_normalisation_removes_depth_difference():
    counts = np.array([[10, 40], [50, 200], [20, 80], [5, 20]])
    norm = normalize_counts(counts, size_factors_median_of_ratios(counts))
    assert norm[:, 0] == pytest.approx(norm[:, 1])


def test_size_factors_need_a_gene_without_zeros():
    with pytest.raises(ValueError):
        size_factors_median_of_ratios(np.array([[0, 1], [2, 0]]))


def test_log2_fold_change_with_pseudocount():
    norm = np.array([[7.0, 7.0, 15.0, 15.0]])
    lfc = log2_fold_change(norm, np.array([False, False, True, True]), pseudocount=1.0)
    assert lfc[0] == pytest.approx(math.log2(16 / 8))


def test_hypergeometric_matches_direct_combinatorics():
    universe = {f"g{i}" for i in range(100)}
    gene_set = {f"g{i}" for i in range(10)}
    selected = {f"g{i}" for i in range(5)} | {f"g{i}" for i in range(50, 60)}
    table = hypergeometric_enrichment(selected, {"s": gene_set}, universe, min_set_size=1)
    # N=100, K=10, n=15, k=5 -> P(X >= 5)
    expected = sum(
        math.comb(10, x) * math.comb(90, 15 - x) for x in range(5, 11)
    ) / math.comb(100, 15)
    assert table.loc[0, "p_value"] == pytest.approx(expected, rel=1e-9)
    assert table.loc[0, "overlap"] == 5


def test_enrichment_respects_the_universe():
    universe = {f"g{i}" for i in range(50)}
    outside = {"x1", "x2", "x3"}
    table = hypergeometric_enrichment({"g0", "g1"} | outside, {"s": {"g0", "g1", "g2"} | outside}, universe, 1)
    assert table.loc[0, "set_size"] == 3  # genes outside the universe are dropped


def test_de_pipeline_recovers_planted_signal():
    exp = simulate_experiment(seed=3)
    de = welch_de_test(exp.counts, exp.is_treated)
    called = (de["q_value"] < 0.05) & (de["log2_fc"].abs() >= 1.0)
    truth = exp.truth.loc[de.index, "is_de"]
    assert (called & truth).sum() >= 0.6 * truth.sum()
    assert (called & ~truth).sum() <= 0.1 * max(called.sum(), 1)
    # Estimated size factors track the simulated library sizes.
    corr = np.corrcoef(np.log(de.attrs["size_factors"]), np.log(exp.true_size_factors))[0, 1]
    assert corr > 0.95


def test_planted_gene_sets_are_the_top_hits():
    exp = simulate_experiment(seed=3)
    de = welch_de_test(exp.counts, exp.is_treated)
    called = set(de.index[(de["q_value"] < 0.05) & (de["log2_fc"].abs() >= 1.0)])
    table = hypergeometric_enrichment(called, exp.gene_sets, set(de.index[de["q_value"].notna()]))
    assert set(table["gene_set"].head(2)) == {"set_00", "set_07"}
    assert (table["q_value"].head(2) < 0.05).all()


def test_welch_requires_replicates():
    counts = pd.DataFrame(np.ones((5, 3)), columns=["a", "b", "c"])
    with pytest.raises(ValueError):
        welch_de_test(counts, np.array([False, True, True]))
