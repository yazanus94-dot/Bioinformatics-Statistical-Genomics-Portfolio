"""End-to-end association pipeline: QC -> allelic tests -> multiple-testing control."""

from __future__ import annotations

import numpy as np
import pandas as pd

from .stats import (
    allelic_test,
    benjamini_hochberg,
    bonferroni,
    genomic_inflation_factor,
    hwe_exact_pvalue,
    minor_allele_frequency,
)

CASE_COLS = ["case_hom_ref", "case_het", "case_hom_alt"]
CTRL_COLS = ["ctrl_hom_ref", "ctrl_het", "ctrl_hom_alt"]


def run_pipeline(
    df: pd.DataFrame,
    hwe_threshold: float = 1e-6,
    maf_threshold: float = 0.01,
    fdr: float = 0.05,
) -> tuple[pd.DataFrame, dict]:
    """Annotate every variant with QC flags, association statistics and adjusted p-values.

    Association is run on *all* variants so that the effect of QC can be shown, but
    adjusted p-values and discoveries are reported both for the naive analysis
    (no QC) and for the QC-filtered analysis.
    """
    out = df.copy()

    out["hwe_p_controls"] = [hwe_exact_pvalue(*row) for row in out[CTRL_COLS].itertuples(index=False)]
    pooled = out[CASE_COLS].to_numpy() + out[CTRL_COLS].to_numpy()
    out["maf"] = [minor_allele_frequency(g) for g in pooled]
    out["pass_qc"] = (out["hwe_p_controls"] >= hwe_threshold) & (out["maf"] >= maf_threshold)

    results = [
        allelic_test(case, ctrl)
        for case, ctrl in zip(
            out[CASE_COLS].itertuples(index=False), out[CTRL_COLS].itertuples(index=False)
        )
    ]
    out["odds_ratio"] = [r.odds_ratio for r in results]
    out["ci_low"] = [r.ci_low for r in results]
    out["ci_high"] = [r.ci_high for r in results]
    out["p_value"] = [r.p_value for r in results]

    # Naive: correct over every variant. QC'd: correct over variants that pass QC only.
    out["q_bh_naive"] = benjamini_hochberg(out["p_value"])
    out["q_bonf_naive"] = bonferroni(out["p_value"])
    out["q_bh_qc"] = np.nan
    out["q_bonf_qc"] = np.nan
    passed = out["pass_qc"].to_numpy()
    out.loc[passed, "q_bh_qc"] = benjamini_hochberg(out.loc[passed, "p_value"])
    out.loc[passed, "q_bonf_qc"] = bonferroni(out.loc[passed, "p_value"])

    naive_hits = out["q_bh_naive"] < fdr
    qc_hits = out["q_bh_qc"] < fdr

    summary = {
        "n_variants": int(len(out)),
        "n_fail_hwe": int((out["hwe_p_controls"] < hwe_threshold).sum()),
        "n_fail_maf": int((out["maf"] < maf_threshold).sum()),
        "n_pass_qc": int(passed.sum()),
        "fdr_threshold": fdr,
        "lambda_gc_naive": genomic_inflation_factor(out["p_value"]),
        "lambda_gc_qc": genomic_inflation_factor(out.loc[passed, "p_value"]),
        "naive_discoveries": int(naive_hits.sum()),
        "qc_discoveries": int(qc_hits.sum()),
        "bonferroni_discoveries_qc": int((out["q_bonf_qc"] < fdr).sum()),
    }
    if {"is_causal", "is_artifact"} <= set(out.columns):
        summary.update(
            {
                "n_causal_truth": int(out["is_causal"].sum()),
                "naive_true_positives": int((naive_hits & out["is_causal"]).sum()),
                "naive_false_positives": int((naive_hits & ~out["is_causal"]).sum()),
                "naive_artifact_hits": int((naive_hits & out["is_artifact"]).sum()),
                "qc_true_positives": int((qc_hits & out["is_causal"]).sum()),
                "qc_false_positives": int((qc_hits & ~out["is_causal"]).sum()),
                "artifacts_removed_by_qc": int((out["is_artifact"] & ~out["pass_qc"]).sum()),
                "n_artifacts_truth": int(out["is_artifact"].sum()),
            }
        )
    return out, summary
