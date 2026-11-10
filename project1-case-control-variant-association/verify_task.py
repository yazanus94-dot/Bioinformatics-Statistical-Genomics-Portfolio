#!/usr/bin/env python3
"""Reproduce the answers in TASK.md."""

from variant_assoc import allelic_test, benjamini_hochberg, hwe_exact_pvalue

DATA = {
    "V1": ((430, 450, 120), (490, 420, 90)),
    "V2": ((440, 440, 120), (500, 460, 40)),
    "V3": ((600, 340, 60), (640, 320, 40)),
    "V4": ((355, 485, 160), (360, 480, 160)),
    "V5": ((780, 200, 20), (810, 180, 10)),
    "V6": ((200, 490, 310), (250, 500, 250)),
    "V7": ((700, 270, 30), (700, 270, 30)),
    "V8": ((140, 470, 390), (160, 480, 360)),
}

hwe = {v: hwe_exact_pvalue(*ctrl) for v, (_, ctrl) in DATA.items()}
assoc = {v: allelic_test(case, ctrl) for v, (case, ctrl) in DATA.items()}
kept = [v for v in DATA if hwe[v] >= 1e-3]

print("a) fails QC:", {v: f"{p:.3g}" for v, p in hwe.items() if p < 1e-3})

q_qc = dict(zip(kept, benjamini_hochberg([assoc[v].p_value for v in kept])))
q_all = dict(zip(DATA, benjamini_hochberg([assoc[v].p_value for v in DATA])))
significant = [v for v in kept if q_qc[v] < 0.05]
print("b) significant after QC:", significant, {v: round(float(q_qc[v]), 5) for v in significant})

lead = min(kept, key=lambda v: assoc[v].p_value)
r = assoc[lead]
print(f"c) lead {lead}: p={r.p_value:.3g} OR={r.odds_ratio:.3f} CI=({r.ci_low:.3f}, {r.ci_high:.3f})")
print(f"d) V3 q without QC={q_all['V3']:.4f}, with QC={q_qc['V3']:.4f}")
