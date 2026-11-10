"""Over-representation analysis (ORA) with the hypergeometric test."""

from __future__ import annotations

from typing import Mapping

import pandas as pd
from scipy import stats


def hypergeometric_enrichment(
    selected: set[str],
    gene_sets: Mapping[str, set[str]],
    universe: set[str],
    min_set_size: int = 5,
) -> pd.DataFrame:
    """Test whether each gene set is over-represented among ``selected`` genes.

    For a set with K genes in a universe of N, and a selected list of n genes of
    which k fall in the set, the p-value is P(X >= k) for a hypergeometric X.
    Both the set and the selected list are restricted to the universe first, which
    matters: using the whole genome as the universe instead of the genes that were
    actually testable inflates significance.
    """
    selected = set(selected) & set(universe)
    n_universe, n_selected = len(universe), len(selected)

    rows = []
    for name, members in gene_sets.items():
        members = set(members) & set(universe)
        if len(members) < min_set_size:
            continue
        overlap = len(members & selected)
        p_value = stats.hypergeom.sf(overlap - 1, n_universe, len(members), n_selected)
        expected = n_selected * len(members) / n_universe
        rows.append(
            {
                "gene_set": name,
                "set_size": len(members),
                "overlap": overlap,
                "expected_overlap": expected,
                "fold_enrichment": overlap / expected if expected > 0 else float("nan"),
                "p_value": p_value,
            }
        )

    table = pd.DataFrame(rows).sort_values("p_value").reset_index(drop=True)
    if not table.empty:
        table["q_value"] = stats.false_discovery_control(table["p_value"].to_numpy(), method="bh")
    return table
