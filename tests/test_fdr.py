import pathlib
import sys

import numpy as np

PROJECT_ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.append(str(PROJECT_ROOT / "IBDReduce"))

from ibdreduce_v3 import compute_fdr_adjusted_pvalues


def test_fdr_adjustment_matches_manual_example():
    empp = np.array([0.05, 0.2, 0.4], dtype=float)
    permutation_pvalues = np.array([
        [0.8, 0.05, 0.6, 0.1],
        [0.3, 0.25, 0.1, 0.2],
        [0.6, 0.35, 0.2, 0.45],
    ])

    adjusted = compute_fdr_adjusted_pvalues(empp, permutation_pvalues)

    expected = []
    for p in empp:
        rstar = np.sum(permutation_pvalues <= p, axis=0)
        rp = np.sum(empp <= p)
        pm = p * permutation_pvalues.shape[0]
        rb = np.percentile(rstar, 95)
        if rp - rb >= pm:
            expected.append(np.mean(rstar / (rstar + rp - pm)))
        else:
            expected.append(np.sum(rstar >= 1) / len(rstar))

    np.testing.assert_allclose(adjusted, np.array(expected, dtype=float))

    observed_rp = [np.sum(empp <= p) for p in empp]
    permutation_first_column_counts = [np.sum(permutation_pvalues[:, 0] <= p) for p in empp]
    assert observed_rp != permutation_first_column_counts
    assert observed_rp[1] == 2  # sanity-check the hand-calculated rejection counts
