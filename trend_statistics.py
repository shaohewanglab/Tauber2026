"""Two-sided Spearman permutation trends for independent model runs or glands."""
from __future__ import annotations

import itertools
import math

import numpy as np
from scipy.stats import rankdata

DEFAULT_RESAMPLES = 999_999
DEFAULT_SEED = 20261018


def permutation_spearman(x, y, *, n_resamples=DEFAULT_RESAMPLES,
                         seed=DEFAULT_SEED, batch=10_000):
    """Permute outcome ranks against fixed predictor ranks, preserving ties.

    The two-sided statistic is abs(rho). Enumerate all permutations for n <= 8
    when the requested budget permits; otherwise use independent Monte Carlo
    permutations with (extreme + 1)/(B + 1), never a reported zero P value.
    """
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if x.ndim != 1 or y.ndim != 1 or len(x) != len(y) or len(x) < 3:
        raise ValueError("Provide equally sized 1D arrays with at least three independent observations.")
    if not (np.isfinite(x).all() and np.isfinite(y).all()):
        raise ValueError("Trend observations must be finite; missing values are not silently removed.")
    if n_resamples < 1 or batch < 1:
        raise ValueError("Permutation count and batch size must be positive.")
    a, b = rankdata(x), rankdata(y)
    a -= a.mean()
    b -= b.mean()
    norm = np.linalg.norm(a) * np.linalg.norm(b)
    if norm == 0:
        raise ValueError("Spearman correlation requires variation in both variables.")
    rho = float(a @ b / norm)
    threshold = abs(rho) - 1e-14  # Count numerical ties in the observed tail.
    exact = len(x) <= 8 and math.factorial(len(x)) <= n_resamples
    extreme = 0
    if exact:
        draws = math.factorial(len(x))
        for permutation in itertools.permutations(b):
            extreme += abs(np.dot(a, permutation) / norm) >= threshold
        p = extreme / draws
    else:
        draws = int(n_resamples)
        rng = np.random.default_rng(seed)
        for start in range(0, draws, batch):
            size = min(batch, draws - start)
            permutations = rng.permuted(np.broadcast_to(b, (size, len(b))), axis=1)
            extreme += int(np.count_nonzero(np.abs(permutations @ a / norm) >= threshold))
        p = (extreme + 1) / (draws + 1)
    return {"method": "Spearman rank correlation with two-sided permutation test",
            "null_hypothesis": "Predictor and outcome are independent; outcomes are exchangeable across predictor values",
            "rho": rho, "n": len(x), "p_raw": float(p),
            "permutation_mode": "exhaustive" if exact else "Monte Carlo",
            "permutations": draws, "extreme_permutations": int(extreme),
            "seed": None if exact else int(seed),
            "two_sided_convention": "Absolute Spearman rho at least as large as observed",
            "monte_carlo_correction": None if exact else "(extreme + 1)/(permutations + 1)",
            "predictor_ties_preserved": True, "outcome_ties_preserved": True,
            "statistical_unit": "One independently simulated run or one gland; never a pooled count treated as multiple replicates"}
