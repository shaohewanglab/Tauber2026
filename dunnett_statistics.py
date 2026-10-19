"""Two-sided, single-step Dunnett comparisons of means against one control.

The multivariate-t tail is integrated deterministically, conditional on the
shared control normal variate and independent pooled chi-square variate.
Directly integrating the union tail avoids subtracting a CDF close to one.
This matters for the very small probabilities in the manuscript data.
"""
from functools import lru_cache

import numpy as np
from scipy import integrate, special, stats

ASSUMPTIONS = ("Independent observations within and between groups; approximately "
               "normal within-group errors with one common variance. Dunnett's "
               "single-step adjustment controls the family of treatment-versus-control "
               "comparisons; no additional Holm correction is applied.")


@lru_cache(maxsize=256)
def _adjusted_tail(threshold, sample_sizes, control_size, df, nodes=256):
    if threshold == 0:
        return 1.0
    if len(sample_sizes) == 1:
        return float(2*stats.t.sf(threshold, df))
    normal, weights = special.roots_hermitenorm(nodes)
    weights /= np.sqrt(2*np.pi)
    ratios = np.asarray(sample_sizes, dtype=float)/control_size
    centers = normal[:, None]*np.sqrt(ratios)
    halfwidth_scales = np.sqrt(1+ratios)

    def integrand(v):
        halfwidth = threshold*np.sqrt(v/df)*halfwidth_scales
        # Conditional probability of exceeding either boundary for each test.
        outside = np.clip(special.ndtr(centers-halfwidth)
                          + special.ndtr(-centers-halfwidth), 0, 1)
        with np.errstate(divide="ignore"):
            any_outside = -np.expm1(np.log1p(-outside).sum(axis=1))
        return float(weights @ any_outside)*stats.chi2.pdf(v, df)

    probability, error = integrate.quad(integrand, 0, np.inf, epsabs=1e-24,
                                      epsrel=1e-9, limit=250)
    if error > max(1e-22, probability*1e-7):
        raise ArithmeticError("Dunnett tail integration did not converge.")
    return float(np.clip(probability, 0, 1))


def dunnett_vs_control(samples, control):
    """Return signed pooled t, marginal P, and family-adjusted Dunnett P.

    Each sample is one treatment group. Zero within-group variance is allowed
    when the pooled variance is positive; a wholly degenerate family is rejected.
    """
    groups = [np.asarray(control, dtype=float),
              *[np.asarray(sample, dtype=float) for sample in samples]]
    if len(groups) < 2 or any(g.ndim != 1 or len(g) < 2 or not np.isfinite(g).all()
                             for g in groups):
        raise ValueError("Control and every treatment require at least two finite observations.")
    sizes = np.array([len(g) for g in groups])
    df = int(sizes.sum()-len(groups))
    pooled_variance = float(sum(np.sum((g-g.mean())**2) for g in groups)/df)
    if pooled_variance <= 0:
        raise ValueError("Dunnett inference requires a positive pooled within-group variance.")
    differences = np.array([g.mean()-groups[0].mean() for g in groups[1:]])
    standard_errors = np.sqrt(pooled_variance*(1/sizes[0]+1/sizes[1:]))
    t = differences/standard_errors
    raw = 2*stats.t.sf(np.abs(t), df)
    adjusted = np.array([_adjusted_tail(float(abs(value)), tuple(sizes[1:]),
                                      int(sizes[0]), df) for value in t])
    # Independent marginal t tails bound the family union probability.
    if np.any(adjusted < raw*(1-1e-7)) or np.any(adjusted > np.minimum(1, len(t)*raw)*(1+1e-7)):
        raise ArithmeticError("Dunnett probability violates marginal/Bonferroni bounds.")
    return dict(statistic=t, p_raw=raw, p_adjusted=adjusted,
                mean_difference=differences, standard_error=standard_errors,
                df=df, pooled_variance=pooled_variance, family_size=len(t),
                assumptions=ASSUMPTIONS,
                calculation="Deterministic conditional normal/chi-square quadrature; 256 normal nodes")
