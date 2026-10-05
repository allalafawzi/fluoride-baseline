#!/usr/bin/env python3
"""
Statistical helpers shared by the analysis scripts.

The module exists so that wilson() is defined once, with one z value. The same
formula used to sit in five scripts under three different roundings of the
97.5 % normal quantile: 1.96, 1.959963985 and 1.959963984540054. The largest
discrepancy between them is 1e-5 percentage point, and the published intervals
are identical to two decimals under all three, so no reported number changes
here. What changes is that copies of one formula can no longer drift apart.

The value kept is z = 1.959963984540054, which is scipy.stats.norm.ppf(0.975)
in double precision, exact at machine precision; the other two are roundings of
it. It is defined here and imported everywhere else.
"""
import math

# 97.5 % normal quantile in double precision: scipy.stats.norm.ppf(0.975).
# The scripts that used 1.96 or 1.959963985 carried roundings of it.
Z95 = 1.959963984540054


def wilson(k, n, z=Z95):
    """Wilson interval for a proportion.

    Returns (p, lower bound, upper bound). For k = 0 the upper bound is
    z**2 / (n + z**2), which is the form used in section 2.7 for the
    dehalogenase upper bound.

    Wilson is preferred over the normal approximation because it stays inside
    [0, 1] and keeps reasonable coverage for extreme counts, where the normal
    approximation returns negative bounds.
    """
    if n == 0:
        return (0.0, 0.0, 0.0)
    p = k / n
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    demi = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (p, max(0.0, centre - demi), min(1.0, centre + demi))


def clopper_pearson_upper_zero(n, alpha=0.05):
    """Exact Clopper-Pearson upper bound for a zero count.

    Equals 1 - (alpha/2)**(1/n). Given in section 2.7 for comparison with the
    Wilson bound: Wilson is what is reported throughout, and this function only
    shows the gap between the two methods.
    """
    if n == 0:
        return 1.0
    return 1.0 - (alpha / 2.0) ** (1.0 / n)


def design_effect(m_eff, icc):
    """Kish design effect: DEFF = 1 + (m_eff - 1) * ICC.

    Kish (1965), Survey Sampling, section 5.4, equations 5.4.1 and 5.4.2.
    m_eff is the effective cluster size, sum(n**2) / sum(n).

    Kish states the approximations himself: his expressions involve "slight
    approximations", and for clusters of unequal size using the mean size gives
    only "serviceable approximations". Substituting N_eff for N is therefore
    approximate, not exact.
    """
    return 1.0 + (m_eff - 1.0) * icc
