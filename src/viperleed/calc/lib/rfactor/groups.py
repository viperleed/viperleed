"""viperleed.calc.lib.rfactor.groups."""

__authors__ = ('Alexandra M. Imre (@alexmiame)',
               'Florian Kraushofer (@fkraushofer)',)
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2026-01-07'
__license__ = 'GPLv3+'

from viperleed.calc.lib import dynamic_numerical_lib as dnl


def group_rfactors(numerators, denominators, groups=None, num_groups=None):
    r"""Calculate grouped R-factor from numerators and denominators.

    R-factors are often calculated over multiple beams, where the
    R-factor for a set $`g`$ of beams is given by
    $`R_g = \sum_{b \in g} \int f(I_{1,b}, I_{2,b}) \, dE \; / \;
    \sum_{b \in g} \int g(I_{1,b}, I_{2,b}) \, dE`$,
    where $`I_{1,b}`$ and $`I_{2,b}`$ are the intensities of beam $`b`$
    for the two data sets being compared and $`f`$ and $`g`$ are
    functions specific to the R-factor being calculated.
    Here, `numerators` and `denominators` are arrays of shape
    (n_beams,) containing the values of $`\int f(I_{1,b}, I_{2,b})
    \, dE`$ and $`\int g(I_{1,b}, I_{2,b}) \, dE`$ for each beam
    $`b`$.
    This function allows grouping these beams into arbitrary sets
    and calculating the R-factor for each set.
    The grouping can be specified in three ways:
    - If `groups` is None, the overall R-factor is calculated.
    - If `groups` is "beam", individual per-beam R-factors are
    calculated.
    - If `groups` is an array-like of integers of shape (n_beams,),
    the beams are grouped according to the group ids specified in
    `groups`, and per-group R-factors are calculated.

    Parameters
    ----------
    numerators : dnl.xp.ndarray
        Array of R-factor numerators of shape (n_beams,)
    denominators : dnl.xp.ndarray
        Array of R-factor denominators of shape (n_beams,)
    groups : None or "beam" or array-like of int, optional
        - None: return overall R-factor = sum(num) / sum(den)
        - "beam": return per-beam R-factors = num / den
        - array of ints: group id per beam; returns per-group R-factors.
            Group ids are assumed to be non-negative integers.
    num_groups : int, optional
        Number of groups, required if groups is array-like of int.
        If not specified, it will be inferred from the
        `groups`. Must be specified for JIT compiled code to ensure
        static array sizes.

    Returns
    -------
    dnl.xp.ndarray
        R-factors per beam of shape (n_beams,) if groups is "beam",
        shape (1,) if groups is None, or shape (n_groups,) if
        groups is an array-like of integers.
        Any invalid groups (e.g. groups with zero denominator) will
        return have R-factor set to nan.
    """
    # check numerators and denominators have the same shape
    if numerators.shape != denominators.shape:
        raise ValueError(
            'Numerators and denominators must have the same shape.'
        )

    if groups is None:
        return dnl.xp.sum(numerators) / dnl.xp.sum(denominators)
    if isinstance(groups, str) and groups == 'beam':
        return numerators / denominators

    # otherwise, assume groups is array-like of integers
    _groups = dnl.xp.asarray(groups, dtype=int)

    if num_groups is None:
        num_groups = dnl.xp.max(groups) + 1
    grouped_numerators = dnl.bincount(
        _groups, weights=numerators, length=num_groups
    )
    grouped_denominators = dnl.bincount(
        _groups, weights=denominators, length=num_groups
    )

    mask = grouped_denominators != 0
    safe_denominators = dnl.xp.where(mask, grouped_denominators, 1.0)
    result = dnl.xp.where(mask, grouped_numerators / safe_denominators, dnl.xp.nan)
    return result
