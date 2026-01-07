"""viperleed.calc.lib.rfactor.groups"""

__authors__ = ('Alexander M. Imre (@amimre)',)
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2026-01-07'
__license__ = 'GPLv3+'

from . import xp


def group_rfactors(numerators, denominators, groups=None, num_groups=None):
    """Calculate grouped R-factor from numerators and denominators.

    Parameters
    ----------
    numerators : xp.ndarray
        Array of R-factor numerators of shape (n_beams,)
    denominators : xp.ndarray
        Array of R-factor denominators of shape (n_beams,)
    groups : None | "beam" | array-like of int, optional
        - None: return overall R-factor = sum(num) / sum(den)
        - "beam": return per-beam R-factors = num / den
        - array of ints: group id per beam; returns per-group R-factors.
          Group ids are assumed to be non-negative integers.
    num_groups : int, optional
        Number of groups, required if groups is array-like of int.
        Must be specified to ensure compatibility with just-in-time
        compilation and static array shapes.

    Returns
    -------
    xp.ndarray
        R-factors per beam of shape (n_beams,) if groups is "beam",
        shape (1,) if groups is None, or shape (n_groups,) if
        groups is an array-like of integers.
    """
    # check numerators and denominators have the same shape
    if numerators.shape != denominators.shape:
        raise ValueError(
            'Numerators and denominators must have the same shape.'
        )

    if groups is None:
        return xp.sum(numerators) / xp.sum(denominators)
    if isinstance(groups, str) and groups == 'beam':
        return numerators / denominators

    # otherwise, assume groups is array-like of integers
    groups = xp.asarray(groups, dtype=int)

    num_groups = xp.max(groups) + 1
    grouped_numerators = xp.bincount(
        groups, weights=numerators, minlength=num_groups
    )
    grouped_denominators = xp.bincount(
        groups, weights=denominators, minlength=num_groups
    )
    return grouped_numerators / grouped_denominators
