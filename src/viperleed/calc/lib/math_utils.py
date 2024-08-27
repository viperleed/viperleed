"""Module math_utils of viperleed.calc.lib.

Defines basic mathematical functions.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-08-26'
__license__ = 'GPLv3+'

import re

import numpy as np


def floor_eps(eps):
    """Return a callable producing floored values after adding `eps`.

    Parameters
    ----------
    eps : float or numpy.ndarray
        The quantity to add to the values before flooring. If an array,
        the values to be floored must have the same shape, or be a
        single float value.

    Returns
    -------
    _floor : callable
        A function returning the floored version of its argument, after
        adding `eps`. The function accepts either a single value or
        a numpy array.
    """
    def _floor(values):
        return np.floor(values + eps)
    return _floor
