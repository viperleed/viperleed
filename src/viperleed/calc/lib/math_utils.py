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


def angle(v1, v2):
    """Return the angle between two 2D vectors."""
    # Use cross product for sine, dot product for cosine
    return np.arctan2(v1[0]*v2[1] - v1[1]*v2[0], v1[0]*v2[0] + v1[1]*v2[1])


def cosvec(x, y):
    """Return the cosine of the angle between two vectors.

    Parameters
    ----------
    x : Sequence
        First vector
    y : Sequence
        Second vector

    Returns
    -------
    float
        Cosine of the angle between the two vectors.
    """
    return np.dot(x, y) / (np.linalg.norm(x) * np.linalg.norm(y))


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


def lcm(a, b):
    """Return the lowest common multiple of two integers."""
    return a * b // np.gcd(a, b)
