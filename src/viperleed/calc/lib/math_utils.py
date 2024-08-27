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

import warnings

import numpy as np


def angle(*vectors):
    """Return the angle (in radians) between two (sequences of) 2D vectors."""
    nargs = 2
    if len(vectors) != nargs:
        _wrong = 'many' if len(vectors) > nargs else 'few'
        raise TypeError(f'Too {_wrong} arguments. Must be exactly two.')
    # Use cross product for sine, dot product for cosine
    _cross = np.cross(*vectors)
    _dot = np.einsum('...i,...i->...', *vectors)
    return np.arctan2(_cross, _dot)


def cosvec(*vectors):
    """Return the cosine of the angle(s) between two (sequences of) vectors."""
    nargs = 2
    if len(vectors) != nargs:
        _wrong = 'many' if len(vectors) > nargs else 'few'
        raise TypeError(f'Too {_wrong} arguments. Must be exactly two.')
    _norms = np.linalg.norm(vectors, axis=-1)
    _dot = np.einsum('...i,...i->...', *vectors)
    with warnings.catch_warnings():  # Catch 0/0 --> norm == 0
        warnings.filterwarnings('error',
                                message=r'.*true_divide',
                                category=RuntimeWarning)
        try:
            return _dot / np.prod(_norms, axis=0)
        except RuntimeWarning:
            raise ValueError('Some vector has zero norm') from None


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


def lcm(*pair_of_integers):
    """Return the lowest common multiple of two (sequences of) integers."""
    nargs = 2
    if len(pair_of_integers) != nargs:
        _wrong = 'many' if len(pair_of_integers) > nargs else 'few'
        raise TypeError(f'Too {_wrong} arguments. Must be exactly two.')
    return np.prod(pair_of_integers, axis=0) // np.gcd(*pair_of_integers)
