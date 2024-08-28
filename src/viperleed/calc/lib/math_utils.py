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

from functools import wraps
import warnings

import numpy as np


def _with_two_args(func):
    """Raise TypeError if `func` is not called with exactly 2 arguments."""
    nargs = 2
    @wraps(func)
    def _wrapper(*args):
        if len(args) != nargs:
            _wrong = 'many' if len(args) > nargs else 'few'
            raise TypeError(f'{func.__name__}: Too {_wrong} '
                            'arguments. Must be exactly two.')
        vec_1, vec_2 = (np.asanyarray(a) for a in args)
        if vec_1.shape != vec_2.shape:
            raise ValueError(f'{func.__name__}: Inconsistent shapes '
                                f'{vec_1.shape} != {vec_2.shape}')
        return func(vec_1, vec_2)
    return _wrapper


@_with_two_args
def angle(*vectors):
    """Return the angle (in radians) between two (sequences of) 2D vectors."""
    vec_1, vec_2 = vectors
    if vec_1.shape and vec_1.shape[-1] != 2:
        raise ValueError('angle: only accepts 2D vectors. '
                         f'Found {vec_1.shape[-1]} ones')
    # Use cross product for sine, dot product for cosine
    # The cross product of 2D arrays via np.cross has been deprecated
    # in NumPy 2.0. Use the one-line solution from NumPy #26620:
    _cross = vec_1[..., 0] * vec_2[..., 1] - vec_1[..., 1] * vec_2[..., 0]
    _dot = np.einsum('...i,...i->...', vec_1, vec_2)
    return np.arctan2(_cross, _dot)


@_with_two_args
def cosvec(*vectors):
    """Return the cosine of the angle(s) between two (sequences of) vectors."""
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


@_with_two_args
def lcm(*pair_of_integers):
    """Return the lowest common multiple of two (sequences of) integers."""
    return np.prod(pair_of_integers, axis=0) // np.gcd(*pair_of_integers)
