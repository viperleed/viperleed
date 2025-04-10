"""Module matrix of viperleed.calc.lib.

Collects exceptions, functions for matrix operations, and various
transformation matrices.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-08-27'
__license__ = 'GPLv3+'

import numpy as np


class NonIntegerMatrixError(ValueError):
    """A matrix that should have integer values does not."""


class SingularMatrixError(ValueError, ZeroDivisionError):
    """A matrix that needs inversion is singular."""


def ensure_integer_matrix(matrix, eps=1e-6):
    """Return a rounded version of `matrix`. Raise if it isn't integer."""
    rounded = np.round(matrix)
    if not np.isfinite(rounded).all():
        raise NonIntegerMatrixError('Matrix contains non-finite '
                                    f'values: {matrix}')
    if not np.allclose(matrix, rounded, atol=eps, rtol=0):
        raise NonIntegerMatrixError('Matrix contains non-integer values '
                                    f'within tolerance of {eps}: {matrix}')
    return rounded


def rotation_matrix(angle, dim=2):
    """Return a matrix for in-plane rotation by `angle`.

    Parameters
    ----------
    angle : float
        The rotation angle in radians. Positive angles correspond to
        counterclockwise rotations when the matrix returned operates
        from the left on column vectors.
    dim : int, optional
        The size of the matrix to return. Default is 2.

    Returns
    -------
    numpy.ndarray
        Shape (`dim`, `dim`). The rotation part is in the top-left 2x2
        corner.
    """
    if dim < 2:  # pylint: disable=magic-value-comparison  # Clear
        raise ValueError('Rotation matrix needs at least dimension 2')
    _cos, _sin = np.cos(angle), np.sin(angle)
    rotation = (_cos, -_sin), (_sin, _cos)
    matrix = np.eye(dim, dtype=float)
    matrix[:2, :2] = rotation
    return matrix


def rotation_matrix_order(order, dim=2):
    """Return a matrix for an `order`-fold in-plane rotation.

    Parameters
    ----------
    order : float
        The rotation order, typically an integer. Positive orders
        correspond to counterclockwise rotations when the matrix
        returned operates from the left on column vectors.
    dim : int, optional
        The size of the matrix to return. Default is 2.

    Returns
    -------
    numpy.ndarray
        Shape (`dim`, `dim`). The rotation part is in the top-left 2x2
        corner.
    """
    angle = 2 * np.pi / order
    return rotation_matrix(angle, dim=dim)
