"""Module coordinates of viperleed.calc.lib.

Defines functions for handling atomic coordinates."""


__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-08-25'
__license__ = 'GPLv3+'

import itertools

import numpy as np


COLLAPSE_EPS = 1e-8  # Default for collapsing fractional coordinates


def add_edges_and_corners(cartesian, fractional, releps, ucell, props=None):
    """Add atoms at edges and corners to the Cartesian coordinates given.

    Notice that edges and corners will be added only for dimensions
    up to the smallest among: shape(fractional)[1], len(releps),
    len(ucell). The only requirements are that `fractional` contains
    at most the same number of elements as `cartesian`, and that the
    second axis of `cartesian` and `ucell` are broadcastable. This
    allows, e.g., to pass in 3D coordinates but only add repeats
    for 2D edges and corners.

    Parameters
    ----------
    cartesian : Sequence
        The Cartesian coordinates to which replicas of atoms at
        edges and corners will be appended. Shape (N, 2) or (N, 3).
    fractional : Sequence
        Fractional coordinates corresponding to `cartesian`. Should
        have shape (M, m), with M <= N.
    releps : Sequence
        Tolerance on fractional coordinates to consider an atom from
        `cartesian` to be close to an edge. Shape (2,) or (3,).
    ucell : Sequence
        Unit cell to be used for the displacement of atoms close to
        edges/corners. Shape (2, 2) or (3, 3), and such that unit
        cell vectors are rows. The second dimension should match
        the second one of `cartesian`.
    props : list or None, optional
        If given, it should have the same length as `cartesian`. It
        is interpreted as a list of properties that should also be
        duplicated whenever an atom is duplicated. Default is None.

    Returns
    -------
    cartesian_extended : numpy.ndarray
        Cartesian atomic coordinates to which extra atoms have been
        appended at the end for all those atoms in `cartesian` that
        were at edges/corners.
    props_extended : list or None
        Corresponding properties of the added atoms. None if `props`
        was not given to begin with.
    """
    for i, pos in enumerate(fractional):
        addlist = []
        for posj, epsj, uvec in zip(pos, releps, ucell):
            if abs(posj) < epsj:        # 'left' edge
                addlist.append(cartesian[i] + uvec)
            elif abs(posj - 1) < epsj:  # 'right' edge
                addlist.append(cartesian[i] - uvec)
        # pylint: disable=magic-value-comparison  # For 2 & 3
        if len(addlist) == 2:
            # 2D corner - add the diagonally opposed one
            addlist.append(sum(addlist) - cartesian[i])
        elif len(addlist) == 3:
            # 3D corner - add all diagonally opposed points
            addlist.extend(
                p1 + p2 - cartesian[i]
                for p1, p2 in itertools.combinations(addlist, 2)
                )
            addlist.append(sum(addlist[:3]) - 2 * cartesian[i])

        if not addlist:
            continue
        cartesian = np.concatenate((cartesian, addlist))
        if props is not None:
            props.extend([props[i]]*len(addlist))
    return cartesian, props


def collapse(cartesians, ucell, ucell_inv=None, method='floor'):
    """Collapse Cartesian coordinates to the base cell.

    Parameters
    ----------
    cartesians : numpy.ndarray
        Cartesian coordinates to be collapsed. Shape (N, 2) or (N, 3).
    ucell : numpy.ndarray
        Basis vectors (rows) of the unit cell to be used for collapsing.
    ucell_inv : numpy.ndarray or None, optional
        The inverse of ucell. If not given or None, it is computed
        on the fly. Default is None.
    method : {'floor', 'round'}
        Which method should be used for collapsing. 'floor'
        returns Cartesian coordinates collapsed to the (0, 0)
        cell, i.e., with all fractional coordinates in range
        [0, 1]. 'round' collapses the coordinates so that they
        are as close to zero as possible. Fractional coordinates
        are then in range [-0.5, 0.5]. The default is 'floor'.

    Returns
    -------
    collapsed_cartesians : numpy.ndarray
        The Cartesian coordinates, collapsed according to method
    fractional_coordinates : numpy.ndarray
        The corresponding fractional coordinates
    """
    if ucell_inv is None:
        ucell_inv = np.linalg.inv(ucell)
    fractional = collapse_fractional(cartesians.dot(ucell_inv),
                                     method=method)
    return fractional.dot(ucell), fractional


def _floor_eps(eps):
    """Return the floored-int version of values after adding eps."""
    def _floor(values):
        return np.floor(values + eps)
    return _floor


def collapse_fractional(coordinates, method='floor',
                        eps=COLLAPSE_EPS, in_place=False):
    """Collapse fractional coordinates to the base cell.

    Parameters
    ----------
    coordinates : numpy.ndarray
        Fractional coordinates to be collapsed. Shape (N, 2) or (N, 3).
    method : {'floor', 'round'}
        Which method should be used for collapsing.  'round' collapses
        the coordinates so that they are as close to zero as possible.
        Fractional coordinates are then in range [-0.5, 0.5]. 'floor'
        returns fractional coordinates collapsed to the (0, 0) cell,
        that is, with all fractional coordinates in range [-`eps`,
        1-`eps`]. This ensures that collapsing a fractional coordinate
        of (essentially) 1 gives (essentially) zero. Default is 'floor'.
    eps : float or Sequence, optional
        Used only if method == 'floor'. The (fractional) tolerance
        for collapsing. If a sequence, it should have as many items
        as the second axis of coordinates. Default is 1e-8.
    in_place : bool, optional
        If True, the function modifies directly `coordinates`. Can be
        used to save some memory. Default is False.

    Returns
    -------
    collapsed : numpy.ndarray
        The collapsed fractional coordinates. Will be the same object
        passed in if in_place is True-thy.

    Raises
    ------
    ValueError
        If method is not one of the valid methods.
    """
    _methods = {'f': _floor_eps(eps), 'r': np.round}
    try:
        round_ = _methods[method[0]]
    except (TypeError, IndexError, KeyError) as exc:
        raise ValueError('collapse_fractional: Unknown '
                         f'method={method}') from exc
    if in_place:
        collapsed = coordinates
        collapsed -= round_(coordinates)
    else:
        collapsed = coordinates - round_(coordinates)
    return collapsed
