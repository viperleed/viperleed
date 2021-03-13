"""
======================================
  ViPErLEED Graphical User Interface
======================================
      *** module guilib.helpers ***

Created: 2021-03-11
Author: Michele Riva

This module collects helper functions used in various places in the GUI
"""


import numpy as np


def two_by_n_array_to_tuples(two_by_n, axis=None):
    """
    Converts a 2xN or Nx2 array into a zip object, that will return
    tuples when iterated over

    Parameters
    ----------
    two_by_n : numpy.ndarray
        Can have shape (N, 2) or (2, N). The output will always be an iterator
        of tuples, independent of the shape of the input
    axis : 0, 1, or None
        Which axis contains the 2-tuples. Use axis=0 if two_by_n[i] is the
        i-th tuple, axis=1 if two_by_n[:, i] is the i-th tuple. If None or
        omitted, the axis is inferred from which of the dimensions is equal
        to two, unless both are. In that case axis is mandatory.

    Returns
    -------
    zip-object
    """
    # Make sure we get an array. It would be possible to also do this with
    # other array-like objects but it's not implemented
    if not isinstance(two_by_n, np.ndarray):
        raise TypeError("Need a numpy array as input")
    if (all(length != 2 for length in two_by_n.shape)
        or len(two_by_n.shape) != 2):
        raise ValueError("Invalid shape of index array. Expected (2, N)"
                             f" or (N, 2), found {two_by_n.shape}")
    if axis is None:
        # determine which one is the right axis
        if two_by_n.shape[0] == 2 and two_by_n.shape[1] != 2:
            axis = 0
        elif two_by_n.shape[0] != 2 and two_by_n.shape[1] == 2:
            axis = 1
        else:
            raise RuntimeError("Cannot determine automatically which is the"
                               "correct axis, as two_by_n.shape is (2, 2). "
                               "Provide the correct axis with the axis "
                               "optional argument")
    if axis == 0:  # Shape == (2, N)
        return zip(two_by_n[0], two_by_n[1])
    # Shape == (N, 2)
    return zip(two_by_n[:, 0], two_by_n[:, 1])


def two_by_two_array_to_tuple(two_by_two):
    """
    Convenience function that converts a 2x2 array to a 2x2 tuple. Does no
    type or shape checking!
    """
    return tuple(map(tuple, two_by_two))


def conventional_angles(theta, phi):
    """
    Given two angles theta and phi in degrees, returns a modified version such
    that 0 <= theta <= 180 and 0 <= phi < 360, and such that they represent
    the same vector as before

    Parameters
    ----------
    theta : float
        polar angle in degrees, negative values affect the value of phi
    phi : float
        azimuthal angle in degrees, returned modulo 360

    Returns
    -------
    (theta, phi) : (float, float)
    """
    phi %= 360
    if theta < 0:
        # need to change the sign of theta, and adapt phi accordingly
        theta *= -1
        phi -= 180
        phi %= 360
    return theta, phi


def remove_duplicates(data, return_type=None):
    """
    Removes duplicates from the input. The function preserves order, but works
    with generator objects only if all the items are hashable.

    Parameters
    ----------
    data : iterable
    return_type : callable, optional (default: same as input or tuple)
        A callable that can instantiate an iterable, which is then returned.
        If not given, the same type as the input is returned. It is generally
        safer to pass return_type, as only a few cases are handled as expected.
        Falls back on returning a tuple in case return_type would raise abs
        TypeError. Does not work as expected for multi-dimensional numpy arrays!

    Returns
    -------
    uniques : return_type if given, type(data) if it doesn't raise errors,
              tuple otherwise
    """
    # try to use an order-preserving O(1) algorithm
    try:
        ret = dict.fromkeys(data)
    except TypeError:
        # fall back to an O(n^2) algorithm in case there are unhashable elements
        ret = [d for i, d in enumerate(data) if d not in data[:i]]

    if return_type is None:
        if isinstance(data, np.ndarray):
            return_type = np.array
        else:
            return_type = type(data)

    if return_type is str:
        return_type = lambda x: "".join(map(str, x))
    elif return_type in (np.array, np.asarray):
        return_type = lambda x: np.fromiter(x, data.dtype)

    try:
        return_type(ret)
    except (TypeError, ValueError):
        return_type = tuple
    return return_type(ret)