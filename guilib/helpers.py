"""Module viperleed.guilib.helpers.

======================================
  ViPErLEED Graphical User Interface
======================================
      *** module guilib.helpers ***

Created: 2021-03-11
Author: Michele Riva

This module collects helper functions used in various places in the GUI
"""

import itertools

import numpy as np


def two_by_n_array_to_tuples(two_by_n, axis=None):
    """Convert a 2xN or Nx2 array into a zip object.

    Iterating over the return value will return tuples.

    Parameters
    ----------
    two_by_n : numpy.ndarray
        Can have shape (N, 2) or (2, N). The output will always be
        an iterator of tuples, independent of the shape of the input
    axis : 0, 1, or None
        Which axis contains the 2-tuples. Use axis=0 if two_by_n[i]
        is the i-th tuple, axis=1 if two_by_n[:, i] is the i-th tuple.
        If None or omitted, the axis is inferred from which of the
        dimensions is equal to two, unless both are. In that case
        axis is mandatory.

    Returns
    -------
    zip-object

    Raises
    ------
    TypeError
        When the input is not a numpy.ndarray
    ValueError
        When the shape of the input is neither (2, N) nor (N, 2)
    RuntimeError
        When shape = (2, 2) and no axis was given
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
    """Convert a 2x2 array to a 2x2 tuple.

    Does no type or shape checking!

    Parameters
    ----------
    two_by_two : numpy.ndarray

    Returns
    -------
    tuple of tuples
    """
    return tuple(map(tuple, two_by_two))


def two_d_iterable_to_array(iterable, dtype=float, shape=None):
    """Return a numpy.ndarray from a 2D iterable.

    This function uses a faster alternative than simply using
    np.asarray(iterable, dtype=dtype). Rough speed tests
    suggest that this is a factor of 1.5-5.5 faster than
    np.asarray when running on a list of tuples, and a
    factor of 3-35 faster when running on a list of
    gl.BeamIndex (likely similar performance on other
    subclasses of tuple).

    No explicit check is done on shapes.

    Parameters
    ----------
    iterable : sequence
        len(shape) == 2
    dtype : numpy.dtype, default=float
        Data type that will be used for the output array.
    shape : tuple, default=None
        If given, the return array is reshaped to shape.

    Returns
    -------
    numpy.ndarray
    """
    if isinstance(iterable, np.ndarray):
        return iterable
    arr = np.fromiter(itertools.chain.from_iterable(iterable),
                      dtype=dtype)
    if shape is None:
        return arr
    return arr.reshape(*shape)


def conventional_angles(theta, phi):
    """Return angles in a well defined convention.

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
    """Remove duplicates from the input.

    The function preserves order, but works with generator
    objects only if all the items are hashable.

    Parameters
    ----------
    data : iterable
    return_type : callable, optional (default: same as input or tuple)
        A callable that can instantiate an iterable, which
        is then returned.  If not given, the same type as
        the input is returned.  It is generally safer to
        pass return_type, as only a few cases are handled
        as expected.  Falls back on returning a tuple in
        case return_type would raise TypeError. Does not
        work as expected for multi-dimensional numpy arrays!

    Returns
    -------
    uniques : return_type if given, type(data) if it doesn't
              raise errors, tuple otherwise
    """
    def __return_string(elements):
        return "".join(map(str, elements))

    def __return_array(elements):
        return np.fromiter(elements, data.dtype)

    # try to use an order-preserving O(1) algorithm
    try:
        ret = dict.fromkeys(data)
    except TypeError:
        # fall back to an O(n^2) algorithm in case there are
        # unhashable elements
        ret = [d for i, d in enumerate(data) if d not in data[:i]]

    if return_type is None:
        if isinstance(data, np.ndarray):
            return_type = np.array
        else:
            return_type = type(data)

    if return_type is str:
        return_type = __return_string
    elif return_type in (np.array, np.asarray):
        return_type = __return_array

    try:
        ret = return_type(ret)
    except (TypeError, ValueError):
        ret = tuple(ret)
    return ret


def single_spaces_only(input_string):
    """Return a string with double-spaces converted to single ones.

    Paramters
    ---------
    input_string : str
        The string to be processed

    Returns
    -------
    str
    """
    while "  " in input_string:
        input_string = input_string.replace("  ", " ")
    return input_string


def array2string(matrix):
    """Return a 1-line string representation of an array."""
    matrix = np.array2string(matrix, separator=',', suppress_small=True)
    return single_spaces_only(matrix).replace('\n', '')


def prime_numbers():
    """Yield an infinite number of prime numbers.

    Yields
    ------
    int
        The next prime number.
    """
    # Algorithm is taken from https://stackoverflow.com/a/10733621/849891
    # and is essentially Erastothenes sieve.
    yield from (2, 3, 5, 7)
    sieve = {}
    primes = prime_numbers()
    prime = next(primes) and next(primes)
    assert prime == 3
    prime_squared = prime*prime
    for i in itertools.count(9, 2):
        if i in sieve:            # composite
            step = sieve.pop(i)
        elif i < prime_squared:   # prime
            yield i
            continue
        else:                     # composite, == prime_squared
            assert i == prime_squared
            step = 2*prime
            prime = next(primes)
            prime_squared = prime*prime
        i += step
        while i in sieve:
            i += step
        sieve[i] = step


def equal_dicts(dict_a, dict_b, ignore_keys=[]):
    """Return dict_a == dict_b, optionally ignoring some keys."""
    if not ignore_keys:
        return dict_a == dict_b
    a_keys = set(dict_a).difference(ignore_keys)
    b_keys = set(dict_b).difference(ignore_keys)
    return a_keys == b_keys and all(dict_a[k] == dict_b[k] for k in a_keys)