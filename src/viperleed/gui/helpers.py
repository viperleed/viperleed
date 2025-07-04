"""Module helpers of viperleed.gui.

This module collects helper functions used in various spots in the GUI.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2021-03-11'
__license__ = 'GPLv3+'

import ast
import itertools
from pathlib import Path
import re
import sys

import numpy as np

import viperleed


def integer_part_length(*numbers):
    return max(len(f"{int(number)}") for number in numbers)


def format_floats(format_specs, *numbers):
    """
    Formats floats to have them all aligned to the point, and according to
    format_spec. format_spec is in the form

    [[fill]align][sign][#][0][minimumwidth][.precision][type]

    - fill and align will be used to pad the overall resulting string
    - sign, #, 0 are discarded
    - minimumwidth is the minimum number of characters used to represent the
      integer parts. Defaults to the minimum number necessary to have all
      numbers aligned to the decimal point
    - precision is the number of decimal digits. Defaults to 5.
    """
    int_len = integer_part_length(*numbers)

    # standard pattern for format_specs
    pattern = (r"^(?P<align>[<>=^]\d+)?"
               r"[\+\- ]?"
               r"[#]?"
               r"[0]?"
               r"(?P<minwidth>\d*)"
               r"([.](?P<prec>\d+))?"
               r"[bdcoxXneEfFgG\%]?$")
    m = re.match(pattern, format_specs)
    if m is None:
        raise TypeError("Incorrect format specifier for type f.")
    # process the format specifications:
    # (1) Number of decimal digits
    prec = int(m['prec'] or '5')
    # (2) Minimum width of integer part
    #    (will left-pad with spaces if needed)
    min_w = int(m['minwidth'] or 0) or int_len
    tot_w = max(min_w + prec, int_len + prec) + 1  # +1 for decimal point

    format = f">{tot_w}.{prec}f"
    raw = ','.join(f"{float(number):{format}}" for number in numbers)
    align = ''
    if m.group('align'):
        align = f"{m.group('align')}"

    return f"{raw:{align}}"


def resources_path(dir_name):
    """Return the correct path to `dir_name`.

    This is useful when building an executable from pyinstaller.
    When built from pyinstaller, it takes the path relative to the
    temporary path in which the "exe" is extracted during execution.

    Parameters
    ----------
    dir_name : str
        Relative path within the viperleed source tree.

    Returns
    -------
    str
        Absolute version of `dir_name` (in the tree where viperleed
        is installed).
    """
    try:
        root = Path(sys._MEIPASS).resolve()  # pyinstaller
    except AttributeError:
        # NB: __file__ points to __init__.py
        root = Path(viperleed.__file__).resolve().parent
    return str(root/dir_name)


def string_matrix_to_numpy(str_matrix, dtype=float, needs_shape=tuple()):
    """
    Takes a string representing a matrix with float values and parses it into a
    numpy array of the right shape. It can also check if the string represents
    a matrix of a specific shape.

    Parameters
    ----------
    - str_matrix: str
                  String to be parsed to extract a matrix. It needs to be
                  formatted as an array-like input acceptable for numpy
    - dtype: data type; default float
             Data type that the numpy array returned will have
    - needs_shape: array-like; default = tuple()
                   if this keyword argument is given, the function checks that
                   the shape of the parsed matrix matches the input

    Returns
    -------
    np.array or None

    If needs_shape is given, the function returns None in case the shapes do not
    match. Otherwise always returns a numpy array.
    """
    if not isinstance(str_matrix, str):
        raise TypeError("string_matrix_to_numpy: str_matrix should be string. "
                        f"Found {type(str_matrix).__name__!r} instead")
    try:
        _ = len(needs_shape)
    except TypeError:
        raise TypeError("string_matrix_to_numpy: needs_shape should be an "
                        f"array-like. Found {type(needs_shape).__name__!r} "
                        "instead") from None
    try:
        matrix = np.asarray(ast.literal_eval(str_matrix), dtype=float)
    except SyntaxError as err:
        raise RuntimeError("Could not extract a matrix from string "
                           f"{str_matrix}. Most likely a bracket or "
                           "a comma is missing") from err

    # if dtype=int, round first, then typecast
    if dtype == int:
        matrix = matrix.round().astype(int)

    if len(needs_shape) > 0 and matrix.shape != tuple(needs_shape):
        raise RuntimeError("string_matrix_to_numpy: matrix has the wrong "
                           f"shape. Expected {tuple(needs_shape)}, found "
                           f"{matrix.shape}")
    return matrix


# Disable due to pylint bug for sequence of allowed values
#  ('axis' parameter). Bug is still there as of pylint 2.9.3
# pylint: disable=missing-param-doc,missing-type-doc
def two_by_n_array_to_tuples(two_by_n, axis=None):
    """Convert a 2xN or Nx2 array into a zip object.

    Iterating over the return value will return tuples.

    Parameters
    ----------
    two_by_n : numpy.ndarray
        Can have shape (N, 2) or (2, N). The output will always be
        an iterator of tuples, independent of the shape of the input
    axis : {None, 0, 1}
        Which axis contains the 2-tuples. Use axis=0 if two_by_n[i]
        is the i-th tuple, axis=1 if two_by_n[:, i] is the i-th tuple.
        If None or omitted, the axis is inferred from which of the
        dimensions is equal to two, unless both are. In that case
        axis is mandatory.

    Returns
    -------
    zip
        Iterator of tuples containing the elements in the
        original array

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

    # Disable, as I find it clearer with the explicit check
    # pylint: disable=compare-to-zero
    if axis == 0:  # Shape == (2, N)
        return zip(two_by_n[0], two_by_n[1])
    # Shape == (N, 2)
    return zip(two_by_n[:, 0], two_by_n[:, 1])
# pylint: enable=missing-param-doc,missing-type-doc


# Disable for performance reasons
# pylint: disable=bad-builtin
def two_by_two_array_to_tuple(two_by_two):
    """Convert a 2x2 array to a 2x2 tuple.

    Does no type or shape checking!

    Parameters
    ----------
    two_by_two : numpy.ndarray

    Returns
    -------
    tuple of tuples
        Converted array
    """
    return tuple(map(tuple, two_by_two))
# pylint: enable=bad-builtin


def two_d_iterable_to_array(iterable, dtype=float, shape=None):
    """Return a numpy.ndarray from a 2D iterable.

    This function uses a faster alternative than simply using
    np.asarray(iterable, dtype=dtype). Rough speed tests suggest
    that this is a factor of 1.5-5.5 faster than np.asarray when
    running on a list of tuples, and a factor of 3-35 faster when
    running on a list of BeamIndex (likely similar performance
    on other subclasses of tuple).

    Parameters
    ----------
    iterable : Sequence
        The sequence to be converted to numpy.ndarray
    dtype : numpy.dtype, optional
        Data type that will be used for the output array.
        Default is float.
    shape : tuple or None, optional
        If given and not None, the return array is reshaped
        to shape. Works correctly only if len(shape) == 2,
        but this is unchecked for speed reasons. Default is None

    Returns
    -------
    numpy.ndarray
        Array version of the iterable passed
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

    Given two angles theta and phi in degrees, returns a modified
    version such that 0 <= theta <= 180 and 0 <= phi < 360, and
    such that they represent the same vector as before

    Parameters
    ----------
    theta : float
        polar angle in degrees, negative values affect the value of phi
    phi : float
        azimuthal angle in degrees, returned modulo 360

    Returns
    -------
    theta : float
        Polar angle in conventional setting, i.e., 0 <= theta <= 180
    phi : float
        Azimuthal angle in conventional setting, i.e., 0 <= phi < 360
    """
    phi %= 360
    if theta < 0:
        # need to change the sign of theta, and adapt phi accordingly
        theta *= -1
        phi -= 180
        phi %= 360
    return theta, phi


# Disable as it does not really make sense to write separate
# functions for different data types. Perhaps one could make
# it a class.
# pylint: disable=too-complex
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
    uniques : object
        Type is return_type if given, type(data) if
        it doesn't raise errors, tuple otherwise
    """
    def __return_string(elements):
        return "".join(str(el) for el in elements)

    def __return_array(elements):
        return np.fromiter(elements, data.dtype)

    # try to use an order-preserving O(1) algorithm
    try:
        ret = dict.fromkeys(data)
    except TypeError as exc:
        if 'unhashable' not in exc.args[0]:
            raise
        # fall back to an O(n^2) algorithm in case there are
        # non-hashable elements
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
        ret_with_type = return_type(ret)
    except (TypeError, ValueError):
        ret_with_type = tuple(ret)
    return ret_with_type
# pylint: enable=too-complex


def single_spaces_only(input_string):
    """Return a string with double-spaces converted to single ones.

    Parameters
    ---------
    input_string : str
        The string to be processed

    Returns
    -------
    str
        A version of input_string that
        contains only single white spaces
    """
    while "  " in input_string:
        input_string = input_string.replace("  ", " ")
    return input_string


def array_to_string(matrix):
    """Return a 1-line string representation of an array."""
    matrix = np.array2string(matrix, separator=',', suppress_small=True)
    matrix = single_spaces_only(matrix).replace('\n', '')
    return matrix.replace('[ ','[').replace(' ]', ']').replace(' ,', ',')


# Disable, as pylint does not detect
# that this is an infinite generator
# pylint: disable=stop-iteration-return
def prime_numbers():
    """Yield an infinite number of prime numbers.

    Yields
    ------
    int
        The next prime number.
    """
    # Algorithm from https://stackoverflow.com/a/10733621/849891
    # and it is essentially Erastothene's sieve.
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
# pylint: enable=stop-iteration-return


def equal_dicts(dict_a, dict_b, ignore_keys=None):
    """Return dict_a == dict_b, optionally ignoring some keys.

    Parameters
    ----------
    dict_a : dict
        The first of the dictionaries to be compared
    dict_b : dict
        The second of the dictionaries to be compared
    ignore_keys : Sequence or None, optional
        Keys to be skipped while comparing the dictionaries.
        If not given or None, all keys are compared. Default
        is None.

    Returns
    -------
    bool
        True if dict_a == dict_b
    """
    if not ignore_keys:
        return dict_a == dict_b
    a_keys = set(dict_a).difference(ignore_keys)
    b_keys = set(dict_b).difference(ignore_keys)
    return a_keys == b_keys and all(dict_a[k] == dict_b[k] for k in a_keys)


def is_integer_matrix(matrix, eps=1e-3):
    """Return whether a matrix contains only integers.

    Will return True also if all elements are floats that
    are close to integer values.

    Parameters
    ----------
    matrix : Sequence
        Matrix to be tested
    eps : float, optional
        Absolute tolerance to determine if a float is
        close to an integer.  Default is 1e-3.

    Returns
    -------
    commensurate : bool
        True if matrix is commensurate
    """
    if matrix is None:
        return False

    if isinstance(matrix, np.ndarray) and matrix.dtype == int:
        return True

    matrix = np.asarray(matrix)

    return bool(np.all(np.abs(matrix - matrix.round()) < eps))


def justify_sequences(sequences, filler=None):
    """Justify sequences to same length.

    Each of the given sequences is left-justified into a sequence
    as long as the longest sequence given, with missing vaules
    replaced by filler.

    Parameters
    ----------
    sequences : Sequence
        Sequences to be justified. These sequences are not altered.
    filler : object, optional
        Inserted object in vacant spaces. Default is None.

    Returns
    -------
    justified : Sequence
        Sequence with justified copies of the input sequences.
        This return value preserves the type of both the outer
        and the inner sequences.
    """
    max_length = max(len(l) for l in sequences)
    just = (l + type(l)((filler,))*(max_length-len(l)) for l in sequences)
    return type(sequences)(just)
