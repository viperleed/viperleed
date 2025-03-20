"""Module string_utils of viperleed.calc.lib.

Collects functions useful for manipulating strings and converting
them from/to other formats. Most of these functions used to be in
calc.lib.base, originally from 2019-06-13.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-06-13'
__license__ = 'GPLv3+'

import itertools
import re

import numpy as np

from viperleed.calc.lib.itertools_utils import batched
from viperleed.calc.lib.itertools_utils import consecutive_groups


def harvard_commas(*items, sep='and'):
    """Return a Harvard-comma-separated string of the items in `sequence`."""
    if len(items) > 2:  # pylint: disable=magic-value-comparison
        commas = ', '.join(str(i) for i in items[:-1])
        return commas + f', {sep} {items[-1]}'
    return f' {sep} '.join(str(i) for i in items)


def parent_name(dotted_name, remove=''):
    """Return a version of `dotted_name` with the last attribute removed.

    Parameters
    ----------
    dotted_name : str
        The dotted name (e.g., the name of an attribute or a module).
    remove : str, optional
        What should be removed from the right side. It may, or may
        not contain a dot. If it does not begin with a leading dot,
        a single dot is prepended to `remove`. Notice that any
        portion of `dotted_name` to the right of `remove` will be
        removed! If not given, the last dotted portion is removed.
        Default is an empty string.

    Returns
    -------
    stripped_name : str
        A version of `dotted_name` with [.]`remove`[rest] removed from
        the right.
    """
    remove = remove or '.'
    if not remove.startswith('.'):
        remove = '.' + remove
    stripped_name, *_ = dotted_name.rsplit(remove, maxsplit=1)
    return stripped_name


def range_to_str(integers, sep=', ', range_sep='-'):
    """Return a sort-and-compressed string version of an integer iterable.

    Parameters
    ----------
    integers : Iterable
        The list of integers to be sorted and compressed. This function
        will consume one-pass iterators.
    sep : str, optional
        The separator to use if multiple non-contiguous ranges are
        found in `integers`. Default is ', '.
    range_sep : str, optional
        The separator to use between the first and last item of
        contiguous ranges in `integers`. Default is '-'.

    Returns
    -------
    range_string : str
        A sorted and compressed form of the integers specified in
        `integers`. For example, when called with [1, 6, 4, 5, 2, 8]
        will return '1-2, 4-6, 8'.

    Raises
    ------
    TypeError
        If any of the items in integers is not an int.
    """
    integers, type_check = itertools.tee(integers, 2)
    not_an_int = next((type(v) for v in type_check if not isinstance(v, int)),
                      None)
    if not_an_int:
        raise TypeError('range_to_str: expected list of int, '
                        f'found type {not_an_int.__name__!r}')
    # Remove duplicates and sort input, as consecutive_groups
    # does not handle these correctly
    sorted_integers = sorted(dict.fromkeys(integers))
    first_and_last_items = ((group[0], group[-1])
                            for group in consecutive_groups(sorted_integers))
    return sep.join(str(start) if start == end else f'{start}{range_sep}{end}'
                    for start, end in first_and_last_items)


def read_int_line(line, width=3):
    """Read an (arbitrary length) line of integers with fixed width.

    Parameters
    ----------
    line : str
        The line to interpret
    width : integer, optional
        The width of each integer. The default is 3.

    Returns
    -------
    tuple of int
    """
    line = line.rstrip()
    return tuple(int(''.join(batch)) for batch in batched(line, width))


def read_int_range(ranges):
    """Return a list of integers from a string.

    Parameters
    ----------
    ranges : str
        The string containing the integers. It may contain multiple,
        space-separated range specifications in the form i1-i2/i1:i2.
        In this case, the bounds are considered inclusive, i.e., the
        return value contains both i1 and i2.

    Returns
    -------
    list of int

    Raises
    ------
    ValueError
        If `ranges` cannot be interpreted as a string
        defining one or more ranges of integers.
    """
    out = []
    for int_or_range in ranges.split():
        try:
            out.append(int(int_or_range))
        except ValueError:
            try:
                start, stop = split_string_range(int_or_range)
            except ValueError as exc:
                raise ValueError(f'Could not interpret {int_or_range!r} '
                                 'as a single integer or a range of '
                                 'integers') from exc
            try:
                out.extend(range(int(start), int(stop)+1))
            except ValueError:
                raise ValueError(
                    f'Could not interpret {start!r} or {stop!r} '
                    f'as integer values in range {int_or_range!r}'
                    ) from None
    integers = list(dict.fromkeys(out))  # Like set, but keep order
    if not integers:
        raise ValueError(f'No integers in {ranges!r}')
    return integers


def read_vector(string, ucell=None):
    """Return a 3-item vector in Cartesian coordinates from a string.

    Parameters
    ----------
    string : str
        The string version of the vector to return. May have form
        'xyz[f1 f2 f3]', 'abc[f1 f2 f3]', or just '[f1 f2 f3]'.
    ucell : numpy.ndarray, optional
        The unit cell to use to transform a fractional vector into a
        Cartesian one. Unit vectors are rows, i.e., a, b, c == ucell
        This argument is only used if `string` has the 'abc[f1 f2 f3]'
        form, and it is mandatory in this case.

    Returns
    -------
    numpy.ndarray

    Raises
    ------
    ValueError
        If `string` has an invalid format.
    TypeError
        If `string` has format 'abc[f1 f2 f3]', but no `ucell` is given.
    """
    # pylint: disable-next=magic-value-comparison  # Clear enough
    _fractional = 'abc' in string
    if _fractional and ucell is None:
        raise TypeError(f'ucell is mandatory for {string!r} '
                        'in fractional form.')
    _three_floats = r'\s+'.join(r'(-?\d+(?:\.\d+)?)' for _ in range(3))
    match_ = re.match(fr'\s*(?:xyz|abc)?\[\s*{_three_floats}\s*\]', string)
    if not match_:
        raise ValueError(
            f'Invalid format for {string!r}. Expected \'[f1 f2 f3]\', '
            '\'xyz[f1 f2 f3]\', or \'abc[f1 f2 f3]\'.'
            )
    try:
        vector = np.array(match_.groups()).astype(float)
    except ValueError:
        raise ValueError(f'Failed to convert {match_.groups()} '
                         'to numeric values.') from None
    return vector.dot(ucell) if _fractional else vector


def rsplit_once(string, sep):
    """Split `string` once from the right at `sep`."""
    try:
        return string.rsplit(sep, maxsplit=1)
    except AttributeError:
        raise TypeError(f'Invalid {type(string).__name__!r}. '
                        'Expected str.') from None


def split_string_range(range_string):
    """Return start and stop as strings from 'start:stop' or 'start-stop'."""
    try:
        sep = next(c for c in '-:' if c in range_string)
    except StopIteration:
        raise ValueError(f'Invalid range string: {range_string}') from None
    range_list = range_string.split(sep)
    return range_list[0], range_list[-1]


def strip_comments(line, strip_whitespaces=True):
    """Return the part of `line` to the left of comments."""
    for comment_char in '!#%':
        try:
            line, *_ = line.split(comment_char)
        except AttributeError:  # Not a string
            raise TypeError(f'Invalid line type {type(line).__name__!r}. '
                            'Expected str.') from None
        if not line.strip():
            break
    return line.strip() if strip_whitespaces else line


UNDERSCORE = '_'

def to_snake_case(other_case):
    """Return a snake_case name from another one."""
    if re.search(r'[^\w\d]', other_case):
        # Contains some funny characters.
        raise ValueError(f'{other_case!r} is not a valid snake_case, '
                         'camelCase, or PascalCase name')
    parts = (p for p in re.split(r'([A-Z]\d*(?:[a-z]+\d*)+)', other_case) if p)
    snake_case = UNDERSCORE.join(p.lower() for p in parts)
    if UNDERSCORE in other_case and UNDERSCORE*2 in snake_case:
        raise NotImplementedError(f'{other_case!r} is in mixed case')
    return snake_case
