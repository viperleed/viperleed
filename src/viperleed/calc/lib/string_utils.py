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


import re


def parent_name(dotted_name, remove=''):
    """Return a version of `dotted_name` with the last attribute removed.

    Parameters
    ----------
    dotted_name : str
        The dotted name (e.g., the name of an attribute or a module).
    remove : str, optional
        What should be removed from the right side. It may, or may
        not contain a dot. If it does not begin with a leading dot,
        a single dot is prepended to `remove`. If not given, the
        last dotted portion is removed. Default is an empty string.

    Returns
    -------
    stripped_name : str
        A version of `dotted_name` with [.]`remove` removed from
        the right.
    """
    remove = remove or '.'
    if not remove.startswith('.'):
        remove = '.' + remove
    stripped_name, *_ = dotted_name.rsplit(remove, maxsplit=1)
    return stripped_name


def range_to_str(il):
    """Return a sort-and-compressed string version of an integer sequence.

    Parameters
    ----------
    il : Sequence
        The list of integers to be sorted and compressed.

    Returns
    -------
    range_string : str
        A sorted and compressed form of the integers specified in
        `il`. For example, when called with [1, 6, 4, 5, 2, 8]
        will return '1-2, 4-6, 8'.

    Raises
    ------
    TypeError
        If any of the items in integers is not an int.
    """
    if not all(isinstance(v, int) for v in il):
        t = next(type(v) for v in il if not isinstance(v, int))
        raise TypeError(
            f'range_to_str: expected list of int, found type {t.__name__}'
            )
    sl = sorted(il, reverse=True)
    prev = sl.pop()
    rmin = prev
    out = str(prev)
    while sl:
        v = sl.pop()
        if v == prev:
            continue
        if v - prev == 1:
            prev = v
            continue
        if prev != rmin:
            out += f'-{prev}, {v}'
        else:
            out += f', {v}'
        prev = v
        rmin = v
    if prev != rmin:
        out += f'-{prev}'
    return out


def readIntLine(line, width=3):                                                 # TODO: Probably better ways with list comprehension
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
    out = []
    while line:
        chunk, line = line[:width], line[width:]
        out.append(int(chunk))
    return tuple(out)


def readIntRange(s):
    """Returns a list of integers from a string.

    Parameters
    ----------
    s : str
        The string containing the integers. It may contain multiple,
        space-separated, range specifications in the form i1-i2/i1:i2.
        In this case, the bounds are considered inclusive, i.e., the
        return value contains both i1 and i2.

    Returns
    -------
    list of int
    """
    out = []
    sl = s.split()
    for ss in sl:
        try:
            out.append(int(ss))
        except ValueError:
            try:
                start, stop = split_string_range(ss)
            except ValueError:
                return []

            try:
                out.extend(range(int(start), int(stop)+1))
            except ValueError:
                return []
    return list(dict.fromkeys(out))  # Like set, but keep order


def split_string_range(range_string):
    """Return start and stop as strings from 'start:stop' or 'start-stop'."""
    range_list = []
    if '-' in range_string:
        range_list = range_string.split('-')
    elif ':' in range_string:
        range_list = range_string.split(':')
    try:
        return range_list[0], range_list[-1]
    except IndexError:
        raise ValueError(f'Invalid range string: {range_string}') from None


def strip_comments(line, strip_whitespaces=True):
    """Return the part of `line` to the left of comments."""
    for comment_char in '!#%':
        try:
            line, *_ = line.split(comment_char)
        except ValueError:  # Nothing left to split
            return ''
    return line.strip() if strip_whitespaces else line


UNDERSCORE = '_'

def to_snake_case(other_case):
    """Return a snake_case name from another one."""
    pascal_parts = re.findall(r'[A-Z]+[a-z\d]*', other_case)
    if not pascal_parts:
        return other_case
    small_caps_match = re.match(r'[a-z_]+', other_case)
    parts = (pascal_parts if not small_caps_match
             else (small_caps_match.group(), *pascal_parts))
    if ''.join(parts) != other_case:
        raise ValueError(f'{other_case!r} is not in PascalCase or camelCase')
    snake_case = UNDERSCORE.join(p.lower() for p in parts)
    if UNDERSCORE in other_case and UNDERSCORE*2 in snake_case:
        raise NotImplementedError(f'{other_case!r} is in mixed case')
    return snake_case
