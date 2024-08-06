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


def readIntRange(s):
    """Takes a string, returns a list of integers. If the string is a single
    int or space-separated ints, the return value is a list containing only
    those int. If the string contains i1-i2 or i1:i2, the list is
    range(i1, i2+1), i.e. contains both i1 and i2."""
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
    """Return start and stop as strings from "start:stop" or "start-stop"."""
    range_list = []
    if '-' in range_string:
        range_list = range_string.split('-')
    elif ':' in range_string:
        range_list = range_string.split(':')
    try:
        return range_list[0], range_list[-1]
    except IndexError:
        raise ValueError(f'Invalid range string: {range_string}') from None
