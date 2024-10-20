"""Module sequence_utils of viperleed.calc.lib.

Defines functions useful to manipulate sequences."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-08-25'
__license__ = 'GPLv3+'

from collections.abc import Sequence
from operator import itemgetter

from viperleed.calc.lib.itertools_utils import pairwise


# Adapted from https://stackoverflow.com/questions/45382174
def conditional_sort(sequence, skip, key=None):
    """Return a sorted list of items in sequence, skipping some.

    Parameters
    ----------
    sequence : Sequence
        The sequence whose items should be sorted.
    skip : callable
        A callable that will be called on each item. Should return True
        for those items that the user wants to stay in their place and
        not be sorted. Only those items for which skip returns False
        are sorted.
    key : callable, optional
        Same meaning as for the sorted built-in function.

    Returns
    -------
    sorted_sequence : type(Sequence)
        Sequence with items selected by `skip` in unchanged position,
        all the others sorted according to `key`.
    """
    if not isinstance(sequence, Sequence):
        raise TypeError('Not a sequence')
    if not callable(skip):
        raise TypeError('skip must be callable')
    cls = type(sequence)

    sorted_items = iter(sorted(
        (item for item in sequence if not skip(item)),
        key=key,
        ))
    return cls(item if skip(item) else next(sorted_items)
               for item in sequence)


def max_diff(sequence):
    """Return the index and value of the largest difference of two items.

    Parameters
    ----------
    sequence : Sequence
        The sequence of items for which the largest difference is
        returned. Must have at least two items.

    Returns
    -------
    ind : int
        The first index such that sequence[ind] - sequence[ind-1]
        is the largest among the differences between two consecutive
        items in `sequence`.
    maxdiff : object
        The value of the largest difference.

    Raises
    ------
    ValueError
        If `sequence` is too short to compute a difference.
    """
    if len(sequence) < 2:  # pylint: disable=magic-value-comparison
        raise ValueError(f'Argument {sequence!r} has only one item')
    ind_and_diff = (
        (i + 1, second - first)
        for i, (first, second) in enumerate(pairwise(sequence))
        )
    return max(ind_and_diff, key=itemgetter(1))


def recombine_items(sequence, com):
    """Merge items that start/end with `com` with the previous/next one."""
    if not sequence:
        return []
    iterable = iter(sequence)
    buffer = next(iterable)
    newlist = []
    for item in iterable:
        if buffer.endswith(com) or item.startswith(com):
            buffer += item
        else:
            newlist.append(buffer)
            buffer = item
    newlist.append(buffer)
    return newlist
