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


def recombineListElements(llist, com):
    """Merge items that start/end with `com` with the next/previous ones."""
    i = 0
    newlist = llist[:]
    while i < len(newlist)-1:
        if newlist[i][-1] == com or newlist[i+1][0] == com:
            newlist[i] += newlist.pop(i+1)
        else:
            i += 1
    return newlist
