"""Module sequence_utils of viperleed.calc.lib.

Defines functions useful to manipulate sequences."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-08-25'
__license__ = 'GPLv3+'


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
