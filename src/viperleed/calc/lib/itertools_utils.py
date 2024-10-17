"""Module itertools_utils of viperleed.calc.lib.

Collects functions working with and returning iterables and iterators.
Extends the itertools standard-library module, and provides backward-
compatible itertools functions for earlier python version.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-08-06'
__license__ = 'GPLv3+'

from collections import deque
from itertools import cycle as iter_cycle
from itertools import groupby
from itertools import islice
import sys

PY_313 = (3, 13)

if sys.version_info >= PY_313:
    # batched is available already in 3.12, but has no 'strict' keyword
    from itertools import batched
else:
    batched = None  # pylint: disable=invalid-name   # OK not UPPERCASE

try:
    from itertools import pairwise   # Python >= 3.10
except ImportError:
    pairwise = None  # Defined later


if not batched:
    # Use a version adapted from the itertools documentation under
    # (by removing the walrus operator in order to support Py3.7 too)
    # https://docs.python.org/3.12/library/itertools.html#itertools.batched
    # pylint: disable-next=invalid-name   # Same signature as Python3.X
    def batched(iterable, n, *, strict=False):
        """Yield `n`-long tuples from iterable."""
        if n < 1:
            raise ValueError('n must be at least one')
        iterator = iter(iterable)
        batch = tuple(islice(iterator, n))
        while batch:
            if strict and len(batch) != n:
                raise ValueError('batched(): incomplete batch')
            yield batch
            batch = tuple(islice(iterator, n))


# The next one is adapted from more-itertools, which uses a recipe from
# Python2: https://docs.python.org/2.6/library/itertools.html#examples
# The more-itertools version also has a second keyword argument to pass
# along a 'ordering' callable. We can pull that version if we need this
# for anything more complex than integers.
def consecutive_groups(iterable):
    """Yield groups of consecutive items.

    Parameters
    ----------
    iterable : Iterable
        The iterable from which items should be returned in a grouped
        fashion. It is assumed that `iterable` is already sorted in
        increasing order.

    Yields
    ------
    consecutive_items : Iterable
        Consecutive items in iterable. It shares it source with
        `iterable`. When an an output group is advanced, the previous
        group is no longer available unless its elements are copied
        (e.g., into a list).

        >>> iterable = [1, 2, 3, 11, 12, 21, 22]
        >>> saved_groups = []
        >>> for group in consecutive_groups(iterable):
        ...     saved_groups.append(list(group))  # Copy group elements
        >>> saved_groups
        [[1, 2, 3], [11, 12], [21, 22]]
    """
    # The trick is that the lambda uses the difference between
    # each item's index and the item itself. For example, with
    # [    1,      2,      3,      11,      12,      21,      22]
    # enumerate gives
    # [(0, 1), (1, 2), (2, 3), (4, 11), (5, 12), (6, 21), (7, 22)]
    # hence, the keys used for grouping are
    # [   -1,     -1,     -1,      -7,      -7,     -15,     -15]
    # so that the groups are split correctly (itertools.groupby
    # makes a new group every time the key changes value).
    grouped = groupby(
        enumerate(iterable),
        key=lambda ind_and_item: ind_and_item[0] - ind_and_item[1]
        )
    groups = (g for _, g in grouped)  # Don't care about the keys
    try:  # pylint: disable=too-many-try-statements
        # We can't have fewer statements, as the thing that raises
        # is the iteration in the 'for'. One could consume the whole
        # 'groups' beforehand, but that's not great as we would drop
        # support for infinite sequences.
        for group in groups:
            # Remove the enumerate indices
            yield tuple(item for _, item in group)
    except TypeError:
        raise NotImplementedError('consecutive_groups supports '
                                  'only integers.') from None


def cycle(iterable, start=0):
    """Return a generator that cycles though `iterable` beginning at start."""
    if not isinstance(start, int):
        raise TypeError(f'start must be int, not {type(start).__name__!r}')
    if iterable:
        try:
            start %= len(iterable)
        except TypeError:
            # Not a sequence. Do it the hard way, i.e., slicing from
            # start. May be not performing well for large start values
            # compared to the number of elements.
            pass
    _cycled = iter_cycle(iterable)
    return islice(_cycled, start, None)


def n_wise(iterable, n_items):
    """Yield `n`-tuples of items from `iterable`."""
    if n_items < 2:  # pylint: disable=magic-value-comparison  # Clear
        raise ValueError('n_wise needs at least n_items=2')
    iterable = iter(iterable)  # In case it's a Sequence
    items = deque(islice(iterable, n_items - 1), n_items)
    for item in iterable:
        items.append(item)
        yield tuple(items)


if not pairwise:
    def pairwise(iterable):
        """Yield pairs of items from `iterable` as tuples."""
        yield from n_wise(iterable, 2)


def threewise(iterable):
    """Yield triplets of items from an `iterable`."""
    yield from n_wise(iterable, 3)
