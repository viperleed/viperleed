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
from itertools import islice
from itertools import cycle as iter_cycle


try:
    from itertools import pairwise   # Python >= 3.10
except ImportError:
    pairwise = None  # Defined later


def cycle(iterable, start=0):
    """Return a generator that cycles though `iterable` beginning at start."""
    if not isinstance(start, int):
        raise TypeError(f'start must be int, not {type(start).__name__!r}')
    if iterable:
        try:
            start %= len(iterable)
        except AttributeError:
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
