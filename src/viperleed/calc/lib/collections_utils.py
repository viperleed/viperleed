"""Module collections_utils of viperleed.calc.lib.

Defines custom collections of items, based on the ABCs of the
collections.abc standard-library module.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-09-17'
__license__ = 'GPLv3+'

from collections.abc import MutableSet


class IdentitySet(MutableSet):
    """A set checking item presence via identity rather than equality.

    Notice that items are also sorted by insertion order. Moreover,
    non-hashable items can also be stored safely.

    Instances of this class are NOT SAFE in a multiprocessing context,
    as objects are typically copied between process instances. Since
    IdentitySet uses identity-based comparisons, the copied objects
    are identified as different. In code:
    >>> set_ = IndentitySet(objects)
    >>> copied_set = deepcopy(set_)
    >>> assert copied_set != set_
    """

    # Identity comparison happens internally by storing a map of
    # {id(item): item} to ensure O(1) item access. Notice that, while
    # normally one should not rely on object ids, it is safe to do so
    # here. The only sensible risk with object ids is that they may be
    # reused IF objects are garbage-collected. However, since the map
    # in this IdentitySet retains hard references to the items, no
    # garbage collection can happen unless an item is removed from
    # this IdentitySet (or the IdentitySet instance is destroyed). This
    # makes usage of ids safe (except multiprocessing).

    __slots__ = '_items', '__weakref__'

    def __init__(self, items=None):
        """Initialize instance from an optional iterable."""
        items = items or tuple()
        self._items = {id(item): item for item in items}

    def __contains__(self, item):
        """Return whether the `item` object is in this IdentitySet."""
        return id(item) in self._items

    def __iter__(self):
        """Return an iterator of items in this IdentitySet."""
        return iter(self._items.values())

    def __len__(self):
        """Return the number of items in this IdentitySet."""
        return len(self._items)

    def __repr__(self):
        """Return a representation string for this IdentitySet."""
        return f'{type(self).__name__}{tuple(self._items.values())}'

    def add(self, value):
        """Add `value` to this IdentitySet."""
        self._items[id(value)] = value

    def discard(self, value):
        """Silently remove `value` from this IdentitySet."""
        self._items.pop(id(value), None)

    def update(self, *others):
        """Add elements from other iterables to this IdentitySet."""
        for iterable in others:
            self |= iterable
