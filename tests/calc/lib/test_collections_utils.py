"""Tests for module collections_utils of viperleed.calc.lib."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-09-17'
__license__ = 'GPLv3+'

from collections import namedtuple
import copy
import gc  # The garbage collector

from pytest_cases import parametrize

from viperleed.calc.lib.collections_utils import IdentitySet

test_obj = object()  # To use the exact same object

# The next one is needed as Python (at least) 3.12 caches tuple ids
MockImmutable = namedtuple('MockImmutable', ('value',))


class TestIdentitySet:
    """Tests for the IdentitySet class."""

    _init = {
        'empty': ((), 0),
        'one item': ((object(),), 1),
        'more items': (tuple(object() for _ in range(5)), 5),
        'same values, immutable': ((MockImmutable(2), MockImmutable(2)), 2),
        'same values, mutable': (([1, 2, 3], [1, 2, 3]), 2),
        'duplicates': ((test_obj, test_obj), 1),
        'duplicates, None': ((None, None), 1),
        }

    @parametrize('items,n_items', _init.values(), ids=_init)
    def test_initialization(self, items, n_items):
        """Check correct initialization from a bunch of items."""
        set_ = IdentitySet(items)
        assert len(set_) == n_items
        assert all(i in set_ for i in items)
        for item in items:
            set_.add(item)  # These are now duplicates
        assert len(set_) == n_items

    @parametrize('items,n_items', _init.values(), ids=_init)
    def test_add(self, items, n_items):
        """Check correct addition of items one at a time."""
        set_ = IdentitySet()
        for item in items:
            set_.add(item)
            assert item in set_
        assert len(set_) == n_items

    def test_add_duplicate(self):
        """Test correct addition of a duplicate."""
        dupe = object()
        set_ = IdentitySet()
        for _ in range(5):
            set_.add(dupe)
        assert len(set_) == 1
        assert dupe in set_

    def test_add_object_again_after_removal(self):
        """Check correct addition of an item after its removal."""
        item = object()
        set_ = IdentitySet([item])
        assert item in set_
        set_.discard(item)
        assert item not in set_
        set_.add(item)
        assert item in set_

    _contains = {  # added, not added
        'different objects': (object(), object(), False),
        'same object': (test_obj, test_obj, True),
        'same object, None': (None, None, True),
        'same object, type': (object, object, True),
        'equal objects': (MockImmutable(2), MockImmutable(2), False),
        }

    @parametrize('add,not_add,contained', _contains.values(), ids=_contains)
    def test_contains(self, add, not_add, contained):
        """Check membership after adding an item and not another one."""
        set_ = IdentitySet([add])
        assert add in set_
        other_contained = not_add in set_
        assert other_contained == contained

    def test_deepcopy_gives_different_objects(self):
        """Check that deep-copying an IdentitySet gives different objects."""
        # This is a proxy test for the behavior one would expect when
        # IdentitySet objects are used, for example, in multiprocessing
        items = object(), object(), object()
        set_ = IdentitySet(items)
        set_copy = copy.deepcopy(set_)
        assert set_copy != set_

    _discard = {  # added, not_added
        'empty': ((), (object(),)),
        'different objects': ((object(), 1), (object(), -5)),
        'same object': ((test_obj,), (test_obj, -5)),
        'same object, None': ((None,), (None, None)),
        'equal objects': ((1, 2, 3, 4), (1, 2, 3, 4)),
        }

    @parametrize('add,not_add', _discard.values(), ids=_discard)
    def test_discard(self, add, not_add):
        """Check correct discarding of existent and non-existent items."""
        set_ = IdentitySet(add)
        len_ = len(add)
        for item in add:
            assert item in set_
            assert len(set_) == len_
            set_.discard(item)
            len_ -= 1
            assert item not in set_
        assert not set_

        for item in not_add:
            assert item not in set_
            set_.discard(item)
            assert item not in set_
        assert not set_

    @parametrize('items,_', _init.values(), ids=_init)
    def test_iteration(self, items, _):
        """Check that iteration gives back the original objects."""
        set_ = IdentitySet(items)
        iter_items = list(set_)
        assert all(i in iter_items for i in items)

    def test_large_number_of_objects(self):
        """Check correct handling of a large number of objects."""
        n_items = 10000
        items = [object() for _ in range(n_items)]
        set_ = IdentitySet(items)
        assert len(set_) == n_items
        for item in items:
            assert item in set_
        set_.discard(items[0])  # Remove one object
        assert items[0] not in set_
        assert len(set_) == n_items - 1

    def test_object_id_reuse(self):
        """Check correct behavior when object ids are reused after deletion."""
        item = 42  # Use a small integer, as CPython caches those
        set_ = IdentitySet([item])
        item_id = id(item)

        # Remove the object and trigger garbage collection
        set_.discard(item)
        del item      # Now reference count should be zero.
        gc.collect()  # Force garbage collection to free the object

        new_item = 42  # Same value, very likely the same id

        assert id(new_item) == item_id
        assert new_item not in set_

    def test_repr(self):
        """Check __repr__ method."""
        items = object(), object()
        set_ = IdentitySet(items)
        expect = f'IdentitySet{items}'
        assert repr(set_) == expect

    _update = {  # The first one is used for __init__
        'duplicates by identity': (([test_obj], [test_obj, object()]), 2),
        'empty, one': (([object()], []), 1),
        'empty, multiple': (([object()], [], (), set()), 1),
        'large iterable': (([], [object() for _ in range(1000)]), 1000),
        'mixed iterables': (
            ([object()], [object()], {object()}, (object(),)),
            4,
            ),
        'multiple iterables': (([object()], [object()], [object()]), 3),
        'overlapping items': (([test_obj], [test_obj, None], [None, 3]), 3),
        'same content, different objects': (([1, 2], [1, 2]), 2),
        }

    @parametrize('others,n_items', _update.values(), ids=_update)
    def test_update(self, others, n_items):
        """Check expected outcome of updating an IdentitySet."""
        set_ = IdentitySet(others[0])
        set_.update(*others[1:])
        items = [i for other in others for i in other]
        assert all(i in set_ for i in items)
        assert len(set_) == n_items

    def test_update_with_self(self):
        """Check correct updating using IdentitySet as an iterable."""
        items = object(), object()
        set_ = IdentitySet(items)
        self.test_update((set_, set_, set_), len(set_))
