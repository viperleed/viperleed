"""Tests for module itertools_utils of viperleed/calc/lib."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-08-06'
__license__ = 'GPLv3+'

from itertools import islice
import sys

import pytest
from pytest_cases import parametrize
from pytest_cases import parametrize_with_cases

from viperleed.calc.lib.itertools_utils import consecutive_groups
from viperleed.calc.lib.itertools_utils import cycle
from viperleed.calc.lib.itertools_utils import n_wise
from viperleed.calc.lib.itertools_utils import pairwise
from viperleed.calc.lib.itertools_utils import threewise

from .cases_itertools_utils import CasesSequence
from .cases_itertools_utils import CasesRaises

PY_310 = (3, 10)
PY_311 = (3, 11)


class ReenteringIterator:
    """An iterator that tracks re-entrance."""

    def __init__(self, reenter_at, iterator_func):
        """Initialize instance.

        Parameters
        ----------
        reenter_at : object
            The item at which this iterator reenters.
        iterator_func : callable
            Will be called with this instance as the only argument.

        Returns
        -------
        None.
        """
        self.count = 0
        self.reenter_at = reenter_at
        self.iterator = iterator_func(self)

    def __iter__(self):
        return self

    def __next__(self):
        self.count +=1
        if self.count in self.reenter_at:
            return next(self.iterator)
        return [self.count]  # new object


class StoppingIterator:
    """An iterator with a maximum number of items."""

    def __init__(self, maxcount):
        """Initialize instance.

        Parameters
        ----------
        maxcount : int
            The number of elements yielded by this iterator.

        Returns
        -------
        None.
        """
        self.count = 0
        self.maxcount = maxcount
        self.iterator = None

    def __call__(self, iterator_func):
        self.count = 0
        self.iterator = iterator_func(self)
        return self.iterator

    def __iter__(self):
        return self

    def __next__(self):
        if self.count >= self.maxcount:
            raise StopIteration
        self.count +=1
        if self.count == 1:
            return next(self.iterator, None)
        return [self.count]  # new object


class TestConsecutiveGroups:
    """Tests for the consecutive_groups iterator."""

    _groups = {
        'list': ([1, 2, 3, 11, 12, 21, 22], [(1, 2, 3), (11, 12), (21, 22)]),
        'single': ([5], [(5,)]),
        'no consecutive': ([1, 3, 5, 7], [(1,), (3,), (5,), (7,)]),
        'all consecutive': (list(range(5)), [tuple(range(5))]),
        'empty': ([], []),
        'large gap': ([1, 2, 50, 51, 52, 100], [(1, 2), (50, 51, 52), (100,)]),
        'unsorted': ([1, 3, 2, 4], [(1,), (3,), (2,), (4,)]),
        'negative': ([-3, -2, -1, 1, 2, 3], [(-3, -2, -1), (1, 2, 3)]),
        }

    @parametrize('iterable,expect', _groups.values(), ids=_groups)
    def test_groups(self, iterable, expect):
        """Check expected outcome for a simple iterable."""
        assert list(consecutive_groups(iterable)) == expect

    def test_infinite(self):
        """Check result of grouping on an infinite iterator."""
        _infinite = cycle([1,2,3], start=1)
        result = list(islice(consecutive_groups(_infinite), 4))
        expect = [(2, 3), (1, 2, 3), (1, 2, 3), (1, 2, 3)]
        assert result == expect

    _fails = {
        # (iter, what you'd intuitively expect, what you actually get)
        # These situations could be handled with custom ordering
        # functions, like more-itertools does.
        'duplicates': ([1, 2, 2, 3, 5, 6, 6, 7],
                       [(1, 2, 2, 3), (5, 6, 6, 7)],
                       [(1, 2), (2, 3), (5, 6), (6, 7)]),
        'floats': ([1.1, 1.2, 1.3, 2.1, 2.2, 3.1],
                   [(1.1, 1.2, 1.3), (2.1, 2.2), (3.1,)],
                   [(1.1,), (1.2,), (1.3,), (2.1,), (2.2,), (3.1,)]),
        }

    @parametrize('iterable,would_like,is_instead', _fails.values(), ids=_fails)
    def test_fails(self, iterable, would_like, is_instead):
        """Check situations that we do not currently support."""
        result = list(consecutive_groups(iterable))
        assert result != would_like
        assert result == is_instead

    _raises = {
        'string': ('abcefg', NotImplementedError),
        'tuples': ([(1, 2, 3), (4, 5, 6)], NotImplementedError),
        }

    @parametrize('iterable,exc', _raises.values(), ids=_raises)
    def test_raises(self, iterable, exc):
        """Check complaints for unsupported situations."""
        with pytest.raises(exc):
            tuple(consecutive_groups(iterable))


class TestCycle:
    """Tests for the extension of itertools' infinite cycle iterator."""

    _simple = {
        'empty': (tuple(), 12, []),
        'range': (range(5), 13, [0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2]),
        'generator': ((2*i for i in range(5)), 8, [0, 2, 4, 6, 8, 0, 2, 4,]),
        'list': ([1, 2, 3], 6, [1, 2, 3, 1, 2, 3]),
        'one item': ([42], 7, [42, 42, 42, 42, 42, 42, 42]),
        'string': ('abc', 8, ['a', 'b', 'c', 'a', 'b', 'c', 'a', 'b']),
        'tuple': (tuple(range(10)),
                  21,
                  [0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
                  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
                  0]),
        }

    @parametrize('sequence,n_items,expect', _simple.values(), ids=_simple)
    def test_simple(self, sequence, n_items, expect):
        """Check a simple call to cycle, without a start."""
        assert list(islice(cycle(sequence), n_items)) == expect

    _start = {
        'start at zero': ([1, 2, 3], 0, [1, 2, 3, 1, 2, 3]),
        'start at one': ([1, 2, 3], 1, [2, 3, 1, 2, 3, 1]),
        'start at two': ([1, 2, 3], 2, [3, 1, 2, 3, 1, 2]),
        'start past length': ([1, 2, 3], 4, [2, 3, 1, 2, 3, 1]),
        'start past length further': ([1, 2, 3], 5, [3, 1, 2, 3, 1, 2]),
        'start past length very large': ([1, 2, 3], 10279, [2, 3, 1, 2, 3, 1]),
        'start negative': ([1, 2, 3], -1, [3, 1, 2, 3, 1, 2]),
        'start negative past length': ([1, 2, 3], -4, [3, 1, 2, 3, 1, 2]),
        'empty': (tuple(), 28, []),
        'generator': (range(4), 37, [1, 2, 3, 0, 1, 2]),
        'one item': ([42], 9237, [42, 42, 42, 42, 42, 42]),
        'string': ('abcd', 2, ['c', 'd', 'a', 'b', 'c', 'd']),
        'tuple': ((1, 2, 3), 1, [2, 3, 1, 2, 3, 1]),
        }

    @parametrize('sequence,start,expect', _start.values(), ids=_start)
    def test_with_start(self, sequence, start, expect):
        """Check correct outcome when a start value is given."""
        assert list(islice(cycle(sequence, start=start), 6)) == expect

    def test_mutable_consumed(self):
        """Check mutation of a sequence after consumption has no impact."""
        seq = [1, 2, 3]
        gen = cycle(seq)
        assert list(islice(gen, 7)) == [1, 2, 3, 1, 2, 3, 1]
        seq.append(4)
        # Notice that it goes on from where it was before
        assert list(islice(gen, 7)) == [2, 3, 1, 2, 3, 1, 2]

    def test_mutable_changed(self):
        """Check mutation of a sequence before consumption has impact."""
        seq = [1, 2, 3]
        gen = cycle(seq, start=1)
        seq.append(4)
        assert list(islice(gen, 7)) == [2, 3, 4, 1, 2, 3, 4]

    @parametrize(start=('a', 2.5))
    def test_raises(self, start):
        """check complaints for wrong start types."""
        with pytest.raises(TypeError):
            list(islice(cycle([1, 2, 3], start=start), 6))


class TestNWise:
    """Tests for the n_wise class."""

    _init = {
        'empty string': ('', 5, []),
        'single letter': ('a', 5, []),
        'two letters': ('ab', 5, []),
        'five letters': ('abcde', 6, []),
        'six letters': ('abcde', 5, [tuple('abcde')]),
        'ten letters': ('abcdefghij', 8,
                         [('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'),
                          ('b', 'c', 'd', 'e', 'f', 'g', 'h', 'i'),
                          ('c', 'd', 'e', 'f', 'g', 'h', 'i', 'j')]),
        'large range': (range(10_000), 6,
                        list(zip(*(range(i, 10_000) for i in range(6))))),
        }

    @parametrize('arg,n_items,expect', _init.values(), ids=_init)
    def test_list(self, arg, n_items, expect):
        """Check expected outcome for list(pairwise(arg))."""
        assert list(n_wise(arg, n_items)) == expect

    _raises = {
        'negative': (-5, ValueError),
        'float': (2.38, ValueError),
        'non numeric': ({}, TypeError),
        }

    @parametrize('n_items,exc', _raises.values(), ids=_raises)
    def test_raises(self, n_items, exc):
        """Check complaints for funny n_wise values."""
        with pytest.raises(exc):
            tuple(n_wise(range(10), n_items))


class TestPairwise:
    """Tests for the pairwise iterator.

    These tests are adapted from the ones of CPython 3.14, as found at
    https://github.com/python/cpython/blob/main/Lib/test/test_itertools.py

    This guarantees that our implementation of pairwise passes the same
    tests as the most recent CPython implementation.
    """

    _init = {
        'empty string': ('', []),
        'single letter': ('a', []),
        'two letters': ('ab', [('a', 'b')]),
        'five letters': ('abcde',
                         [('a', 'b'), ('b', 'c'), ('c', 'd'), ('d', 'e')]),
        'large range': (range(10_000),
                        list(zip(range(10_000), range(1, 10_000)))),
        }

    @parametrize('arg,expect', _init.values(), ids=_init)
    def test_list(self, arg, expect):
        """Check expected outcome for list(pairwise(arg))."""
        assert list(pairwise(arg)) == expect

    _raises = {
        'too few arguments': ((), {}, TypeError),
        'too many arguments': (('abc', 10), {}, TypeError),
        # The next one is in CPython, but we can't really do it
        # in a backwards-compatible manner, as positional-only
        # arguments have been introduced in Python 3.8, and we
        # don't really want to keep such a version granularity
        # if possible.
        # 'keyword arguments': ((), {'iterable': 'abc'}, TypeError),
        }

    @parametrize('args,kwargs,exc', _raises.values(), ids=_raises)
    def test_raises(self, args, kwargs, exc):
        """Check complaints when given *args, **kwargs."""
        with pytest.raises(exc):
            pairwise(*args, **kwargs)

    _raises_iter = {
        **_raises,
        # The CPython implementation raises the following exceptions
        # already at call time. We cannot, but at least we can check
        # that we do once the pairwise is iterated over.
        'non-iterable argument': ((None,), {}, TypeError),
        }

    @parametrize('args,kwargs,exc', _raises_iter.values(), ids=_raises_iter)
    def test_raises_when_consumed(self, args, kwargs, exc):
        """Check complaints when given *args, **kwargs and is iterated."""
        with pytest.raises(exc):
            tuple(pairwise(*args, **kwargs))

    xfail_reason = (
        'Looks like in CPython < 3.10 the recursive call to '
        'StoppingIterator/ReenteringIterator.__next__ that comes from '
        'returning next(self.iterator) causes a ValueError(generator '
        'already executing). This seems like a problem with our current '
        'implementation of pairwise that is not really equivalent to '
        'the one of CPython. This mechanics seems too complicated for us '
        'to actually use it: do not worry about the inconsistency for now.'
        )
    _reenter = {
        'one': ({1}, [(([2], [3]), [4]),
                      ([4], [5])]),
        'two': ({2}, [([1], ([1], [3])),
                      (([1], [3]), [4]),
                      ([4], [5])]),
        'three': ({3}, [([1], [2]),
                        ([2], ([2], [4])),
                        (([2], [4]), [5]),
                        ([5], [6])]),
        'one and two': ({1, 2}, [((([3], [4]), [5]), [6]),
                                 ([6], [7])]),
        'one and three': ({1, 3}, [(([2], ([2], [4])), [5]),
                                   ([5], [6])]),
        'one and four': ({1, 4}, [(([2], [3]), (([2], [3]), [5])),
                                  ((([2], [3]), [5]), [6]),
                                  ([6], [7])]),
        'two and three': ({2, 3}, [([1], ([1], ([1], [4]))),
                                   (([1], ([1], [4])), [5]),
                                   ([5], [6])]),
        }

    @parametrize('reenter_at,expect', _reenter.values(), ids=_reenter)
    @pytest.mark.xfail(sys.version_info < PY_310, reason=xfail_reason)
    def test_reenter(self, reenter_at, expect):
        """Check that reentering pairwise after reenter_at yields expected."""
        iterator = ReenteringIterator(reenter_at, pairwise).iterator
        for item in expect:
            assert next(iterator) == item

    _reenter_finite = {
        1: [],
        2: [],
        3: [],
        4: [(([2], [3]), [4])],
        }

    @parametrize('maxcount,expect',
                 _reenter_finite.items(),
                 ids=_reenter_finite)
    @pytest.mark.xfail(sys.version_info < PY_310, reason=xfail_reason)
    def test_reenter_finite_length(self, maxcount, expect):
        """Check pairwise result when reentering with a max number of items."""
        iterable = StoppingIterator(maxcount)
        iterable(pairwise)
        assert list(iterable.iterator) == expect

    _sequences = (
        '123',
        '',
        range(1000),
        ('do', 1.2),
        range(2000,2200,5),
        )

    @parametrize(seq=_sequences)
    @parametrize_with_cases('seq_type', cases=CasesSequence)
    def test_sequence(self, seq, seq_type):
        """Check correct result of pairwise for non-builtin sequences."""
        as_list = list(seq_type(seq))
        expect = list(zip(as_list, as_list[1:]))
        actual = list(pairwise(seq_type(seq)))
        assert actual == expect

    _raises_sequences = (  # Empty string does not raise always
        '123',
        range(1000),
        ('do', 1.2),
        range(2000,2200,5),
        )

    @parametrize(seq=_raises_sequences)
    @parametrize_with_cases('seq_type', cases=CasesRaises)
    def test_sequence_raises(self, seq, seq_type):
        """Check complaints when trying to iterate over a funny sequence."""
        with pytest.raises(seq_type.exc):
            # The CPython implementation raises even without
            # consuming pairwise. Ours does not, but that's
            # enough.
            tuple(pairwise(seq_type(seq)))


# We don't need many tests for this one, as
# it's only a simple wrapper around n_wise
class TestThreeWise:  # pylint: disable=too-few-public-methods
    """Tests for the threewise function."""

    _init = {
        'empty string': ('', []),
        'single letter': ('a', []),
        'two letters': ('ab', []),
        'three letters': ('abc', [('a', 'b', 'c')]),
        'five letters': ('abcde',
                         [('a', 'b', 'c'), ('b', 'c', 'd'), ('c', 'd', 'e')]),
        'large range': (
            range(10_000),
            list(zip(range(10_000), range(1, 10_000), range(2, 10_000)))
            ),
        }

    @parametrize('arg,expect', _init.values(), ids=_init)
    def test_list(self, arg, expect):
        """Check expected outcome for list(pairwise(arg))."""
        assert list(threewise(arg)) == expect
