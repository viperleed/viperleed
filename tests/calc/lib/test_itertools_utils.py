"""Tests for module itertools_utils of viperleed/calc/lib."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-08-06'
__license__ = 'GPLv3+'

import sys

import pytest
from pytest_cases import parametrize
from pytest_cases import parametrize_with_cases

from viperleed.calc.lib.itertools_utils import pairwise

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
