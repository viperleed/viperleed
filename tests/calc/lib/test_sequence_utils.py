"""Tests for module sequence_utils of viperleed.calc.lib."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-08-25'
__license__ = 'GPLv3+'

import pytest
from pytest_cases import parametrize

from viperleed.calc.lib.sequence_utils import conditional_sort
from viperleed.calc.lib.sequence_utils import max_diff
from viperleed.calc.lib.sequence_utils import recombine_items


class TestConditionalSort:
    """Tests for the conditional_sort function."""

    _valid = {
        'empty': ([], lambda x: False, None, []),
        'mixed types': (['cherry', 'apple', 8, 'banana', 2],
                        lambda x: isinstance(x, str),
                        None,
                        ['cherry', 'apple', 2, 'banana', 8]),
        'no skipping': ([5, 3, 4, 2, 1], lambda x: False, None,
                        [1, 2, 3, 4, 5]),
        # pylint: disable-next=magic-value-comparison
        'skip one': ((5, 3, 4, 2, 1), lambda x: x == 3, None,
                     (1, 3, 2, 4, 5)),
        'with key': (
            # pylint: disable-next=magic-value-comparison
            ['apple_long', 'banana', 'cherry'], lambda x: x == 'banana', len,
            ['cherry', 'banana', 'apple_long']),
        }

    @parametrize('seq,skip,key,expect', _valid.values(), ids=_valid)
    def test_valid(self, seq, skip, key, expect):
        """Check correctness of result with valid arguments."""
        seq_sorted = conditional_sort(seq, skip, key=key)
        assert seq_sorted == expect

    _raises = {
        'not a sequence': (None, lambda x: False),
        'skip not callable': ([1, 2, 3], 'not_callable'),
        }

    @parametrize(args=_raises.values(), ids=_raises)
    def test_raises(self, args):
        """Check complaints for invalid arguments."""
        with pytest.raises(TypeError):
            conditional_sort(*args)


class TestMaxDiff:
    """Tests for the max_diff function."""

    _valid = {
        'large': (list(range(10_000)), (1, 1)),
        'negative': ([5, -5, 0, -10, 10], (4, 20)),
        'no diff': ([1, 1, 1, 1], (1, 0)),
        'same diff multiple': ([1, 2, 4, 5, 7, 8], (2, 2)),
        'sorted': ([1, 3, 6, 10], (3, 4)),
        'two items': ([3, 5], (1, 2)),
        }

    @parametrize('seq,expect', _valid.values(), ids=_valid)
    def test_valid(self, seq, expect):
        """Check correct result with acceptable arguments."""
        assert max_diff(seq) == expect

    _raises = {
        'empty': ([], ValueError),
        'single item': ([5], ValueError),
        }

    @parametrize('seq,exc', _raises.values(), ids=_raises)
    def test_raises(self, seq, exc):
        """Check complaints with inappropriate arguments."""
        with pytest.raises(exc):
            max_diff(seq)


class TestRecombineItems:
    """Tests for the recombine_items function."""

    _valid = {
        'at end': (['apple', 'banana-', '-cherry'], '-',
                   ['apple', 'banana--cherry']),
        'at_start': (['-apple', 'banana', 'cherry'], '-',
                     ['-apple', 'banana', 'cherry']),
        'empty': ({}, '-', []),
        'left': (['a', 'b', '-c', 'd'], '-', ['a', 'b-c', 'd']),
        'right': (['a', 'b-', 'c', 'd'], '-', ['a', 'b-c', 'd']),
        'mixed': (['a', 'b-', '-c', 'd-', '-e'], '-', ['a', 'b--c', 'd--e']),
        'multiple': (('a', 'b-', '-c', '-d', 'e'), '-', ['a', 'b--c-d', 'e']),
        'no combinations': (['apple', 'banana', 'cherry'], '-',
                            ['apple', 'banana', 'cherry']),
        'one item': (['a'], '-', ['a']),
        'only combinations': (['-', '-', '-'], '-', ['---']),
        }

    @parametrize('seq,sep,expect', _valid.values(), ids=_valid)
    def test_valid(self, seq, sep, expect):
        """Check correct outcome with acceptable arguments."""
        assert recombine_items(seq, sep) == expect

    def test_raises_not_string(self):
        """Check complaints with non-string items."""
        with pytest.raises(AttributeError):
            recombine_items(['a', 1], '-')
