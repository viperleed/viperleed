"""Tests for module sequence_utils of viperleed.calc.lib."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-08-25'
__license__ = 'GPLv3+'

import pytest
from pytest_cases import parametrize

from viperleed.calc.lib.sequence_utils import recombineListElements


class TestRecombineListElements:
    """Tests for the recombineListElements function."""

    _valid = {
        'at end': (['apple', 'banana-', '-cherry'], '-',
                   ['apple', 'banana--cherry']),
        'at_start': (['-apple', 'banana', 'cherry'], '-',
                     ['-apple', 'banana', 'cherry']),
        'empty': ([], '-', []),
        'left': (['a', 'b', '-c', 'd'], '-', ['a', 'b-c', 'd']),
        'right': (['a', 'b-', 'c', 'd'], '-', ['a', 'b-c', 'd']),
        'mixed': (['a', 'b-', '-c', 'd-', '-e'], '-', ['a', 'b--c', 'd--e']),
        'multiple': (['a', 'b-', '-c', '-d', 'e'], '-', ['a', 'b--c-d', 'e']),
        'no combinations': (['apple', 'banana', 'cherry'], '-',
                            ['apple', 'banana', 'cherry']),
        'one item': (['a'], '-', ['a']),
        'only combinations': (['-', '-', '-'], '-', ['---']),
        }

    @parametrize('seq,sep,expect', _valid.values(), ids=_valid)
    def test_valid(self, seq, sep, expect):
        """Check correct outcome with acceptable arguments."""
        assert recombineListElements(seq, sep) == expect

    def test_raises_not_string(self):
        """Check complaints with non-string items."""
        with pytest.raises(TypeError):
            recombineListElements(['a', 1], '-')
