"""Tests for module search_cull of viperleed.calc.classes.rparams.special."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-12-14'
__license__ = 'GPLv3+'

import pytest
from pytest_cases import fixture, parametrize

from viperleed.calc.classes.rparams.special.search_cull import SearchCull
from viperleed.calc.classes.rparams.special.search_cull import SearchCullType


@fixture(name='value_and_type')
def factory_value_and_type():
    """Return the value and type_ attributes of a SearchCull."""
    def _make(cull):
        return cull.value, cull.type_
    return _make


class TestSearchCull:
    """Collection of tests for SearchCull objects."""

    valid = {  # args, expected
        'int': ((5,), (5, SearchCullType.UNSET)),
        'int, genetic': ((5, 'genetic'), (5, SearchCullType.GENETIC)),
        'float, clone': ((0.5, 'clone'), (0.5, SearchCullType.CLONE)),
        }
    invalid = {  # args, exception
        'neg': ((-1,), ValueError),
        'not a float': (('invalid',), TypeError),
        'invalid type_': ((5, 'invalid'), ValueError),
        'too large float': ((2.3,), ValueError),
        }

    @parametrize('args,expected', valid.values(), ids=valid)
    def test_valid(self, args, expected, value_and_type):
        """Check correct initialization outcome."""
        cull = SearchCull(*args)
        assert value_and_type(cull) == expected

    @parametrize('args,exc', invalid.values(), ids=invalid)
    def test_invalid(self, args, exc):
        """Check complaints for invalid initialization values."""
        with pytest.raises(exc):
            SearchCull(*args)

    frac_and_type = {  # is_fractional, has_type
        'int': (False, False),
        'int, genetic': (False, True),
        'float, clone': (True, True)
        }

    @parametrize('key,expected', frac_and_type.items(), ids=frac_and_type)
    def test_is_fractional_has_type(self, key, expected):
        """Check correct result of is_fractional and has_type."""
        cull = SearchCull(*self.valid[key][0])
        assert cull.is_fractional is expected[0]
        assert cull.has_type is expected[1]

    bool_ = {True: 2, False: 0}

    @parametrize('result,val', bool_.items(), ids=bool_)
    def test_bool(self, val, result):
        """Check correct truth value of SearchCull."""
        cull = SearchCull(val)
        assert bool(cull) is result

    def test_repr(self):
        """Check correct contents of repr(cull)."""
        cull = SearchCull(3, 'random')
        repr_ = repr(cull)
        assert repr_.startswith('SearchCull(3')
        assert repr_.endswith('.RANDOM)')

    individuals = {  # value, result_out_of_100
        int: (5, 5),
        float: (0.3, 30),
        }
    individuals_invalid = {  # cull, population, exc
        'too many': (20, 19, ValueError),
        'non int population': (20, 'invalid', ValueError),
        'neg population': (0.5, -5, ValueError),
        }

    @parametrize('val,expected', individuals.values(), ids=individuals)
    def test_individuals(self, val, expected):
        """Check correctness of nr_individuals."""
        cull = SearchCull(val)
        assert cull.nr_individuals(100) == expected

    @parametrize('val,population,exc', individuals_invalid.values(),
                 ids=individuals_invalid)
    def test_individuals_fails(self, val, population, exc):
        """Check complaints for invalid number of culled individuals."""
        cull = SearchCull(val)
        with pytest.raises(exc):
            cull.nr_individuals(population)


class TestSearchCullType:
    """Collection of tests for the SearchCullType enum."""

    def test_properties(self):
        """Check correctness of GENETIC cull type."""
        type_ = SearchCullType.GENETIC
        assert type_.is_genetic is True
        assert type_.is_clone is False
        assert type_.is_random is False

    def test_unset(self):
        """Check correctness of UNSET cull type."""
        type_ = SearchCullType.UNSET
        assert type_.is_genetic is False
        assert type_.is_clone is False
        assert type_.is_random is False

    def test_comparison(self):
        """Check selected comparisons."""
        assert SearchCullType.GENETIC == SearchCullType.GENETIC
        assert SearchCullType.GENETIC != SearchCullType.CLONE
        # pylint: disable=magic-value-comparison
        assert SearchCullType.GENETIC != 'invalid_type'
