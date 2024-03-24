"""Tests for module layer_cuts of viperleed.calc.classes.rparams.special."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-10-22'
__license__ = 'GPLv3+'

from dataclasses import dataclass

import pytest
from pytest_cases import parametrize, fixture

from viperleed.calc.classes.rparams.special.layer_cuts import (
    LayerCutToken as Cut,
    LayerCutTokenType as CutType,
    LayerCuts
    )

from .....helpers import InfoBase


class TestLayerCutToken:
    """Collection of tests for LayerCutToken objects."""

    valid_tokens = {
        'numeric': (CutType.NUMERIC, 0.73),
        'order >': (CutType.ORDERING, '>'),
        'order <': (CutType.ORDERING, '<'),
        'auto dc': (CutType.AUTO_DC, 1.5),
        'auto dz': (CutType.AUTO_DZ, 1.4),
        'auto w num bounds': (CutType.AUTO_DC, 1.37, 0.2, 0.8),
        }
    invalid_tokens = {
        'numeric not-a-real': ((CutType.NUMERIC, 'abcd'), TypeError),
        'numeric too large': ((CutType.NUMERIC, 1.01), ValueError),
        'numeric too small': ((CutType.NUMERIC, -0.01), ValueError),
        'auto not-a-real': ((CutType.AUTO_DC, 'abcd'), TypeError),
        'order not <>': ((CutType.ORDERING, '!'), ValueError),
        'non auto w bounds': ((CutType.NUMERIC, 0.57, 0.3, 0.4), ValueError),
        'auto w wrong bounds': ((CutType.AUTO_DC, 1.3, '0.2', 0.8), TypeError),
        }

    @fixture(name='make_token', scope='session')
    def fixture_make_token(self):
        """Return a token from a name."""
        def _make(name):
            return Cut(*self.valid_tokens[name])
        return _make

    @parametrize('name,args', valid_tokens.items(), ids=valid_tokens)
    def test_creation(self, name, args, make_token):
        """Check successful creation of a LayerCutToken."""
        token = make_token(name)
        assert token.value == args[1]
        assert token.type_ == args[0]

    @parametrize('token_args,exc', invalid_tokens.values(), ids=invalid_tokens)
    def test_creation_fails(self, token_args, exc):
        """Check instantiation complaints with invalid values."""
        with pytest.raises(exc):
            Cut(*token_args)

    invalid_str = {
        'numeric not-a-real': 'abcd',
        'numeric too large': '1.01',
        'numeric too small': '-0.01',
        'auto not-a-real': 'dz(1.2.3)',
        'order not <>': '!',
        }

    @parametrize(string=invalid_str.values(), ids=invalid_str)
    def test_from_string_invalid(self, string):
        """Check that an INVALID token is created."""
        token = Cut.from_string(string)
        assert token.type_ is CutType.INVALID
        assert str(token) == 'INVALID'

    attrs_info = {'numeric': {'is_auto_cut': False,
                              'is_ordering': False,
                              'is_numeric': True},
                  'order >': {'is_auto_cut': False,
                              'is_ordering': True,
                              'is_numeric': False},
                  'auto dc': {'is_auto_cut': True,
                              'is_auto_cut_dc': True,
                              'is_ordering': False,
                              'is_numeric': False},
                  'auto dz': {'is_auto_cut': True,
                              'is_auto_cut_dc': False,
                              'is_ordering': False,
                              'is_numeric': False},}

    @parametrize(name=attrs_info)
    def test_attrs(self, name, make_token):
        """Check that properties of a token are as expected."""
        token = make_token(name)
        attrs = self.attrs_info[name]
        for attr, value in attrs.items():
            assert getattr(token, attr) == value

    valid_bounds = {  # bounds, extra_lower, expected
        'both':  ((0.12, 0.87), (), (0.12, 0.87)),
        'upper': ((None, 0.38), (), (0, 0.38)),
        'lower': ((0.09, None), (), (0.09, 1)),
        'both extra below':  ((0.12, 0.87), (0.01, 0.05), (0.12, 0.87)),
        'both extra above':  ((0.12, 0.87), (0.01, 0.23), (0.23, 0.87)),
        'lower extra below':  ((0.09, None), (0.01, 0.05), (0.09, 1)),
        'lower extra above':  ((0.09, None), (0.01, 0.17), (0.17, 1)),
        }

    @parametrize('bounds,extra,expect', valid_bounds.values(),
                 ids=valid_bounds)
    def test_get_bounds_valid(self, bounds, extra, expect, make_token):
        """Check correctness of the bounds of an auto-cut token."""
        token = make_token('auto dc')
        with_bounds = token.with_bounds(*bounds)
        assert with_bounds.get_bounds(extra) == expect

    invalid_bounds = {  # token_name, bounds, extra, exc_set, exc_get
        'min > max': ('auto dc', (0.35, 0.2), (), None, RuntimeError),
        'extra too large':  ('auto dc', (None, 0.38), (0.39, 0.86),
                             None, ValueError),
        'wrong token': ('numeric', (0.1, 0.2), (), TypeError, AttributeError),
        }

    @parametrize('name,bounds,extra,exc_set,exc_get', invalid_bounds.values(),
                 ids=invalid_bounds)
    def test_get_bounds_invalid(self, name, bounds, extra,
                                exc_set, exc_get, make_token):
        """Check correctness of the bounds of an auto-cut token."""
        token = make_token(name)
        if exc_set is None:
            token = token.with_bounds(*bounds)
        else:
            with pytest.raises(exc_set):
                token.with_bounds(*bounds)
        with pytest.raises(exc_get):
            token.get_bounds(extra)

    def test_set_bounds_twice_fails(self, make_token):
        """Check that it is not possible to set bounds twice."""
        token = make_token('auto dc')
        bounded = token.with_bounds(0.1, None)
        assert bounded is not token
        with pytest.raises(RuntimeError):
            # Must fail also when trying to set the other bound
            bounded.with_bounds(None, 0.2)

    def test_set_bounds_wrong_type(self, make_token):
        """Check that it is not possible to set bounds twice."""
        token = make_token('auto dc')
        wrong_bound = make_token('auto dz')
        with pytest.raises(TypeError):
            # Must fail also when trying to set the other bound
            token.with_bounds(wrong_bound, None)

    equals = {
        'numbers eq': (Cut.make_numeric(0.39), Cut.make_numeric(0.39), True),
        'numbers neq': (Cut.make_numeric(0.39), Cut.make_numeric(0.15), False),
        'numbers eq string': (Cut.make_numeric(0.39), '0.39', True),
        'numbers eq float': (Cut.make_numeric(0.23), 0.23, True),
        'ordering eq': (Cut(CutType.ORDERING, '>'),
                        Cut(CutType.ORDERING, '>'),
                        True),
        'ordering neq': (Cut(CutType.ORDERING, '>'),
                         Cut(CutType.ORDERING, '<'),
                         False),
        'dc no bounds': (Cut(CutType.AUTO_DC, 1.34),
                         Cut(CutType.AUTO_DC, 1.34),
                         True),
        'dc one bound': (
            Cut(CutType.AUTO_DC, 1.34, Cut.make_numeric(0.1), None),
            Cut(CutType.AUTO_DC, 1.34),
            False
            ),
        'dc both bound': (
            Cut(CutType.AUTO_DC, 1.34, Cut.make_numeric(0.1), None),
            Cut(CutType.AUTO_DC, 1.34, None, Cut.make_numeric(0.75)),
            False
            ),
        }

    @parametrize('first,second,result', equals.values(), ids=equals)
    def test_equal(self, first, second, result):
        """Check the expected result of an equality test."""
        assert (first == second) is result


@dataclass(repr=False)
class LayerCutInfo(InfoBase):
    """Expected data for a given LayerCuts input."""

    value : object     # To be used to make a LayerCuts
    str_tokens : str   # str(cuts._tokens)
    str_self : str     # str(cuts)
    format_4f: str     # format(cuts, '.4f'), or f'{cuts:.4f}'
    len_ : int         # len(cuts)


class TestLayerCuts:
    """Collection of tests for LayerCuts objects."""

    @staticmethod
    def _make_tokens_string(cuts):
        """Return a string version of cuts._tokens."""
        # pylint: disable=protected-access
        # We really do want to access that, as this is the best
        # way to test correct interpretation (and we don't need
        # a public interface for this).
        return ' '.join(str(t) for t in cuts._tokens)

    def _validate_expected_info(self, cuts, info, subtests):
        """Check that cuts has the values expected from info."""
        with subtests.test('bool True'):
            assert cuts
        with subtests.test('str(cuts._tokens)'):
            assert self._make_tokens_string(cuts) == info.str_tokens
        with subtests.test('str(cuts)'):
            assert str(cuts) == info.str_self
        with subtests.test('format(cuts, .4f)'):
            assert f'{cuts:.4f}' == info.format_4f
        with subtests.test('len(cuts)'):
            assert len(cuts) == info.len_

    str_valid = {
        'numeric': LayerCutInfo('0.5 0.3', '0.5 0.3', '0.5 0.3',
                                '0.5000 0.3000', 2),
        'dc only': LayerCutInfo('dc(3)', 'dc(3.0)', 'dc(3.0)', 'dc(3.0)', 1),
        'number + dc': LayerCutInfo('0.25 < dc(1.4)', '0.25 dc(1.4)',
                                    '0.25 < dc(1.4)', '0.2500 < dc(1.4)', 2),
        'dz + number': LayerCutInfo('dz(0.9)>0.3', '0.3 dz(0.9)',
                                    'dz(0.9) > 0.3', 'dz(0.9) > 0.3000', 2),
        'multiple': LayerCutInfo(
            '0.15 0.1 0.2<dz(0.3) <0.4 0.5 <dc(0.1)<0.6',
            '0.15 0.1 0.2 dz(0.3) 0.4 0.5 dc(0.1) 0.6',
            '0.15 0.1 0.2 < dz(0.3) < 0.4 0.5 < dc(0.1) < 0.6',
            '0.1500 0.1000 0.2000 < dz(0.3) < 0.4000 '
            '0.5000 < dc(0.1) < 0.6000',
            8
            ),
        }
    str_invalid = {
        'both <>':               '0.1 < dc(0.3) > 0.25',
        '< at end':              '< dc(0.3)',
        'two <<':                '0.25 << dz(1.3)',
        'missing <, number':     '0.25 dc(0.3)',
        'missing <, dc num dz':  'dc(1.4) 0.4 dz(1.3)',
        'missing <, num dc num': '0.12 dc(1.4) 0.4',
        'missing <, dc dz':      'dc(1.4) dz(1.3)',
        'missing num':           'dc(1.2) < dz(1.5)',
        'invalid':               'dc(1.2',
        'not float cutoff':      'dz(1.2.3)',
        '< no dc':               '0.9 > 0.8 > 0.5',
        'bounds swapped':        '0.9 < dc(1.3) < 0.4',
        'numeric within bounds': '0.3 0.1 0.2<dz(1)<0.8 0.9 0.75',
        }

    @parametrize('info', str_valid.values(), ids=str_valid)
    def test_from_string(self, info, subtests):
        """Check successful creation of a LayerCuts from string."""
        cuts = LayerCuts.from_string(info.value)
        self._validate_expected_info(cuts, info, subtests)

    @parametrize(string=str_invalid.values(), ids=str_invalid)
    def test_from_string_invalid(self, string):
        """Check complaints when creating from a string with syntax errors."""
        with pytest.raises(ValueError):
            LayerCuts.from_string(string)

    def test_bool_false(self):
        """Check False-ness of an empty LayerCuts."""
        assert not LayerCuts()

    def test_iteration_less_than(self):
        """Check that iterating over a LayerCuts gives the expected result."""
        tokens = [Cut(CutType.NUMERIC, 0.4),
                  Cut(CutType.ORDERING, '<'),
                  Cut(CutType.AUTO_DC, 1.5)]
        expected = (tokens[0], tokens[-1].with_bounds(0.4, None))
        layer_cuts = LayerCuts(*tokens)
        assert tuple(layer_cuts) == expected

    def test_iteration_larger_than(self):
        """Check that iterating over a LayerCuts gives the expected result."""
        tokens = [Cut(CutType.NUMERIC, 0.4),
                  Cut(CutType.ORDERING, '>'),
                  Cut(CutType.AUTO_DC, 1.5)]
        layer_cuts = LayerCuts(*tokens)
        expected = (tokens[-1].with_bounds(None, 0.4), tokens[0])
        assert tuple(layer_cuts) == expected

    def test_creation_wrong_types(self):
        """Check complaints when trying to create with non-LayerCutToken."""
        with pytest.raises(TypeError):
            LayerCuts(0.1, 0.2, 'dc(1.5)')

    def test_update_from_sequence(self):
        """Check that filling LayerCuts from a Sequence works as expected."""
        sequence = (0.27, '<', 'dc(1.54)', '<', 0.85, Cut.make_numeric(0.97))
        properties = ('is_numeric', 'is_auto_cut_dc',
                      'is_numeric', 'is_numeric')
        values = (0.27, 1.54, 0.85, 0.97)
        cuts = LayerCuts()
        cuts.update_from_sequence(sequence)
        for prop_getter, value, item in zip(properties, values, cuts):
            assert getattr(item, prop_getter)
            assert item.value == value
        assert self._make_tokens_string(cuts) == '0.27 dc(1.54) 0.85 0.97'
        assert str(cuts) == '0.27 < dc(1.54) < 0.85 0.97'
        assert len(cuts) == 4

    def test_update_from_sequence_fails(self):
        """Check that filling LayerCuts from a Sequence works as expected."""
        sequence = (0.27, '<', 'dc(1.54)', '<', 0.85,
                    Cut.make_numeric(0.97), [1.3, 4.5])
        cuts = LayerCuts()
        with pytest.raises(TypeError):  # last item
            cuts.update_from_sequence(sequence)

    as_cuts_valid = {
        **str_valid,
        'sequence' : LayerCutInfo(
            (0.27, '<', 'dc(1.54)', '<', 0.85, Cut.make_numeric(0.97)),
            '0.27 dc(1.54) 0.85 0.97', '0.27 < dc(1.54) < 0.85 0.97',
            '0.2700 < dc(1.54) < 0.8500 0.9700', 4
            ),
        'layer_cuts': LayerCutInfo(LayerCuts.from_string('dz(1.2)'), 'dz(1.2)',
                                   'dz(1.2)', 'dz(1.2)', 1),
        }
    as_cuts_invalid = {
        **{k: (v, ValueError) for k, v in str_invalid.items()},
        'seq invalid': ((0.27, '<', 'dc(1.54)', '<', 0.85,
                         Cut.make_numeric(0.97), [1.3, 4.5]),
                        TypeError),
        'wrong type': (CutType.NUMERIC, TypeError),
        }

    @parametrize(info=as_cuts_valid.values(), ids=as_cuts_valid)
    def test_as_layer_cuts_valid(self, info, subtests):
        """Check correct creation from LayerCuts, string, or sequence."""
        cuts = LayerCuts.as_layer_cuts(info.value)
        with subtests.test('is layer cut'):
            assert isinstance(cuts, LayerCuts)
        self._validate_expected_info(cuts, info, subtests)

    @parametrize('arg,exc', as_cuts_invalid.values(), ids=as_cuts_invalid)
    def test_as_layer_cuts_invalid(self, arg, exc):
        """Check correct creation from LayerCuts, string, or sequence."""
        with pytest.raises(exc):
            LayerCuts.as_layer_cuts(arg)
