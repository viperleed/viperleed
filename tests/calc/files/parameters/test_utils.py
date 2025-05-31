"""Tests for module utils of viperleed.calc.files.parameters."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2023-10-15'
__license__ = 'GPLv3+'

import pytest
from pytest_cases import parametrize

from viperleed.calc.files.parameters.utils import Assignment
from viperleed.calc.files.parameters.utils import NumericBounds as Bounds

from ....helpers import not_raises


_FLOAT_ZERO_FIVE = Bounds(range_=(0, 5))
_INT_FIVE_TEN_MODULO = Bounds(type_=int, range_=(5, 10),
                              out_of_range_event='modulo')


class TestNumericBounds:
    """Collection of tests for the NumericBounds class."""

    _out_of_range_fmt = {
        'both_limits': (_FLOAT_ZERO_FIVE, 10,
                        'Value 10 is not in range [0, 5]'),
        'lower_limit': (Bounds(type_=int, range_=(-1, None)), -20,
                        'Value -20 is less than -1'),
        'lower_limit_open': (
            Bounds(type_=int, range_=(-1, None), accept_limits=(False, False)),
            -20, 'Value -20 is less than or equal to -1'
            ),
        'modulo': (_INT_FIVE_TEN_MODULO, 3, 'Value 3 is not in range [5, 10]'),
        'upper_limit': (Bounds(range_=(None, 5)), 10,
                        'Value 10 is larger than 5'),
        'upper_limit_open': (
            Bounds(range_=(None, 5), accept_limits=(False, False)),
            10, 'Value 10 is larger than or equal to 5'
            ),
        }

    @parametrize('bounds, value, message',
                 _out_of_range_fmt.values(),
                 ids=_out_of_range_fmt)
    def test_format_out_of_range(self, bounds, value, message):
        """Check correctness of out-of-range message."""
        assert bounds.format_out_of_range(value) == message

    _in_range = {  # (bounds, value, result)
        'both_ok': (_FLOAT_ZERO_FIVE, 3, (True, True)),
        'above_upper': (Bounds(type_=int, range_=(1, 10)), 15, (True, False)),
        }

    @parametrize('bounds,value,result', _in_range.values(), ids=_in_range)
    def test_is_in_range(self, bounds, value, result):
        """Check correct identification of in-/out-of-range conditions."""
        assert bounds.is_in_range(value) == result

    _make_in_range = {  # bounds, value, result
        'in_range': (_FLOAT_ZERO_FIVE, 2.3, 2.3),
        'coerce_high': (
            Bounds(type_=int, range_=(1, 10), out_of_range_event='coerce'),
            15, 10
            ),
        'coerce low, upper unbound' : (
            Bounds(type_=int, range_=(12, None), out_of_range_event='coerce'),
            5, 12
            ),
        'modulo_in_range': (_INT_FIVE_TEN_MODULO, 7, 7),
        'modulo_too_high': (_INT_FIVE_TEN_MODULO, 16, 6),
        'modulo_at_max': (_INT_FIVE_TEN_MODULO, 10, 10),                        # TODO: should treat better the cases at bounds! Must depend on accept_limits.
        'modulo_at_min': (_INT_FIVE_TEN_MODULO, 5, 5),                          # TODO: should treat better the cases at bounds! Must depend on accept_limits.
        'modulo_too_low': (_INT_FIVE_TEN_MODULO, 3, 8),
        'float_modulo': (Bounds(range_=(5, 10), out_of_range_event='modulo'),
                         2.4, 7.4),
        }

    @parametrize('bounds, value, result', _make_in_range.values(),
                 ids=_make_in_range)
    def test_make_in_range(self, bounds, value, result):
        """Check correctness of forcing a value to be in-range."""
        assert bounds.make_in_range(value) == result

    _range_settings = {  # (bounds, accept_limits, (min, max), closed_range)
        'int_only_low': (
            Bounds(type_=int, range_=(5, 10), accept_limits=(True, False)),
            (True, False), (5, 10), (5, 9)
            ),
        'int_both': (
            Bounds(type_=int, range_=(1, 10)),
            (True, True), (1, 10), (1, 10)
            ),
        'float_only_high': (
            Bounds(range_=(5, 10), accept_limits=(False, True)),
            (False, True), (5, 10), (5.0001, 10)
            ),
        'float_both': (_FLOAT_ZERO_FIVE, (True, True), (0, 5), (0, 5)),
        'unbound': (Bounds(), (True, True), (None, None), (None, None)),
        }

    @parametrize('bounds,limits,minmax,closed_range',
                 _range_settings.values(), ids=_range_settings)
    def test_range(self,  bounds, limits, minmax, closed_range):
        """Check that range-related attributes have the expected values."""
        assert bounds.accept_limits == limits
        assert (bounds.min, bounds.max) == minmax
        assert bounds.closed_range == closed_range

    def test_range_max_only(self):
        """Check attributes for a lower-unbound NumericBounds."""
        bounds = Bounds(type_=int, range_=(None, 10))
        assert not bounds.unlimited
        assert bounds.max == 10

    def test_range_min_only(self):
        """Check attributes for an upper-unbound NumericBounds."""
        bounds = Bounds(type_=float, range_=(5, None))
        assert not bounds.unlimited
        assert bounds.min == 5

    _true_property = {    # (bounds, attribute to be tested)
        'coerce': (
            Bounds(type_=int, range_=(1, 10), out_of_range_event='coerce'),
            'coerce'
            ),
        'fail': (Bounds(type_=float, range_=(0, 5), out_of_range_event='fail'),
                 'fail'),
        'unbound_float' : (Bounds(type_=float, range_=(None, None)),
                             'unlimited'),
        'unbound_int': (Bounds(type_=int), 'unlimited'),
        }

    @parametrize('bounds, attr', _true_property.values(), ids=_true_property)
    def test_property_is_true(self, bounds, attr):
        """Check that a specific attribute/property is True-thy."""
        assert getattr(bounds, attr)


class TestNumericBoundsRaises:
    """Tests for failing conditions occurring with NumericBounds."""

    _invalid_args = {
        'type': {'type_': str, 'range_': (0, 5)},
        'event': {'range_': (1, 10), 'out_of_range_event': 'invalid_event'},
        'modulo_unlimited': {'type_': int, 'range_': (None, 5),
                             'out_of_range_event': 'modulo'},
        'modulo_accept_lims': {'range_': (0, 5),
                               'out_of_range_event': 'modulo',
                               'accept_limits': (False, False)},
        }

    @parametrize(kwargs=_invalid_args.values(), ids=_invalid_args)
    def test_invalid_init_args(self, kwargs):
        """Check complaints with invalid initialization arguments."""
        with pytest.raises(ValueError):
            Bounds(**kwargs)

    def test_make_in_range_fail(self):
        """Check complaints when making in-range with fail-type event."""
        bounds = _FLOAT_ZERO_FIVE
        with pytest.raises(RuntimeError):
            bounds.make_in_range(20)


class TestAssignment:
    """Collection of tests for Assignment objects."""

    def test_creation(self):
        """Check attributes of a successfully created Assignment."""
        assignment = Assignment(values_str='value1 value2', parameter='PARAM')
        assert assignment.parameter == 'PARAM'
        assert assignment.value == 'value1'
        assert assignment.values == ('value1', 'value2')

    _empty_values = {
        'string values': {'values_str': '', 'parameter': 'P'},
        'sequence values': {'values_str': [], 'parameter': 'P'},
        'spaces value': {'values_str': '   ', 'parameter': 'P'},
        }

    @parametrize(empty_kwargs=_empty_values.values(), ids=_empty_values)
    def test_creation_empty_is_valid(self, empty_kwargs):
        """Check that passing some empty inputs is acceptable."""
        with not_raises(Exception):
            assignment = Assignment(**empty_kwargs)
        assert not assignment.values

    def test_flag(self):
        """Check correct interpretation of the first flag."""
        assignment = Assignment(values_str='value',
                                parameter='PARAM',
                                flags_str='--flag1 --flag2')
        assert assignment.flag == '--flag1'

    def test_other_flags(self):
        """Check correct interpretation of flags beyond the first one."""
        assignment = Assignment(
            values_str='value',
            parameter='PARAM',
            flags_str=('--flag1', '--flag2', '--flag3'),
            )
        assert assignment.other_flags == ('--flag2', '--flag3')

    def test_other_values(self):
        """Check correct interpretation of values beyond the first one."""
        assignment = Assignment(values_str=('value1', 'value2', 'value3'),
                                parameter='PARAM')
        assert assignment.other_values == ('value2', 'value3')

    _values_to_unpack = {
        'str': ('v1 v2 v3', ('v1', 'v2', 'v3')),
        'sequence': (['--flag1', 'flag2', 'f3'], ('--flag1', 'flag2', 'f3')),
        }

    @parametrize('to_unpack,expected', _values_to_unpack.values(),
                 ids=_values_to_unpack)
    def test_unpack_side(self, to_unpack, expected):
        """Check correct unpacking of the left/right side of the '=' sign."""
        # pylint: disable-next=protected-access
        assert Assignment._unpack_assignment_side(to_unpack) == expected


class TestAssignmentRaises:
    """Collection of tests for error conditions foe Assignment objects."""

    def test_creation_no_parameter(self):
        """Check complaints when initializing with one empty argument."""
        with pytest.raises(ValueError):
            Assignment(values_str='value1 value2', parameter='')

    _wrong_type = {
        'values': {'values_str': 123, 'parameter': 'P'},
        'flags': {'values_str': 'value', 'parameter': 'P', 'flags_str': 123},
        'parameter': {'values_str': 'value', 'parameter': 123},
        }

    @parametrize(kwargs=_wrong_type.values(), ids=_wrong_type)
    def test_assignment_creation_non_string_values(self, kwargs):
        """Check complaints when initializing with one wrong type."""
        with pytest.raises(TypeError):
            Assignment(**kwargs)

    def test_unpack_assignment_side_type_error(self):
        """Check complaints when trying to unpack the wrong type."""
        with pytest.raises(TypeError):
            # pylint: disable-next=protected-access
            Assignment._unpack_assignment_side(123)
