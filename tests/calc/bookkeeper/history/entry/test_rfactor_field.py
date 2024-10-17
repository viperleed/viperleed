"""Tests for rfactor_field of viperleed.calc.bookkeeper.history.entry."""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2024-09-04'
__license__ = 'GPLv3+'

import pytest
from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.bookkeeper.history.entry.field import DefaultMessage
from viperleed.calc.bookkeeper.history.entry.field import EmptyField
from viperleed.calc.bookkeeper.history.entry.field import MissingField
from viperleed.calc.bookkeeper.history.entry.rfactor_field import RFactorField
from viperleed.calc.bookkeeper.history.entry.rfactor_field import RRefField
from viperleed.calc.bookkeeper.history.entry.rfactor_field import RSuperField
from viperleed.calc.bookkeeper.history.errors import EntrySyntaxError

from .test_field import _TestFieldUtils


@fixture(name='rfactor_field')
@parametrize(field_cls=(RFactorField, RRefField, RSuperField))
def factory_rfactor_field(field_cls, make_concrete_field_instance):
    """Return an instance of an RRefField (sub)class with a value."""
    return make_concrete_field_instance(field_cls)


class TestRFactorField(_TestFieldUtils):
    """Base class for tests of RFactorField and its subclasses."""

    _init = {
        'empty': (
            '',
            {'value': EmptyField, 'was_understood': False, '_value_str': None,
             'is_empty': True},
            ),
        'components': (
            '1.2345 (0.1234 / 0.5678)',
            {'value': '1.2345 (0.1234 / 0.5678)', 'was_understood': True,
             '_value_str': '1.2345 (0.1234 / 0.5678)'},
            ),
        'extra space, components': (
            '  1.2345 ( 0.1234   / 0.5678 )  ',
            {'value': '1.2345 ( 0.1234   / 0.5678 )', 'was_understood': True,
             '_value_str': '1.2345 ( 0.1234   / 0.5678 )'},
            ),
        'extra space, total only': (
            '  1.2345  ',
            {'value': 1.2345, 'was_understood': True, '_value_str': '1.2345'},
            ),
        'float': (
            1.2345,
            {'value': 1.2345, 'was_understood': True, '_value_str': '1.2345'},
            ),
        'float large': (
            2.1,
            {'value': 2.1, 'was_understood': False, '_value_str': None},
            ),
        'fractional out of range': (
            '1.2345 (0.1234 / 2.5678)',
            {'value': '1.2345 (0.1234 / 2.5678)', 'was_understood': False,
             '_value_str': None},
            ),
        'integer': (
            1,
            {'value': 1, 'was_understood': False, '_value_str': None},
            ),
        'invalid': (
            'invalid',
            {'value': 'invalid', 'was_understood': False},
            ),
        'missing': (
            MissingField,
            {'value': MissingField, 'is_empty': False, 'is_missing': True,
             'was_understood': True},  # Is not mandatory
            ),
        'multiple components': (
            '1.2345 (0.1234 / 0.5678) (0.2345 / 0.6789)',
            {'value': '1.2345 (0.1234 / 0.5678) (0.2345 / 0.6789)',
             'was_understood': False},
            ),
        'no slash': (
            '1.2345 (0.1234 0.5678)',
            {'value': '1.2345 (0.1234 0.5678)', 'was_understood': False},
            ),
        'none': (
            None,
            {'value': EmptyField, 'was_understood': False, 'is_empty': True},
            ),
        'unmatched parentheses': (
            '1.2345 (0.1234 / 0.5678',
            {'value': '1.2345 (0.1234 / 0.5678', 'was_understood': False},
            ),
        }

    @parametrize('value,attrs', _init.values(), ids=_init)
    def test_init(self, value, attrs, rfactor_field):
        """Check result of initialization and value checking."""
        self.check_attrs(rfactor_field, attrs, value)

    _bounds = {
        'two, float': (
            2.0,
            {'value': 2.0, '_value_str': '2.0000'},
            ),
        'two, components': (
            '2.000(2.0  /2.00)  ',
            {'value': '2.000(2.0  /2.00)', '_value_str': '2.000(2.0  /2.00)'},
            ),
        'two, total only': (
            '2.000',
            {'value': 2.0, '_value_str': '2.000'},
            ),
        'zero, float': (
            0.0,
            {'value': 0.0, '_value_str': '0.0000'},
            ),
        'zero, components': (
            '0.0 (0.00/0.0)',
            {'value': '0.0 (0.00/0.0)', '_value_str': '0.0 (0.00/0.0)'},
            ),
        'zero, total only': (
            '0.00',
            {'value': 0.0, '_value_str': '0.00'},
            ),
        'below zero, float': (
            -1.1e-5,
            {'value': -1.1e-5, 'was_understood': False},
            ),
        'above two, float': (
            2.000011,
            {'value': 2.000011, 'was_understood': False},
            ),
        'above two, components': (
            '1.234 (2.000011, 2.0)',
            {'value': '1.234 (2.000011, 2.0)', 'was_understood': False},
            ),
        }

    @parametrize('value,attrs', _bounds.values(), ids=_bounds)
    def test_bounds(self, value, attrs, rfactor_field):
        """Check values at the bounds of validity."""
        self.test_init(value, attrs, rfactor_field)

    @parametrize(invalid=(1, '123', []))
    def test_check_float_not_float(self, invalid, rfactor_field):
        """Check complaints when _check_float_value gets a non-float."""
        field = rfactor_field(invalid)
        reason = DefaultMessage.NOT_FLOAT.value
        with pytest.raises(EntrySyntaxError, match=reason):
            # pylint: disable-next=protected-access       # OK in tests
            field._check_float_value()
