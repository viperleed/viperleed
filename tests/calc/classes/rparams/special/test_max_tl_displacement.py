"""Tests for viperleed.calc.classes.rparams.special.max_tl_displacement."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-06-23'
__license__ = 'GPLv3+'

import pytest
from pytest_cases import parametrize

from viperleed.calc.classes.rparams.defaults import NO_VALUE
from viperleed.calc.classes.rparams.special.max_tl_displacement import (
    MaxTLAction,
    MaxTLDisplacement,
    )


class TestMaxTLAction:
    """Tests for the MaxTLAction class."""

    def test_enum_values(self):
        assert MaxTLAction.IGNORE.value == 'ignore'
        assert MaxTLAction.STOP.value == 'stop'
        assert MaxTLAction.REFCALC.value == 'refcalc'

    def test_has_options_property(self):
        assert not MaxTLAction.IGNORE.has_options
        assert not MaxTLAction.STOP.has_options
        assert MaxTLAction.REFCALC.has_options

    def test_enum_identity(self):
        assert MaxTLAction('ignore') is MaxTLAction.IGNORE
        assert MaxTLAction('stop') is MaxTLAction.STOP
        assert MaxTLAction('refcalc') is MaxTLAction.REFCALC


class TestMaxTLDisplacementValid:
    """Collection of tests for valid inputs to MaxTLDisplacement objects."""

    valid_creation = {
        'two floats': (0.15, 0.2),
        'as string': ('1.0', '0.5'),
        'no value': (1.1, NO_VALUE),
        }
    valid_float = {
        'one': (0.15,),
        'two': (0.15, 0.2),
        }
    valid_action = {
        'ignore': (('ignore',), NO_VALUE),
        'stop': (('stop',), NO_VALUE),
        'no duration': (('refcalc',), NO_VALUE),
        'pure float': (('refcalc', '600'), 600),
        'pure float, lt': (('refcalc', '<', '600'), 600),
        'seconds': (('refcalc', '30s'), 30),
        'minutes': (('refcalc', '5m'), 300),
        'hours': (('refcalc', '2h'), 7200),
        }

    def test_valid_initialization(self):
        obj = MaxTLDisplacement(geo=0.5, _vib=0.3)
        assert isinstance(obj, MaxTLDisplacement)
        assert obj.geo == 0.5
        assert obj.vib == 0.3

    def test_default_vib_fallback_to_geo(self):
        obj = MaxTLDisplacement(geo=0.8)
        assert obj.vib == 0.8

    @parametrize('geo,vib', valid_creation.values(), ids=valid_creation)
    def test_type_conversion_and_range(self, geo, vib):
        obj = MaxTLDisplacement(geo=geo, _vib=vib)
        assert isinstance(obj.geo, float)
        if vib is not NO_VALUE:
            assert isinstance(obj.vib, float)

    @parametrize(args=valid_float.values(), ids=valid_float)
    def test_from_value(self, args):
        obj = MaxTLDisplacement.from_value(args)
        assert isinstance(obj, MaxTLDisplacement)
        assert obj.geo == args[0]
        if len(args) > 1:
            assert obj.vib == args[1]

    @parametrize(args=valid_float.values(), ids=valid_float)
    def test_assign_float_values(self, args):
        obj = MaxTLDisplacement(geo=0.5)
        obj.assign_float_values(args)
        assert obj.geo == args[0]
        if len(args) > 1:
            assert obj.vib == args[1]

    def test_assign_named_geo(self):
        obj = MaxTLDisplacement(geo=1.0)
        obj.assign_single_value('geo', 0.7)
        assert obj.geo == 0.7

    def test_assign_named_vib(self):
        obj = MaxTLDisplacement(geo=1.0)
        obj.assign_single_value('vib', 0.6)
        assert obj.vib == 0.6

    @parametrize('values,expect', valid_action.values(), ids=valid_action)
    def test_assign_refcalc(self, values, expect):
        obj = MaxTLDisplacement(geo=1.0)
        obj.assign_action(*values)
        assert obj.action is MaxTLAction[values[0].upper()]
        if obj.action is MaxTLAction.REFCALC:
            assert obj.max_duration == expect


class TestMaxTLDisplacementInvalid:
    """Collection of tests for invalid inputs to MaxTLDisplacement objects."""
    invalid = {
        'not numeric': ('a', TypeError),
        'negative': (-1, ValueError),
        'zero': (0, ValueError),
        }

    @parametrize('value,err', invalid.values(), ids=invalid)
    def test_invalid(self, value, err):
        with pytest.raises(err):
            MaxTLDisplacement(value)


class TestMaxTLDisplacementTooFar:
    """Tests for whether displacement thresholds have been reached."""

    cases = {   # limit, geo.dist., vib.dist., expect
        'geometry': (.1, 1., 0., True),
        'vibration': (.1, 0., 1., True),
        'in range': (.1, 0.05, 0.05, False),
        }

    @parametrize('limit,dgeo,dvib,expect', cases.values(), ids=cases)
    def test_is_too_far(self, limit, dgeo, dvib, expect, mocker):
        atom = mocker.MagicMock()
        atom.distance.return_value = dgeo
        atom.oriState = mocker.MagicMock()
        atom.site.vibamp = {'X': 0.1}
        atom.site.oriState.vibamp = {'X': 0.1+dvib}
        obj = MaxTLDisplacement(limit)
        assert obj.is_too_far(atom) is expect
