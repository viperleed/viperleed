"""Tests for viperleed.calc.classes.rparams.special.r_factor_type."""

__authors__ = ('Florian Kraushofer (@fkraushofer)',)
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2026-10-08'
__license__ = 'GPLv3+'

import pytest
from pytest_cases import parametrize

from viperleed.calc.classes.rparams.rparams import Rparams
from viperleed.calc.classes.rparams.special.r_factor_type import RFactorType

_ALIASES = {
    RFactorType.PENDRY: ('pendry', 'pendry_r', 'rp', 'pend', 'p',
                         'r_pendry', 'pendry r-factor', 'r_p',
                         'PENDRY', ' Pendry R ', 'pendry-r', '1', 1),
    RFactorType.R2: ('R2', 'r2', 'r squared', 'R²', 'r^2', 'R_2',
                     'r2 factor', 'R SQUARED', '2', 2),
    RFactorType.ZJ: ('zj', 'zanazzi jona', 'zanazzi-jona', 'zanazzi_jona',
                     'zannazi', 'zannazi jona', 'zannazi-jona', 'z',
                     'zj factor', 'ZJ', 'Zanazzi Jona', 'zanazzijona',
                     '3', 3),
    RFactorType.SMOOTH: ('smooth', 'smoothed', 'smooth pendry',
                         'smooth-pendry', 'schmid', 's', 'rs', 'r_s',
                         'r_smooth', 'Smooth_Pendry',
                         'SMOOTH', ' s ', '4', 4),
    }
_VALID = [(alias, member)
          for member, aliases in _ALIASES.items()
          for alias in aliases]
_TENSERLEED_CODES = {  # The integer TensErLEED uses for each type
    'pendry': 1,
    'r2': 2,
    'zj': 3,
    'smooth': 4,
    }


@parametrize('value,expect', _VALID)
def test_from_alias(value, expect):
    """Check that `value` is interpreted as R-factor type `expect`."""
    assert RFactorType(value) is expect


@parametrize(value=('nonsense', 'r1', 'pendryy', '', '0', '5', '²', 0, 5,
                    None))
def test_invalid(value):
    """Check complaints for an invalid R-factor type."""
    with pytest.raises(ValueError, match='not a valid RFactorType'):
        RFactorType(value)


@parametrize('name,code', _TENSERLEED_CODES.items())
def test_tensorleed_code(name, code):
    """Check the integer code passed to TensErLEED for each type."""
    assert RFactorType(name).value == code


def test_member_is_its_own_value():
    """Check that an RFactorType is turned into itself."""
    assert all(RFactorType(member) is member for member in RFactorType)


def test_rparams_default():
    """Check that the default R_FACTOR_TYPE is the Pendry R-factor."""
    assert Rparams().R_FACTOR_TYPE is RFactorType.PENDRY
