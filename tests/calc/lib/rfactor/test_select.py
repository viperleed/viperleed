"""Tests for module select of viperleed.calc.lib.rfactor."""

__authors__ = ('Alexandra M. Imre (@alexmiame)',
               'Florian Kraushofer (@fkraushofer)')
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2026-10-08'
__license__ = 'GPLv3+'

import pytest
from pytest_cases import parametrize

from viperleed.calc.classes.rparams.special.r_factor_type import RFactorType
from viperleed.calc.lib.rfactor import pendry
from viperleed.calc.lib.rfactor import r_1
from viperleed.calc.lib.rfactor import r_2
from viperleed.calc.lib.rfactor import select
from viperleed.calc.lib.rfactor import smooth
from viperleed.calc.lib.rfactor import zannazi_jona

_FUNCTIONS = {
    RFactorType.PENDRY: pendry.r_pendry,
    RFactorType.R2: r_2.r_2,
    RFactorType.ZJ: zannazi_jona.r_zjj,
    RFactorType.SMOOTH: smooth.r_s,
    }

# All the R_FACTOR_TYPE aliases of each R-factor type
_ALIASES = {
    RFactorType.PENDRY: ('pendry', 'pendry_r', 'rp', 'pend', 'p',
                         'r_pendry', 'pendry r-factor', '1'),
    RFactorType.R2: ('R2', 'r2', 'r squared', 'r2 factor', '2'),
    RFactorType.ZJ: ('zj', 'zanazzi jona', 'zanazzi-jona', 'zanazzi_jona',
                     'zannazi', 'zannazi jona', 'z', 'zj factor', '3'),
    RFactorType.SMOOTH: ('smooth', 'smoothed', 'smooth pendry',
                         'smooth-pendry', 'schmid', 's', 'rs', 'r_smooth',
                         '4'),
    }

# The names that from_ase.rfactor_from_csv accepted for which_r
# before RFactorType became the only source of R-factor names
_FROM_ASE_NAMES = {
    pendry.r_pendry: ('pendry', 'r_p', 'r_pendry', 'rp',
                      'pendry r-factor', 'p'),
    r_1.r_1: ('r1', 'r_1', 'r1 factor'),
    r_2.r_2: ('r2', 'r_2', 'r2 factor'),
    smooth.r_s: ('s', 'rs', 'r_s', 'r_smooth', 'smooth', 'schmid'),
    zannazi_jona.r_zjj: ('zj', 'zj factor', 'zannazi', 'zannazi jona',
                         'zannazi-jona'),
    }
_R1_NAMES = _FROM_ASE_NAMES[r_1.r_1]


def _cases(aliases):
    """Yield (name, key) pairs from a {key: names} dictionary."""
    for key, names in aliases.items():
        for name in names:
            yield name, key


class TestSelectRfactor:
    """Tests for the select_rfactor function."""

    @parametrize(r_factor_type=RFactorType)
    def test_member(self, r_factor_type):
        """Check the function selected for an RFactorType."""
        selected = select.select_rfactor(r_factor_type)
        assert selected is _FUNCTIONS[r_factor_type]

    @parametrize('name,r_factor_type', _cases(_ALIASES))
    def test_alias(self, name, r_factor_type):
        """Check the function selected for an R_FACTOR_TYPE alias."""
        selected = select.select_rfactor(name)
        assert selected is _FUNCTIONS[r_factor_type]

    @parametrize(name=_R1_NAMES)
    def test_no_r1(self, name):
        """Check that the R1 R factor cannot be selected."""
        with pytest.raises(ValueError, match='not a valid RFactorType'):
            select.select_rfactor(name)


class TestSelectRfactorByName:
    """Tests for the select_rfactor_by_name function."""

    @parametrize('name,func', _cases(_FROM_ASE_NAMES))
    def test_from_ase_names(self, name, func):
        """Check that the names from_ase accepted before still work."""
        assert select.select_rfactor_by_name(name) is func

    @parametrize('name,r_factor_type',
                 [(n, t) for n, t in _cases(_ALIASES) if not n.isdigit()])
    def test_alias(self, name, r_factor_type):
        """Check the function selected for an R_FACTOR_TYPE alias."""
        selected = select.select_rfactor_by_name(name)
        assert selected is _FUNCTIONS[r_factor_type]

    @parametrize(name=('  RP  ', ' R1 ', 'R_1'))
    def test_strips_and_lowercases(self, name):
        """Check that names are case-insensitive and stripped."""
        expected = r_1.r_1 if '1' in name else pendry.r_pendry
        assert select.select_rfactor_by_name(name) is expected

    @parametrize(name=('1', '2', '3', '4', ' 1 '))
    def test_rejects_codes(self, name):
        """Check that integer codes are not taken as R-factor names."""
        with pytest.raises(ValueError, match='Unknown R-factor name'):
            select.select_rfactor_by_name(name)

    def test_unknown(self):
        """Check complaints for an unknown R-factor name."""
        with pytest.raises(ValueError, match='Unknown R-factor name'):
            select.select_rfactor_by_name('definitely-not-an-rfactor')
