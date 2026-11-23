"""Tests for module selection of viperleed.calc.lib.rfactor."""

__authors__ = ('Alexandra Mia Imre (@alexmiame)',
               'Florian Kraushofer (@fkraushofer)')
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-10-08'
__license__ = 'GPLv3+'

import pytest
from pytest_cases import parametrize

from viperleed.calc.classes.rparams.special.r_factor_type import _ALIASES
from viperleed.calc.classes.rparams.special.r_factor_type import RFactorType
from viperleed.calc.lib.rfactor import pendry
from viperleed.calc.lib.rfactor import r_1
from viperleed.calc.lib.rfactor import r_2
from viperleed.calc.lib.rfactor import selection
from viperleed.calc.lib.rfactor import smooth
from viperleed.calc.lib.rfactor import zannazi_jona

_FUNCTIONS = {
    RFactorType.PENDRY: pendry.r_pendry,
    RFactorType.R2: r_2.r_2,
    RFactorType.ZJ: zannazi_jona.r_zjj,
    RFactorType.SMOOTH: smooth.r_s,
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
_CODES = ('1', '2', '3', '4')


def _cases(aliases):
    """Yield (name, key) pairs from a {key: names} dictionary."""
    for key, names in aliases.items():
        for name in names:
            yield name, key


# That all the R_FACTOR_TYPE aliases (_ALIASES) parse correctly is
# checked in test_r_factor_type. Here we only check that both selection
# functions accept all of them.
@parametrize(select_func=(selection.select_rfactor,
                          selection.select_rfactor_by_name))
@parametrize('name,rfactor_type', _cases(_ALIASES))
def test_alias(select_func, name, rfactor_type):
    """Check the function selected for an R_FACTOR_TYPE alias."""
    assert select_func(name) is _FUNCTIONS[rfactor_type]


class TestSelectRfactor:
    """Tests for the select_rfactor function."""

    @parametrize(rfactor_type=RFactorType)
    def test_member(self, rfactor_type):
        """Check the function selected for an RFactorType."""
        selected = selection.select_rfactor(rfactor_type)
        assert selected is _FUNCTIONS[rfactor_type]

    @parametrize(code=_CODES)
    def test_code(self, code):
        """Check the function selected for a TensErLEED integer code."""
        selected = selection.select_rfactor(code)
        assert selected is _FUNCTIONS[RFactorType(int(code))]

    @parametrize(name=_R1_NAMES)
    def test_no_r1(self, name):
        """Check that the R1 R factor cannot be selected."""
        with pytest.raises(ValueError, match='not a valid RFactorType'):
            selection.select_rfactor(name)


class TestSelectRfactorByName:
    """Tests for the select_rfactor_by_name function."""

    @parametrize('name,func', _cases(_FROM_ASE_NAMES))
    def test_from_ase_names(self, name, func):
        """Check that the names from_ase accepted before still work."""
        assert selection.select_rfactor_by_name(name) is func

    @parametrize('name,expect', (('  RP  ', pendry.r_pendry),
                                 (' R1 ', r_1.r_1),
                                 ('R_1', r_1.r_1)))
    def test_strips_and_lowercases(self, name, expect):
        """Check that names are case-insensitive and stripped."""
        assert selection.select_rfactor_by_name(name) is expect

    @parametrize(name=(*_CODES, ' 1 '))
    def test_rejects_codes(self, name):
        """Check that integer codes are not taken as R-factor names."""
        with pytest.raises(ValueError, match='Unknown R-factor name'):
            selection.select_rfactor_by_name(name)

    def test_unknown(self):
        """Check complaints for an unknown R-factor name."""
        with pytest.raises(ValueError, match='Unknown R-factor name'):
            selection.select_rfactor_by_name('definitely-not-an-rfactor')
