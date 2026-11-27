"""Tests for the LMAX checks of Rparams.updateDerivedParams.

The i-th phase shift of each block is the one of angular momentum
L = i. These tests check that both the automatic LMAX from
PHASESHIFT_EPS and the warning about large phase shifts beyond a
user-given LMAX use this convention (see #203).
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-24'
__license__ = 'GPLv3+'

from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.classes.rparams.rparams import Rparams
from viperleed.calc.classes.rparams.special.l_max import LMax
from viperleed.calc.lib.leedbase import HARTREE_TO_EV
from viperleed.calc.lib.version import Version

N_VALUES = 21  # Up to L = 20, beyond the supported LMAX of 18
LARGE_PHASESHIFT = 0.5
WARNING = 'large phaseshift values being ignored'


def _phaseshifts_up_to(l_largest):
    """Return phase shifts of two blocks, large up to `l_largest`."""
    site_ps = ([LARGE_PHASESHIFT] * (l_largest + 1)
               + [0.0] * (N_VALUES - l_largest - 1))
    small_ps = [0.0] * N_VALUES
    return [(100 / HARTREE_TO_EV, [small_ps, site_ps])]


@fixture(name='rpars')
def fixture_rpars(tmp_path):
    """Return Rparams for updateDerivedParams, with phase shifts."""
    rpars = Rparams()
    rpars.TENSOR_INDEX = 1
    rpars.paths.tensorleed = tmp_path
    rpars.TL_VERSION = Version('2.0.0')
    rpars.V0_REAL = [1.0, 0.0, 0.0, 0.0]
    rpars.fileLoaded['PHASESHIFTS'] = True
    return rpars


@parametrize(l_largest=(7, 10))
def test_automatic_lmax(rpars, l_largest):
    """Check that LMAX is the largest L with a large phase shift."""
    rpars.phaseshifts = _phaseshifts_up_to(l_largest)
    rpars.updateDerivedParams()
    assert rpars.LMAX.max == l_largest


@parametrize(l_largest=(18, 19, 20))
def test_automatic_lmax_clipped(rpars, l_largest):
    """Check that the automatic LMAX is at most 18."""
    rpars.phaseshifts = _phaseshifts_up_to(l_largest)
    rpars.updateDerivedParams()
    assert rpars.LMAX.max == 18


@parametrize('l_largest,warns', ((8, False), (9, True)),
             ids=('L=LMAX', 'L=LMAX+1'))
def test_warn_large_phaseshift_ignored(rpars, l_largest, warns, caplog):
    """Check the warning for large phase shifts beyond a given LMAX."""
    rpars.LMAX = LMax(6, 8)
    rpars.phaseshifts = _phaseshifts_up_to(l_largest)
    rpars.updateDerivedParams()
    assert (WARNING in caplog.text) is warns
    if warns:
        assert 'L=9 in phase-shift block 2' in caplog.text
        assert rpars.halt >= 1


def test_warn_largest_ignored_value(rpars, caplog):
    """Check that the warning reports the largest ignored value."""
    rpars.LMAX = LMax(6, 8)
    small_ps = [0.0] * N_VALUES
    first_ps, second_ps = small_ps.copy(), small_ps.copy()
    first_ps[9] = 0.3
    second_ps[10] = LARGE_PHASESHIFT
    rpars.phaseshifts = [(100 / HARTREE_TO_EV, [first_ps, second_ps])]
    rpars.updateDerivedParams()
    expect = (f'highest value: {LARGE_PHASESHIFT} for L=10 '
              'in phase-shift block 2')
    assert expect in caplog.text
