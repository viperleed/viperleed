"""Tests for viperleed.calc section refcalc.

This module contains tests for the selection of the LMAX values
of a reference calculation (_select_lmax).
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-16'
__license__ = 'GPLv3+'

import logging

from pytest_cases import fixture
from pytest_cases import parametrize

from viperleed.calc.classes.rparams.rparams import Rparams
from viperleed.calc.classes.rparams.special.l_max import LMax
from viperleed.calc.lib.leedbase import HARTREE_TO_EV
from viperleed.calc.lib.version import Version
from viperleed.calc.sections.refcalc import _select_lmax

ENERGIES = 50.0, 150.0, 250.0
PHASESHIFT_EPS = 0.01


def _phaseshifts(*n_nonzero):
    """Return phase shifts of one site with `n_nonzero` large values.

    The energies of the phase shifts are slightly above those of
    ENERGIES, so that the i-th energy uses the i-th phase shifts.
    """
    n_values = 19
    return [
        ((energy + 1) / HARTREE_TO_EV,
         [[0.1 * (i + 1) for i in range(nonzero)]
          + [0.0] * (n_values - nonzero)])
        for energy, nonzero in zip(ENERGIES, n_nonzero)
        ]


@fixture(name='rpars')
def fixture_rpars():
    """Return Rparams for a parallel calculation with an LMAX range."""
    rpars = Rparams()
    rpars.N_CORES = 4
    rpars.LMAX = LMax(6, 12)
    rpars.TL_VERSION = Version('1.7.3')
    rpars.PHASESHIFT_EPS = PHASESHIFT_EPS
    # Largest L above PHASESHIFT_EPS: 3, 8, 15
    rpars.phaseshifts = _phaseshifts(4, 9, 16)
    return rpars


_ONE_LMAX = {
    'one core': {'N_CORES': 1},
    'single LMAX': {'LMAX': LMax(10)},
    'TensErLEED 1.6': {'TL_VERSION': Version('1.6.0')},
    }


@parametrize('attrs', _ONE_LMAX.values(), ids=_ONE_LMAX)
def test_one_lmax(rpars, attrs):
    """Check that the largest LMAX is used for all energies."""
    for attr, value in attrs.items():
        setattr(rpars, attr, value)
    lmax = _select_lmax(rpars, ENERGIES)
    assert lmax == dict.fromkeys(ENERGIES, rpars.LMAX.max)


def test_lmax_per_energy(rpars):
    """Check that LMAX is chosen per energy within the LMAX range."""
    lmax = _select_lmax(rpars, ENERGIES)
    assert lmax == dict(zip(ENERGIES, (6, 9, 12)))


def test_lmax_per_energy_supported_max(rpars, mocker, caplog):
    """Check that LMAX per energy is at most the supported maximum."""
    rpars.LMAX = LMax(6, 18)
    mocker.patch.object(rpars, 'get_limits', return_value=(1, 10))
    with caplog.at_level(logging.INFO):
        lmax = _select_lmax(rpars, ENERGIES)
    assert lmax == dict(zip(ENERGIES, (6, 9, 10)))
    expect_log = 'LMAX was set to 10'
    assert expect_log in caplog.text


def test_lmax_per_energy_all_small(rpars):
    """Check that LMAX.min is used if no phase shift is large."""
    rpars.phaseshifts = _phaseshifts(0, 0, 0)
    lmax = _select_lmax(rpars, ENERGIES)
    assert lmax == dict.fromkeys(ENERGIES, rpars.LMAX.min)


def test_lmax_per_energy_repeated_value(rpars):
    """Check that LMAX comes from the largest L of repeated values."""
    n_values = len(rpars.phaseshifts[0][1][0])
    n_large = 8
    site_ps = [0.5] * n_large + [0.0] * (n_values - n_large)
    rpars.phaseshifts = [(energy, [site_ps])
                         for energy, _ in rpars.phaseshifts]
    lmax = _select_lmax(rpars, ENERGIES)
    assert lmax == dict.fromkeys(ENERGIES, n_large)


_ABOVE_PHASESHIFTS = {  # Highest phase-shift energy is ENERGIES[-1] + 1
    'within 1 eV': (ENERGIES[-1] + 1.5, False),
    'more than 1 eV': (ENERGIES[-1] + 3.0, True),
    }


@parametrize('energy,warns', _ABOVE_PHASESHIFTS.values(),
             ids=_ABOVE_PHASESHIFTS)
def test_lmax_above_phaseshifts(rpars, energy, warns, caplog):
    """Check the warning for energies above those of the phase shifts."""
    lmax = _select_lmax(rpars, (energy,))
    assert lmax == {energy: rpars.LMAX.max}  # From the last phase shifts
    expect_log = 'No appropriate phaseshifts found'
    assert (expect_log in caplog.text) is warns
