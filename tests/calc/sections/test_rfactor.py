"""Tests for the choice of R factor in section Rfactor."""

__authors__ = ('Florian Kraushofer (@fkraushofer)',)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-20'
__license__ = 'GPLv3+'

import pytest
from pytest_cases import parametrize

from viperleed.calc.classes.rparams.rparams import Rparams
from viperleed.calc.classes.rparams.special.energy_range import TheoEnergies
from viperleed.calc.files import parameters
from viperleed.calc.sections import rfactor as rfactor_section

_RUNNERS = {
    'run_legacy_rfactor': 'legacy',
    'run_new_rfactor': 'python',
    }
_LEGACY_CASES = {
    'not given': ('', 'legacy'),
    'true': ('R_FACTOR_LEGACY = True\n', 'legacy'),
    'false': ('R_FACTOR_LEGACY = False\n', 'python'),
    }


def test_r_factor_legacy_default():
    """Check that the TensErLEED R factor is the default."""
    assert Rparams().R_FACTOR_LEGACY is True


@pytest.fixture(name='called_runners')
def fixture_called_runners(monkeypatch):
    """Replace the R-factor runners with ones recording their call."""
    called = []
    for runner, label in _RUNNERS.items():
        def _record(*_, label=label, **__):
            called.append(label)
        monkeypatch.setattr(rfactor_section, runner, _record)
    return called


@parametrize('line,expected', _LEGACY_CASES.values(), ids=_LEGACY_CASES)
def test_rfactor_dispatch(line, expected, called_runners, tmp_path):
    """Check which R factor is run for PARAMETERS with `line`."""
    parameters_file = tmp_path / 'PARAMETERS'
    parameters_file.write_text(f'RUN = 0 1\n{line}', encoding='utf-8')
    rpars = parameters.read(parameters_file)
    parameters.interpret(rpars)
    rpars.THEO_ENERGIES = TheoEnergies.from_sorted_grid((50.0, 60.0, 70.0))
    rpars.refcalc_fdout = 'not empty, no need to fetch spectra'
    rfactor_section.rfactor(None, rpars, index=11)
    assert called_runners == [expected]
