"""Cross-check the Python R factor against frozen TensErLEED values."""

__authors__ = ('Florian Kraushofer (@fkraushofer)',)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-23'
__license__ = 'GPLv3+'

from functools import lru_cache
import json

import pytest

from viperleed.calc.classes.rparams.special.energy_range import IVShiftRange
from viperleed.calc.files.beams import readOUTBEAMS

from ....helpers import TEST_DATA
from .conftest import make_rparams
from .conftest import relabel_as_fractional
from .conftest import run_new_rfactor_with_correspondence


DATA_DIR = TEST_DATA / 'rfactor'

# What this cross-check catches: gross deviations of the Python R
# factor from TensErLEED (wrong formula, normalisation, beam averaging,
# V0r scan, grouping). Fine-grained regressions, such as an R that
# moves by less than ~0.002 or a Pendry V0i off by ~10-20 %, stay
# within the tolerances; they are left to the frozen Python values of
# the other R-factor tests. The systems complement each other: Fe2O3
# catches small errors in the Pendry V0i, Ag(100) only large ones.

# Largest accepted difference between the R factor of the Python code
# and that of TensErLEED. The two differ by construction: TensErLEED
# interpolates the theory with local cubic polynomials and takes
# derivatives by finite differences, the Python code uses splines for
# both. TensErLEED also ignores the contribution of Y values below
# 1e-7, and Zanazzi-Jona is evaluated on a grid ten times finer in
# TensErLEED. The tolerances hold for beams spanning several theory
# steps (at least ~30 eV).
R_TOLERANCE = {
    'pendry': 0.005,
    'smooth': 0.005,
    'r2': 0.005,
    'zj': 0.01,
    }

# R of the integer- and fractional-order groups: TensErLEED sets Y to
# zero where both the intensity and its derivative are below 1e-7. On
# weak beams this moves Pendry and smooth R of a group by more than
# the overall R: in the fractional Fe2O3 case, the fractional-order
# Pendry R differs by +0.0064, and by +0.0002 with the cutoff copied
# into the Python Pendry R.
_GROUP_TOLERANCE = {
    **R_TOLERANCE,
    'pendry': 0.01,
    'smooth': 0.01,
    }

# TensErLEED's R2 carries a factor 1/2 that the Python R2 (the usual
# definition) does not have: the Python value is compared with twice
# the TensErLEED one.
R_SCALE = {'r2': 2.0}

# Best V0r (eV) may differ by one step of the shifts
_V0R_TOLERANCE_STEPS = 1
_EPS = 1e-9

# Position of the R of each beam group in Rparams.stored_R
_GROUP_INDEX = {'integer': 1, 'fractional': 2}


@lru_cache(maxsize=None)
def _expected_for(folder_name):
    """Return the frozen reference data of a system."""
    return json.loads((DATA_DIR / folder_name / 'expected.json').read_text())


def _all_cases():
    """Return (system, theory file, R type) of all references."""
    cases = []
    for folder in sorted(DATA_DIR.iterdir()):
        expected = _expected_for(folder.name)
        for theo_name, references in expected['theory'].items():
            cases.extend((folder.name, theo_name, r_type)
                         for r_type in references)
    return cases


_CASES = _all_cases()
_FRACTIONAL_CASES = [case for case in _CASES
                     if _expected_for(case[0])['fractional']]


@lru_cache(maxsize=None)
def _run_python(folder_name, theo_name, r_type):
    """Return the Rparams of a Python R-factor run on a system."""
    folder = DATA_DIR / folder_name
    expected = _expected_for(folder_name)
    beamcorr_file = folder / expected['beam_correspondence_file']
    correspondence = tuple(
        json.loads(beamcorr_file.read_text())['beam_correspondence']
        )
    expbeams = readOUTBEAMS(folder / expected['exp_file'])
    theobeams = readOUTBEAMS(folder / theo_name)
    if expected['fractional']:
        relabel_as_fractional(expbeams)
    rpars = make_rparams({'v0i': expected['v0_imag']},
                          r_type,
                          expbeams,
                          theobeams)
    rpars.IV_SHIFT_RANGE = IVShiftRange(*expected['iv_shift_range'])
    with pytest.MonkeyPatch.context() as patch:
        run_new_rfactor_with_correspondence(patch,
                                            correspondence,
                                            rpars,
                                            theobeams)
    return rpars


@pytest.mark.parametrize(('folder_name', 'theo_name', 'r_type'), _CASES)
def test_settings_match_tensorleed(folder_name, theo_name, r_type):
    """Check the Python run uses the settings of the TensErLEED run."""
    expected = _expected_for(folder_name)
    rpars = _run_python(folder_name, theo_name, r_type)
    assert expected['intpol_deg'] == rpars.INTPOL_DEG
    assert expected['r_factor_smoothing'] == rpars.R_FACTOR_SMOOTHING
    assert expected['theo_energies'] == rpars.THEO_ENERGIES
    assert expected['v0_imag'] == rpars.V0_IMAG
    assert expected['iv_shift_range'] == rpars.IV_SHIFT_RANGE


@pytest.mark.parametrize(('folder_name', 'theo_name', 'r_type'), _CASES)
def test_r_factor_and_v0r(folder_name, theo_name, r_type):
    """Check overall R and best V0r against TensErLEED."""
    expected = _expected_for(folder_name)
    reference = expected['theory'][theo_name][r_type]
    rpars = _run_python(folder_name, theo_name, r_type)
    r_python = rpars.stored_R['refcalc'][0]
    scale = R_SCALE.get(r_type, 1.0)
    assert abs(r_python - scale * reference['r']) <= R_TOLERANCE[r_type]
    step = expected['iv_shift_range'][2]
    assert abs(rpars.best_v0r - reference['best_v0r']) <= (
        _V0R_TOLERANCE_STEPS * step + _EPS
        )


@pytest.mark.parametrize('group', _GROUP_INDEX)
@pytest.mark.parametrize(('folder_name', 'theo_name', 'r_type'),
                         _FRACTIONAL_CASES)
def test_group_r(folder_name, theo_name, r_type, group):
    """Check R of integer or fractional beams against TensErLEED."""
    reference = _expected_for(folder_name)['theory'][theo_name][r_type]
    rpars = _run_python(folder_name, theo_name, r_type)
    r_python = rpars.stored_R['refcalc'][_GROUP_INDEX[group]]
    r_reference = R_SCALE.get(r_type, 1.0) * reference[f'r_{group}']
    assert abs(r_python - r_reference) <= _GROUP_TOLERANCE[r_type]
