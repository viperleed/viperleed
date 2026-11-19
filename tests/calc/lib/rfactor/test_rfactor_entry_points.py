"""Tests for the entry points of the Python R factor, on real data."""

__authors__ = ('Florian Kraushofer (@fkraushofer)',)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-19'
__license__ = 'GPLv3+'

import json

import numpy as np
import pytest

from viperleed.calc import from_ase
from viperleed.calc.classes.rparams.rparams import Rparams
from viperleed.calc.classes.rparams.special.energy_range import IVShiftRange
from viperleed.calc.classes.rparams.special.energy_range import TheoEnergies
from viperleed.calc.classes.rparams.special.r_factor_type import RFactorType
from viperleed.calc.files.beams import readOUTBEAMS
from viperleed.calc.files.iorfactor import beamlist_to_array
from viperleed.calc.sections import rfactor as rfactor_section

from .conftest import BEAM_DATA_DIR


# Values of R, best V0r and per-beam R at the commit before the
# shared R-factor core was introduced (plan 03). Both entry points
# must keep producing them bit for bit.
_EXPECTED_FILE = (BEAM_DATA_DIR
                  / 'rfactor-entry-points-fe2o3-012-1x1.expected.json')
_SHIFT_RANGE = (-3.0, 3.0)
_STEP = 0.5
_TYPES = tuple(r_type.name.lower() for r_type in RFactorType)


def _run_section(dataset, r_type, monkeypatch):
    """Return the results of run_new_rfactor for `dataset`."""
    monkeypatch.setattr(rfactor_section.leedbase,
                        'getBeamCorrespondence',
                        lambda *_: dataset['beam_correspondence'])
    expbeams = readOUTBEAMS(dataset['exp_file'])
    theobeams = readOUTBEAMS(dataset['theo_file'])
    rpars = Rparams()
    rpars.R_FACTOR_TYPE = RFactorType[r_type.upper()]
    rpars.V0_IMAG = dataset['v0i']
    rpars.PLOT_IV['plot'] = False
    rpars.expbeams = expbeams
    theo_grid, *_ = beamlist_to_array(theobeams)
    rpars.THEO_ENERGIES = TheoEnergies.from_sorted_grid(theo_grid)
    rpars.IV_SHIFT_RANGE = IVShiftRange(*_SHIFT_RANGE, _STEP)
    r_beams = rfactor_section.run_new_rfactor(None,
                                              rpars,
                                              for_error=False,
                                              name='refcalc',
                                              theobeams=theobeams,
                                              expbeams=expbeams)
    overall, integer, fractional = rpars.stored_R['refcalc']
    return {
        'best_v0r': float(rpars.best_v0r),
        'r_overall': float(overall),
        'r_integer': float(integer),
        'r_fractional': float(fractional),
        'r_beams': np.asarray(r_beams, dtype=float).tolist(),
        }


def _run_from_csv(dataset, r_type):
    """Return the results of rfactor_from_csv for `dataset`."""
    best_r, best_v0r, r_beams = from_ase.rfactor_from_csv(
        (dataset['exp_file'], dataset['theo_file']),
        dataset['v0i'],
        v0r_shift_range=_SHIFT_RANGE,
        intpol_step=_STEP,
        which_r=r_type,
        beam_correspondence=dataset['beam_correspondence'],
        )
    return {
        'best_v0r': float(best_v0r),
        'r_overall': float(best_r),
        'r_beams': np.asarray(r_beams, dtype=float).tolist(),
        }


def _load_expected(entry_point, r_type):
    """Return the frozen values for `entry_point` and `r_type`."""
    return json.loads(_EXPECTED_FILE.read_text())[entry_point][r_type]


@pytest.mark.parametrize('r_type', _TYPES)
def test_run_new_rfactor(r_type, fe2o3_012_1x1_dataset, monkeypatch):
    """Check run_new_rfactor bit for bit against frozen values."""
    actual = _run_section(fe2o3_012_1x1_dataset, r_type, monkeypatch)
    # assert_equal is exact, and NaN (no fractional beams) equals NaN
    np.testing.assert_equal(actual,
                            _load_expected('run_new_rfactor', r_type))


@pytest.mark.parametrize('r_type', _TYPES)
def test_rfactor_from_csv(r_type, fe2o3_012_1x1_dataset):
    """Check rfactor_from_csv bit for bit against frozen values."""
    actual = _run_from_csv(fe2o3_012_1x1_dataset, r_type)
    np.testing.assert_equal(actual,
                            _load_expected('rfactor_from_csv', r_type))
