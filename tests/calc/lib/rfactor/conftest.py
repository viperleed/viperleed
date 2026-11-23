"""Shared fixtures and helpers for rfactor tests."""

__authors__ = ('Alexandra Mia Imre (@alexmiame)',
               'Florian Kraushofer (@fkraushofer)')
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2024-02-27'
__license__ = 'GPLv3+'

from fractions import Fraction
import json

import numpy as np
import pytest
from scipy.interpolate import CubicSpline

from viperleed.calc.classes.rparams.rparams import Rparams
from viperleed.calc.classes.rparams.special.energy_range import TheoEnergies
from viperleed.calc.classes.rparams.special.r_factor_type import RFactorType
from viperleed.calc.files.beams import readOUTBEAMS
from viperleed.calc.files.iorfactor import beamlist_to_array
from viperleed.calc.lib.rfactor import pendry
from viperleed.calc.lib.rfactor import r_1
from viperleed.calc.lib.rfactor import r_2
from viperleed.calc.lib.rfactor import smooth
from viperleed.calc.lib.rfactor import zannazi_jona
from viperleed.calc.lib.rfactor.utils import average_beam_array
from viperleed.calc.lib.spline_interpolation import interpolate_ragged_array
from viperleed.calc.sections import rfactor as rfactor_section

from ....helpers import TEST_DATA


_MOCK_ROUGH_E_AXIS = np.arange(15.0, 300.0, 3.0)
_MOCK_FINE_E_STEP = 0.5
_MOCK_FINE_E_AXIS = np.arange(15.0, 300.0, _MOCK_FINE_E_STEP)

_MOCK_SIN_INTENSITY = np.sin(_MOCK_ROUGH_E_AXIS / (3.0 * np.pi)) ** 2
_MOCK_COS_INTENSITY = np.cos(_MOCK_ROUGH_E_AXIS / (3.0 * np.pi)) ** 2
_MOCK_CONST_INTENSITY = np.ones_like(_MOCK_ROUGH_E_AXIS)

FE2O3_012_DIR = TEST_DATA / 'rfactor' / 'fe2o3-012'
EXP_BEAMS_FILE = FE2O3_012_DIR / 'EXPBEAMS.csv'
THEO_BEAMS_FILE = FE2O3_012_DIR / 'THEOBEAMS.csv'
BEAMCORR_FILE = FE2O3_012_DIR / 'beamcorr.json'
EXPECTED_FILE = FE2O3_012_DIR / 'rfactor-regression.expected.json'
ENTRY_POINTS_EXPECTED_FILE = (FE2O3_012_DIR
                              / 'rfactor-entry-points.expected.json')


def pytest_addoption(parser):
    parser.addoption(
        '--update-rfactor-expected',
        action='store_true',
        default=False,
        help=(
            'Regenerate rfactor-regression.expected.json of the '
            'Fe2O3(012) data. The frozen values of the entry-point '
            'tests (rfactor-entry-points.expected.json) are not '
            'regenerated.'
        ),
    )


@pytest.fixture
def sin_spline():
    return CubicSpline(_MOCK_ROUGH_E_AXIS, _MOCK_SIN_INTENSITY)


@pytest.fixture
def cos_spline():
    return CubicSpline(_MOCK_ROUGH_E_AXIS, _MOCK_COS_INTENSITY)


@pytest.fixture
def const_spline():
    return CubicSpline(_MOCK_ROUGH_E_AXIS, _MOCK_CONST_INTENSITY)


@pytest.fixture
def synthetic_pair(sin_spline, cos_spline):
    return {
        'v0i': 3.0,
        'energy_step': _MOCK_FINE_E_STEP,
        'energies': _MOCK_FINE_E_AXIS,
        'data_spline_1': sin_spline,
        'data_spline_2': cos_spline,
    }


def evaluate_spline_triplet(data_spline, energy_grid):
    """Return a spline and its first two derivatives on energy_grid."""
    deriv_1 = data_spline.derivative()
    deriv_2 = deriv_1.derivative()
    return (
        data_spline(energy_grid),
        deriv_1(energy_grid),
        deriv_2(energy_grid),
    )


def _jsonify_hk(hk_values):
    return [list(hk) for hk in hk_values]


def _load_beam_correspondence():
    data = json.loads(BEAMCORR_FILE.read_text())
    return tuple(data['beam_correspondence'])


@pytest.fixture(scope='session')
def fe2o3_012_1x1_dataset():
    """Real beam data with fixed beam correspondence applied to THEO beams."""
    exp_out_beams = readOUTBEAMS(EXP_BEAMS_FILE)
    theo_out_beams = readOUTBEAMS(THEO_BEAMS_FILE)

    exp_grid, _, _, exp_beams = beamlist_to_array(exp_out_beams)
    theo_grid, _, _, theo_beams = beamlist_to_array(theo_out_beams)

    exp_hk = [beam.hk for beam in exp_out_beams]
    theo_hk = [beam.hk for beam in theo_out_beams]

    beam_corr = _load_beam_correspondence()

    if len(beam_corr) != len(exp_hk):
        raise ValueError(
            f'Beam correspondence length {len(beam_corr)} does not match '
            f'number of experimental beams {len(exp_hk)}'
        )
    if len(beam_corr) != len(theo_hk):
        raise ValueError(
            f'Beam correspondence length {len(beam_corr)} does not match '
            f'number of theoretical beams {len(theo_hk)}'
        )

    theo_beams_effective = average_beam_array(theo_beams, beam_corr)
    assert theo_beams_effective.shape[1] == len(exp_hk)

    exp_spline = interpolate_ragged_array(exp_grid, exp_beams)
    theo_spline = interpolate_ragged_array(theo_grid, theo_beams_effective)

    emin = max(np.nanmin(exp_grid), np.nanmin(theo_grid))
    emax = min(np.nanmax(exp_grid), np.nanmax(theo_grid))
    energy_step = 0.5
    energies = np.arange(emin, emax + 0.25 * energy_step, energy_step)

    exp_triplet = evaluate_spline_triplet(exp_spline, energies)
    theo_triplet = evaluate_spline_triplet(theo_spline, energies)

    return {
        'label': 'fe2o3-012-1x1',
        'exp_file': EXP_BEAMS_FILE,
        'theo_file': THEO_BEAMS_FILE,
        'expected_entry_points': ENTRY_POINTS_EXPECTED_FILE,
        'v0i': 3.0,
        'energy_step': energy_step,
        'energies': energies,
        'hk': exp_hk,
        'beam_correspondence': beam_corr,
        'exp_spline': exp_spline,
        'theo_spline': theo_spline,
        'theo_grid': theo_grid,
        'theo_beams': theo_beams_effective,
        'exp_triplet': exp_triplet,
        'theo_triplet': theo_triplet,
    }


@pytest.fixture(scope='session')
def regression_payload(fe2o3_012_1x1_dataset):
    """Return the R factors of the Fe2O3(012) data, as stored."""
    dataset = fe2o3_012_1x1_dataset
    funcs = {
        'pendry': pendry.r_pendry,
        'r1': r_1.r_1,
        'r2': r_2.r_2,
        'rs': smooth.r_s,
        'zj': zannazi_jona.r_zjj,
    }

    payload = {
        'label': dataset['label'],
        'v0i': dataset['v0i'],
        'energy_step': dataset['energy_step'],
        'n_beams': len(dataset['hk']),
        'hk': _jsonify_hk(dataset['hk']),
        'beam_correspondence': list(dataset['beam_correspondence']),
    }

    base_kwargs_spline = {
        'v0_imag': dataset['v0i'],
        'energy_step': dataset['energy_step'],
        'energy_grid': dataset['energies'],
        'data_spline_1': dataset['exp_spline'],
        'data_spline_2': dataset['theo_spline'],
    }
    base_kwargs_arrays = {
        'v0_imag': dataset['v0i'],
        'energy_step': dataset['energy_step'],
        'energy_grid': dataset['energies'],
        'data_and_derivatives_1': dataset['exp_triplet'],
        'data_and_derivatives_2': dataset['theo_triplet'],
    }

    for name, func in funcs.items():
        overall = func(**base_kwargs_arrays)
        per_beam = func(**base_kwargs_arrays, groups='beam')
        overall_from_splines = func(**base_kwargs_spline)
        per_beam_from_splines = func(**base_kwargs_spline, groups='beam')

        payload[name] = {
            'overall': float(np.asarray(overall)),
            'per_beam': np.asarray(per_beam, dtype=float).tolist(),
            'overall_from_splines': float(np.asarray(overall_from_splines)),
            'per_beam_from_splines': np.asarray(
                per_beam_from_splines, dtype=float
            ).tolist(),
        }

    return payload


def make_rparams(dataset, r_type, expbeams, theobeams):
    """Return an Rparams ready for run_new_rfactor.

    Parameters
    ----------
    dataset : dict
        The beam data. Only the 'v0i' entry is used.
    r_type : str
        Name of the R factor.
    expbeams : Sequence of Beam
        The experimental beams.
    theobeams : Sequence of Beam
        The theoretical beams.

    Returns
    -------
    Rparams
        Parameters with the R-factor type, V0_IMAG, experimental beams
        and theory energies set.
    """
    rpars = Rparams()
    rpars.R_FACTOR_TYPE = RFactorType[r_type.upper()]
    rpars.V0_IMAG = dataset['v0i']
    rpars.PLOT_IV['plot'] = False
    rpars.expbeams = expbeams
    theo_grid, *_ = beamlist_to_array(theobeams)
    rpars.THEO_ENERGIES = TheoEnergies.from_sorted_grid(theo_grid)
    return rpars


def relabel_as_fractional(beams):
    """Shift the h index of every third beam by a half."""
    for beam in beams[::3]:
        beam.hk = (beam.hk[0] + Fraction(1, 2), beam.hk[1])


def run_new_rfactor_with_correspondence(monkeypatch,
                                        correspondence,
                                        rpars,
                                        theobeams,
                                        **kwargs):
    """Call run_new_rfactor with a fixed beam correspondence.

    Parameters
    ----------
    monkeypatch : pytest.MonkeyPatch
        Used to replace getBeamCorrespondence.
    correspondence : Sequence
        The beam correspondence to use.
    rpars : Rparams
        Parameters of the calculation.
    theobeams : Sequence of Beam
        The theoretical beams.
    **kwargs
        Passed on to run_new_rfactor. 'for_error' defaults to False
        and 'name' to 'refcalc'.

    Returns
    -------
    numpy.ndarray
        The R factor of each experimental beam.
    """
    monkeypatch.setattr(rfactor_section.leedbase,
                        'getBeamCorrespondence',
                        lambda *_: correspondence)
    kwargs.setdefault('for_error', False)
    kwargs.setdefault('name', 'refcalc')
    return rfactor_section.run_new_rfactor(None,
                                           rpars,
                                           theobeams=theobeams,
                                           **kwargs)
