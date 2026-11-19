"""Tests for the entry points of the Python R factor, on real data."""

__authors__ = ('Florian Kraushofer (@fkraushofer)',)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-19'
__license__ = 'GPLv3+'

from fractions import Fraction
import json

import numpy as np
import pytest

from viperleed.calc import from_ase
from viperleed.calc.classes.rparams.rparams import Rparams
from viperleed.calc.classes.rparams.special.energy_range import IVShiftRange
from viperleed.calc.classes.rparams.special.energy_range import TheoEnergies
from viperleed.calc.classes.rparams.special.r_factor_type import RFactorType
from viperleed.calc.files.beams import readOUTBEAMS
from viperleed.calc.files import iorfactor
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

# Shift ranges (start, stop) and steps for both entry points. The
# frozen values also cover the shift at the edge of the range, a fixed
# shift (best shift of 0 or not), and a range that is not symmetric.
_RANGES = {
    'best-at-bound': ((-3, 0.5), 0.5),
    'fixed-zero': ((0, 0), 0.5),
    'fixed-nonzero': ((1, 1), 0.5),
    'asymmetric': ((-1, 6), 0.5),
    'step-quarter': ((-2, 2), 0.25),
    'step-one': ((-3, 3), 1.0),
    # Shifts that are not exactly representable in binary floating point
    'step-tenths': ((-1.5, 1.5), 0.3),
    }
# Cases only reachable via run_new_rfactor: an undefined step, a
# second calculation for error curves after a normal one, and
# fractional-order experimental beams.
_SECTION_ONLY = {
    'undefined-step': {'step': None},
    'for-error-step-one': {'step': 1.0, 'for_error': True},
    'for-error-undefined-step': {'step': None, 'for_error': True},
    'fractional-beams': {'fractional': True},
    }
_PDF_ENTRY_POINTS = ('run_new_rfactor', 'rfactor_from_csv')


def _make_rparams(dataset, r_type, expbeams, theobeams):
    """Return an Rparams ready for run_new_rfactor."""
    rpars = Rparams()
    rpars.R_FACTOR_TYPE = RFactorType[r_type.upper()]
    rpars.V0_IMAG = dataset['v0i']
    rpars.PLOT_IV['plot'] = False
    rpars.expbeams = expbeams
    theo_grid, *_ = beamlist_to_array(theobeams)
    rpars.THEO_ENERGIES = TheoEnergies.from_sorted_grid(theo_grid)
    return rpars


def _relabel_as_fractional(beams):
    """Shift the h index of every third beam by a half."""
    for beam in beams[::3]:
        beam.hk = (beam.hk[0] + Fraction(1, 2), beam.hk[1])


def _run_section(dataset, r_type, monkeypatch, **options):
    """Return the results of run_new_rfactor for `dataset`.

    Parameters
    ----------
    dataset : dict
        The beam data.
    r_type : str
        Name of the R factor.
    monkeypatch : pytest.MonkeyPatch
        Used to fix the beam correspondence.
    **options
        shift_range : tuple of two float
            Start and stop of the V0r shifts. Default is _SHIFT_RANGE.
        step : float or None
            Step of the shifts. Default is _STEP.
        for_error : bool
            Whether the results are those of a calculation for error
            curves, run after a normal one. Default is False.
        fractional : bool
            Whether some experimental beams are fractional-order.
            Default is False.
        plot : bool
            Whether the PDF writer is called. Default is False.

    Returns
    -------
    dict
        Results of the calculation.
    """
    shift_range = options.get('shift_range', _SHIFT_RANGE)
    step = options.get('step', _STEP)
    for_error = options.get('for_error', False)
    monkeypatch.setattr(rfactor_section.leedbase,
                        'getBeamCorrespondence',
                        lambda *_: dataset['beam_correspondence'])
    expbeams = readOUTBEAMS(dataset['exp_file'])
    theobeams = readOUTBEAMS(dataset['theo_file'])
    if options.get('fractional', False):
        _relabel_as_fractional(expbeams)
    rpars = _make_rparams(dataset, r_type, expbeams, theobeams)
    rpars.PLOT_IV['plot'] = options.get('plot', False)
    shift = shift_range if step is None else (*shift_range, step)
    rpars.IV_SHIFT_RANGE = IVShiftRange(*shift)
    name = 'refcalc'
    if for_error:
        rfactor_section.run_new_rfactor(None,
                                        rpars,
                                        for_error=False,
                                        name=name,
                                        theobeams=theobeams,
                                        expbeams=expbeams)
        name = 'superpos'
    r_beams = rfactor_section.run_new_rfactor(None,
                                              rpars,
                                              for_error=for_error,
                                              name=name,
                                              theobeams=theobeams,
                                              expbeams=expbeams)
    overall, integer, r_fractional = rpars.stored_R[name]
    return {
        'best_v0r': float(rpars.best_v0r),
        'r_overall': float(overall),
        'r_integer': float(integer),
        'r_fractional': float(r_fractional),
        'r_beams': np.asarray(r_beams, dtype=float).tolist(),
        }


def _run_from_csv(dataset, r_type,
                  shift_range=_SHIFT_RANGE, step=_STEP, plot_iv=None):
    """Return the results of rfactor_from_csv for `dataset`."""
    best_r, best_v0r, r_beams = from_ase.rfactor_from_csv(
        (dataset['exp_file'], dataset['theo_file']),
        dataset['v0i'],
        v0r_shift_range=shift_range,
        intpol_step=step,
        which_r=r_type,
        beam_correspondence=dataset['beam_correspondence'],
        plot_iv=plot_iv,
        )
    return {
        'best_v0r': float(best_v0r),
        'r_overall': float(best_r),
        'r_beams': np.asarray(r_beams, dtype=float).tolist(),
        }


def _sum_per_beam(intensities):
    """Return the sum over energies of each beam (a compact digest)."""
    return np.nansum(intensities, axis=0).tolist()


def _capture_pdf_writer(monkeypatch):
    """Replace the PDF writer and return the list of its calls."""
    calls = []

    def _fake_writer(rpars, theo_spline, exp_spline, out_grid, r_beams,
                     **kwargs):
        calls.append({
            'exp_labels': [beam.label for beam in rpars.expbeams],
            'v0_imag': rpars.V0_IMAG,
            'plot_iv': dict(rpars.PLOT_IV),
            'r_factor_type': rpars.R_FACTOR_TYPE.name,
            'out_grid_ends': [out_grid[0], out_grid[-1]],
            'n_out_grid': len(out_grid),
            'theo_sums': _sum_per_beam(theo_spline(out_grid)),
            'exp_sums': _sum_per_beam(exp_spline(out_grid)),
            'r_beams': np.asarray(r_beams).tolist(),
            'kwargs': kwargs,
            })

    # The writer may be bound by name in from_ase or looked up in
    # iorfactor: replace it in both places.
    monkeypatch.setattr(iorfactor,
                        'write_Rfactorpdf_from_splines',
                        _fake_writer)
    monkeypatch.setattr(from_ase,
                        'write_Rfactorpdf_from_splines',
                        _fake_writer,
                        raising=False)
    return calls


def _run_pdf(entry_point, dataset, monkeypatch):
    """Return the only call to the PDF writer made by `entry_point`."""
    calls = _capture_pdf_writer(monkeypatch)
    if entry_point == 'run_new_rfactor':
        _run_section(dataset, 'pendry', monkeypatch, plot=True)
    else:
        _run_from_csv(dataset, 'pendry', plot_iv={'plot': True})
    assert len(calls) == 1
    return calls[0]


def _load_expected(entry_point, r_type, scenario=None):
    """Return the frozen values for `entry_point` and `r_type`."""
    expected = json.loads(_EXPECTED_FILE.read_text())
    if scenario is None:
        return expected[entry_point][r_type]
    return expected['scenarios'][scenario][entry_point][r_type]


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


@pytest.mark.parametrize('r_type', _TYPES)
@pytest.mark.parametrize('scenario', _RANGES)
def test_run_new_rfactor_range(scenario, r_type,
                               fe2o3_012_1x1_dataset, monkeypatch):
    """Check run_new_rfactor for other shift ranges and steps."""
    shift_range, step = _RANGES[scenario]
    actual = _run_section(fe2o3_012_1x1_dataset, r_type, monkeypatch,
                          shift_range=shift_range,
                          step=step)
    np.testing.assert_equal(
        actual,
        _load_expected('run_new_rfactor', r_type, scenario),
        )


@pytest.mark.parametrize('r_type', _TYPES)
@pytest.mark.parametrize('scenario', _RANGES)
def test_rfactor_from_csv_range(scenario, r_type, fe2o3_012_1x1_dataset):
    """Check rfactor_from_csv for other shift ranges and steps."""
    shift_range, step = _RANGES[scenario]
    actual = _run_from_csv(fe2o3_012_1x1_dataset,
                           r_type,
                           shift_range=shift_range,
                           step=step)
    np.testing.assert_equal(
        actual,
        _load_expected('rfactor_from_csv', r_type, scenario),
        )


@pytest.mark.parametrize('r_type', _TYPES)
@pytest.mark.parametrize('scenario', _SECTION_ONLY)
def test_run_new_rfactor_special(scenario, r_type,
                                 fe2o3_012_1x1_dataset, monkeypatch):
    """Check run_new_rfactor for cases without a csv counterpart."""
    kwargs = {'step': _STEP, **_SECTION_ONLY[scenario]}
    actual = _run_section(fe2o3_012_1x1_dataset, r_type, monkeypatch,
                          **kwargs)
    np.testing.assert_equal(
        actual,
        _load_expected('run_new_rfactor', r_type, scenario),
        )


@pytest.mark.parametrize('entry_point', _PDF_ENTRY_POINTS)
def test_pdf_inputs(entry_point, fe2o3_012_1x1_dataset, monkeypatch):
    """Check what the entry points hand to the PDF writer."""
    actual = _run_pdf(entry_point, fe2o3_012_1x1_dataset, monkeypatch)
    expected = _load_expected(entry_point, 'pendry', 'pdf')
    np.testing.assert_equal(actual, expected)
