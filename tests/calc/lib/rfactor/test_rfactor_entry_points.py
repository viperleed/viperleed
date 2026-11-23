"""Tests for the entry points of the Python R factor, on real data."""

__authors__ = ('Florian Kraushofer (@fkraushofer)',)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-19'
__license__ = 'GPLv3+'

import json
import logging

import numpy as np
import pytest

from viperleed.calc import from_ase
from viperleed.calc.classes.rparams.special.energy_range import IVShiftRange
from viperleed.calc.classes.rparams.special.r_factor_type import RFactorType
from viperleed.calc.files import iorfactor
from viperleed.calc.files.beams import readOUTBEAMS
from viperleed.calc.files.iorfactor import RfactorError
from viperleed.calc.files.parameters.errors import ParameterConflictError
from viperleed.calc.lib.matplotlib_utils import CAN_PLOT
from viperleed.calc.sections import rfactor as rfactor_section

from .conftest import make_rparams
from .conftest import relabel_as_fractional
from .conftest import run_new_rfactor_with_correspondence


# Values of R, best V0r and per-beam R recorded before the R-factor
# core was shared between the entry points. Both entry points must
# keep producing them bit for bit.
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


# Energy offset (eV) applied to the experimental beams in a calculation
# for error curves, so the best V0r would differ from the fixed one.
_FOR_ERROR_DATA_SHIFT = 2.0


def _shift_beams(beams, delta):
    """Add `delta` (eV) to all energies of `beams`."""
    for beam in beams:
        beam.intens = {energy + delta: intensity
                       for energy, intensity in beam.intens.items()}


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
    expbeams = readOUTBEAMS(dataset['exp_file'])
    theobeams = readOUTBEAMS(dataset['theo_file'])
    if options.get('fractional', False):
        relabel_as_fractional(expbeams)
    rpars = make_rparams(dataset, r_type, expbeams, theobeams)
    rpars.PLOT_IV['plot'] = options.get('plot', False)
    shift = shift_range if step is None else (*shift_range, step)
    rpars.IV_SHIFT_RANGE = IVShiftRange(*shift)
    name = 'refcalc'
    normal_pass = None
    if for_error:
        run_new_rfactor_with_correspondence(
            monkeypatch,
            dataset['beam_correspondence'],
            rpars,
            theobeams,
            name=name,
            )
        normal_pass = rpars.best_v0r
        _shift_beams(rpars.expbeams, _FOR_ERROR_DATA_SHIFT)
        name = 'superpos'
    r_beams = run_new_rfactor_with_correspondence(
        monkeypatch,
        dataset['beam_correspondence'],
        rpars,
        theobeams,
        for_error=for_error,
        name=name,
        )
    overall, integer, r_fractional = rpars.stored_R[name]
    if normal_pass is not None:
        # The data changed such that the optimum would move. Instead,
        # R is evaluated at the fixed V0r of the normal pass.
        assert rpars.best_v0r == normal_pass
        reference = make_rparams(dataset, r_type, expbeams, theobeams)
        reference.IV_SHIFT_RANGE = IVShiftRange(normal_pass,
                                                normal_pass,
                                                rpars.IV_SHIFT_RANGE.step)
        run_new_rfactor_with_correspondence(
            monkeypatch,
            dataset['beam_correspondence'],
            reference,
            theobeams,
            name=name,
            )
        np.testing.assert_allclose(rpars.stored_R[name],
                                   reference.stored_R[name],
                                   rtol=1e-12)
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


def _load_expected(dataset, entry_point, r_type, scenario=None):
    """Return the frozen values for `entry_point` and `r_type`."""
    expected = json.loads(dataset['expected_entry_points'].read_text())
    if scenario is None:
        return expected[entry_point][r_type]
    return expected['scenarios'][scenario][entry_point][r_type]


@pytest.mark.parametrize('r_type', _TYPES)
def test_run_new_rfactor(r_type, fe2o3_012_1x1_dataset, monkeypatch):
    """Check run_new_rfactor bit for bit against frozen values."""
    actual = _run_section(fe2o3_012_1x1_dataset, r_type, monkeypatch)
    # assert_equal is exact, and NaN (no fractional beams) equals NaN
    np.testing.assert_equal(actual,
                            _load_expected(fe2o3_012_1x1_dataset,
                                           'run_new_rfactor',
                                           r_type))


@pytest.mark.parametrize('r_type', _TYPES)
def test_rfactor_from_csv(r_type, fe2o3_012_1x1_dataset):
    """Check rfactor_from_csv bit for bit against frozen values."""
    actual = _run_from_csv(fe2o3_012_1x1_dataset, r_type)
    np.testing.assert_equal(actual,
                            _load_expected(fe2o3_012_1x1_dataset,
                                           'rfactor_from_csv',
                                           r_type))


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
        _load_expected(fe2o3_012_1x1_dataset,
                       'run_new_rfactor',
                       r_type,
                       scenario),
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
        _load_expected(fe2o3_012_1x1_dataset,
                       'rfactor_from_csv',
                       r_type,
                       scenario),
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
        _load_expected(fe2o3_012_1x1_dataset,
                       'run_new_rfactor',
                       r_type,
                       scenario),
        )


@pytest.mark.parametrize('entry_point', _PDF_ENTRY_POINTS)
def test_pdf_inputs(entry_point, fe2o3_012_1x1_dataset, monkeypatch):
    """Check what the entry points hand to the PDF writer."""
    actual = _run_pdf(entry_point, fe2o3_012_1x1_dataset, monkeypatch)
    expected = _load_expected(fe2o3_012_1x1_dataset,
                               entry_point,
                               'pendry',
                               'pdf')
    np.testing.assert_equal(actual, expected)


def test_run_new_rfactor_rejects_intpol_deg(fe2o3_012_1x1_dataset):
    """Check that INTPOL_DEG != 3 set in code is refused."""
    expbeams = readOUTBEAMS(fe2o3_012_1x1_dataset['exp_file'])
    theobeams = readOUTBEAMS(fe2o3_012_1x1_dataset['theo_file'])
    rpars = make_rparams(fe2o3_012_1x1_dataset,
                         'pendry',
                         expbeams,
                         theobeams)
    rpars.INTPOL_DEG = 5
    with pytest.raises(ParameterConflictError, match='R_FACTOR_LEGACY'):
        rfactor_section.run_new_rfactor(None,
                                        rpars,
                                        for_error=False,
                                        name='refcalc',
                                        theobeams=theobeams)


# A V0r (eV) inside _SHIFT_RANGE, but away from the best V0r of every
# R factor for the fe2o3-012-1x1 data (0.5 or 1.5 eV).
_FIXED_V0R = -2.0


def _run_section_at(dataset, r_type, monkeypatch, shift_range, **kwargs):
    """Return the Rparams after run_new_rfactor and the scan it made.

    Parameters
    ----------
    dataset : dict
        The beam data.
    r_type : str
        Name of the R factor.
    monkeypatch : pytest.MonkeyPatch
        Used to fix the beam correspondence and to capture the scan.
    shift_range : tuple of two float
        Start and stop of the V0r shifts, with step _STEP.
    **kwargs
        best_v0r : float, optional
            Value of best_v0r before the calculation.
        for_error : bool, optional
            Passed on to run_new_rfactor. Default is False.

    Returns
    -------
    rpars : Rparams
        The parameters after the calculation.
    scan : V0rScanResult
        The result of the only call to scan_v0r.
    """
    scans = []
    scan_v0r = rfactor_section.v0r_scan.scan_v0r

    def _record_scan(*args, **scan_kwargs):
        scans.append(scan_v0r(*args, **scan_kwargs))
        return scans[-1]

    monkeypatch.setattr(rfactor_section.v0r_scan, 'scan_v0r', _record_scan)
    theobeams = readOUTBEAMS(dataset['theo_file'])
    rpars = make_rparams(dataset,
                         r_type,
                         readOUTBEAMS(dataset['exp_file']),
                         theobeams)
    rpars.IV_SHIFT_RANGE = IVShiftRange(*shift_range, _STEP)
    rpars.best_v0r = kwargs.get('best_v0r')
    run_new_rfactor_with_correspondence(
        monkeypatch,
        dataset['beam_correspondence'],
        rpars,
        theobeams,
        for_error=kwargs.get('for_error', False),
        )
    assert len(scans) == 1
    return rpars, scans[0]


@pytest.mark.parametrize('r_type', _TYPES)
def test_run_new_rfactor_for_error_uses_best_v0r(r_type,
                                                 fe2o3_012_1x1_dataset,
                                                 monkeypatch):
    """Check that R for error curves is that at the fixed best V0r."""
    dataset = fe2o3_012_1x1_dataset
    rpars, scan = _run_section_at(dataset,
                                  r_type,
                                  monkeypatch,
                                  _SHIFT_RANGE,
                                  best_v0r=_FIXED_V0R,
                                  for_error=True)
    fixed, _ = _run_section_at(dataset,
                               r_type,
                               monkeypatch,
                               (_FIXED_V0R, _FIXED_V0R))
    best, _ = _run_section_at(dataset, r_type, monkeypatch, _SHIFT_RANGE)
    assert best.best_v0r != _FIXED_V0R  # Otherwise the test is moot
    assert rpars.best_v0r == _FIXED_V0R
    assert scan.best_v0r == _FIXED_V0R
    assert rpars.last_R == fixed.last_R
    assert rpars.last_R > best.last_R


@pytest.mark.parametrize('plot', [False, True], ids=['no-plot', 'plot'])
@pytest.mark.parametrize('for_error', [False, True],
                         ids=['normal', 'for-error'])
@pytest.mark.parametrize('r_type', _TYPES)
# pylint: disable-next=too-many-arguments,too-many-positional-arguments
def test_run_new_rfactor_no_warning(r_type,
                                    for_error,
                                    plot,
                                    fe2o3_012_1x1_dataset,
                                    monkeypatch,
                                    caplog):
    """Check that run_new_rfactor logs no warning for clean input."""
    calls = _capture_pdf_writer(monkeypatch)  # Fast, no PDF files
    caplog.set_level(logging.WARNING)
    _run_section(fe2o3_012_1x1_dataset,
                 r_type,
                 monkeypatch,
                 for_error=for_error,
                 plot=plot)
    assert not caplog.records
    assert bool(calls) is plot


@pytest.mark.skipif(not CAN_PLOT, reason='Needs matplotlib')
@pytest.mark.parametrize('r_type', _TYPES)
def test_run_new_rfactor_writes_pdfs(r_type,
                                     fe2o3_012_1x1_dataset,
                                     monkeypatch,
                                     tmp_path):
    """Check the PDFs written by run_new_rfactor with plotting on."""
    monkeypatch.chdir(tmp_path)
    results = _run_section(fe2o3_012_1x1_dataset,
                           r_type,
                           monkeypatch,
                           plot=True)
    assert results['r_overall'] == pytest.approx(
        _load_expected(fe2o3_012_1x1_dataset,
                       'run_new_rfactor',
                       r_type)['r_overall'],
        )
    assert (tmp_path / 'Rfactor_plots_refcalc.pdf').is_file()
    has_analysis = RFactorType[r_type.upper()] in (RFactorType.PENDRY,
                                                   RFactorType.SMOOTH)
    analysis_file = tmp_path / 'Rfactor_analysis_refcalc.pdf'
    assert analysis_file.is_file() == has_analysis


_PDF_R_TYPES = {
    'smooth': ('SMOOTH', True),
    'pendry': ('PENDRY', True),
    'r2': ('R2', False),
    'zj': ('ZJ', False),
    'r1': ('PENDRY', False),
    }


@pytest.mark.parametrize(
    ('which_r', 'expected'),
    _PDF_R_TYPES.items(),
    ids=_PDF_R_TYPES,
)
def test_rfactor_from_csv_pdf_r_type(which_r, expected,
                                     fe2o3_012_1x1_dataset, monkeypatch):
    """Check that the analysis PDF is for the R factor in `which_r`."""
    r_type_name, has_analysis = expected
    calls = _capture_pdf_writer(monkeypatch)
    _run_from_csv(fe2o3_012_1x1_dataset, which_r, plot_iv={'plot': True})
    assert len(calls) == 1
    assert calls[0]['r_factor_type'] == r_type_name
    assert bool(calls[0]['kwargs']['analysis_file']) is has_analysis


_CSV_HEADER = 'E, ( 1| 0), ( 0| 1)\n'
_FIRST_ENERGY = 30
_N_BEAMS = 3


def _csv_beams(energies):
    """Return the content of a CSV file of two beams at `energies`."""
    rows = (f'{energy:.2f}, {1 + i % 3}.0, {2 + i % 2}.0\n'
            for i, energy in enumerate(energies))
    return _CSV_HEADER + ''.join(rows)


_GOOD_CSV = _csv_beams(range(30, 50))
_BAD_CSV_INPUTS = {
    'no exp energies': (_CSV_HEADER, _GOOD_CSV),
    'one exp energy': (_csv_beams([30]), _GOOD_CSV),
    'no theo energies': (_GOOD_CSV, _CSV_HEADER),
    'one theo energy': (_GOOD_CSV, _csv_beams([30])),
    }


@pytest.mark.parametrize('contents', _BAD_CSV_INPUTS.values(),
                         ids=_BAD_CSV_INPUTS)
def test_rfactor_from_csv_rejects_bad_input(contents):
    """Check that empty or single-energy beams are refused."""
    with pytest.raises(ValueError, match=r'Need at least 2, found [01]'):
        from_ase.rfactor_from_csv(contents,
                                  1.0,
                                  beams_file_is_content=(True, True))


def test_rfactor_from_csv_rejects_no_theo_energies():
    """Check the message for theoretical beams without energies."""
    contents = _BAD_CSV_INPUTS['no theo energies']
    with pytest.raises(ValueError, match=r'theoretical energies.*found 0'):
        from_ase.rfactor_from_csv(contents,
                                  1.0,
                                  beams_file_is_content=(True, True))


def test_rfactor_from_csv_single_point_beam_is_skipped():
    """Check that a one-point beam does not affect the overall R."""
    header = 'E, ( 1| 0), ( 0| 1), ( 1| 1)\n'
    rows = [(f'{energy:.2f}, {1 + i % 3}.0, {2 + i % 2}.0, '
             f'{3 - i % 4}.0\n') for i, energy in enumerate(range(30, 50))]
    two_beams = _csv_beams(range(30, 50))
    # The third beam has data only at the first energy
    one_point = header + rows[0] + ''.join(
        row.rsplit(',', 1)[0] + ',\n' for row in rows[1:]
    )
    kwargs = {'beams_file_is_content': (True, True), 'which_r': 'pendry'}
    r_two, *_ = from_ase.rfactor_from_csv((two_beams, two_beams), 1.0,
                                          **kwargs)
    r_three, _, r_beams = from_ase.rfactor_from_csv((one_point, one_point),
                                                    1.0, **kwargs)
    assert r_three == pytest.approx(r_two)
    assert len(r_beams) == _N_BEAMS
    assert np.isnan(r_beams[2])


def _csv_with_short_beam(phase, n_points):
    """Return CSV content with a third beam of n_points energies.

    The third beam is left out if n_points is 0.
    """
    header = 'E, ( 1| 0), ( 0| 1)' + (', ( 1| 1)\n' if n_points else '\n')
    rows = []
    for i, energy in enumerate(range(_FIRST_ENERGY, 50)):
        row = (f'{energy:.2f}, {1 + (i + phase) % 4}.0, '
               f'{2 + (i * phase) % 3}.0')
        if n_points:
            row += f', {1 + (i + phase) % 3}.0' if i < n_points else ','
        rows.append(row + '\n')
    return header + ''.join(rows)


@pytest.mark.parametrize('n_points', [2, 3])
def test_rfactor_from_csv_uses_short_beams(n_points):
    """Check that beams with two or three points are not dropped."""
    # No shift, so that the overlap of a short beam cannot vanish
    kwargs = {
        'beams_file_is_content': (True, True),
        'which_r': 'pendry',
        'v0r_shift_range': (0, 0),
    }
    r_without, *_ = from_ase.rfactor_from_csv(
        (_csv_with_short_beam(0, 0), _csv_with_short_beam(1, 0)),
        1.0,
        **kwargs,
    )
    r_with, _, r_beams = from_ase.rfactor_from_csv(
        (_csv_with_short_beam(0, n_points), _csv_with_short_beam(1, n_points)),
        1.0,
        **kwargs,
    )
    assert np.isfinite(r_beams[2])
    assert r_with != pytest.approx(r_without)


def _csv_with_empty_beam(phase):
    """Return CSV content with a third beam that is all NaN."""
    header, *rows = _csv_with_short_beam(phase, 0).splitlines()
    return '\n'.join([header + ', ( 1| 1)',
                      *(row + ', NaN' for row in rows)]) + '\n'


@pytest.mark.parametrize('empty_side', ['exp', 'theo'])
def test_rfactor_from_csv_empty_beam_is_skipped(empty_side):
    """Check that an all-NaN beam does not affect the overall R."""
    kwargs = {'beams_file_is_content': (True, True), 'which_r': 'pendry'}
    r_without, *_ = from_ase.rfactor_from_csv(
        (_csv_with_short_beam(0, 0), _csv_with_short_beam(1, 0)),
        1.0,
        **kwargs,
    )
    full = (_csv_with_short_beam(0, 20), _csv_with_short_beam(1, 20))
    contents = (
        (_csv_with_empty_beam(0), full[1]) if empty_side == 'exp'
        else (full[0], _csv_with_empty_beam(1))
    )
    r_with, _, r_beams = from_ase.rfactor_from_csv(contents, 1.0, **kwargs)
    assert r_with == pytest.approx(r_without)
    assert len(r_beams) == _N_BEAMS
    assert np.isnan(r_beams[2])


_SKIPPED_MSG = 'Leaving out of the R factor'


def _csv_with_skipped_beams(phase, exp_side):
    """Return CSV content with two full and three short beams.

    The ( 1| 1) beam is all NaN and the ( 2| 0) beam has one point
    in the experiment; the ( 0| 2) beam has one point in theory.
    """
    header, *rows = _csv_with_short_beam(phase, 0).splitlines()
    rows = [
        row + (', NaN' if exp_side else f', {i + 1}.0')
        + (f', {i + 2}.0' if not exp_side or not i else ', NaN')
        + (f', {i + 3}.0' if exp_side or not i else ', NaN')
        for i, row in enumerate(rows)
    ]
    return '\n'.join([header + ', ( 1| 1), ( 2| 0), ( 0| 2)',
                      *rows]) + '\n'


def test_rfactor_from_csv_logs_skipped_beams_once(caplog):
    """Check that each beam too short for the R factor is named once."""
    caplog.set_level(logging.DEBUG)
    contents = (_csv_with_skipped_beams(0, exp_side=True),
                _csv_with_skipped_beams(1, exp_side=False))
    _, _, r_beams = from_ase.rfactor_from_csv(
        contents,
        1.0,
        beams_file_is_content=(True, True),
        which_r='pendry',
        v0r_shift_range=(-3, 3),  # Several shifts
    )
    assert np.isnan(r_beams[2:]).all()
    assert not [r for r in caplog.records if r.levelno >= logging.INFO]
    skipped = [r for r in caplog.records if _SKIPPED_MSG in r.getMessage()]
    assert len(skipped) == 1
    assert skipped[0].levelno == logging.DEBUG
    text = ''.join(r.getMessage() for r in skipped)
    expected = {'( 1| 1)': '0 experimental energies',
                '( 2| 0)': '1 experimental energy',
                '( 0| 2)': '1 theoretical energy'}
    for label, reason in expected.items():
        assert text.count(label) == 1
        assert f'{label}: {reason}' in text
    full_beams = _CSV_HEADER.strip().split(', ')[1:]
    assert not any(label in text for label in full_beams)


def test_rfactor_from_csv_logs_no_skipped_beams(caplog):
    """Check that no beam is reported as skipped for full beams."""
    caplog.set_level(logging.DEBUG)
    from_ase.rfactor_from_csv(
        (_csv_with_short_beam(0, 20), _csv_with_short_beam(1, 20)),
        1.0,
        beams_file_is_content=(True, True),
        which_r='pendry',
    )
    assert not [r for r in caplog.records
                if _SKIPPED_MSG in r.getMessage()]


_UNMATCHED_MSG = 'beams without'


@pytest.mark.parametrize(('extra_theo', 'level'), [(True, logging.DEBUG),
                                                (False, logging.WARNING)])
def test_rfactor_from_csv_logs_unmatched_beams(extra_theo, level, caplog):
    """Check the level of the record naming beams without a mate."""
    caplog.set_level(logging.DEBUG)
    two_beams = _csv_with_short_beam(0, 0)
    three_beams = _csv_with_short_beam(1, 20)
    contents = ((two_beams, three_beams) if extra_theo
                else (three_beams, two_beams))
    from_ase.rfactor_from_csv(contents,
                              1.0,
                              beams_file_is_content=(True, True))
    unmatched = [r for r in caplog.records
                 if _UNMATCHED_MSG in r.getMessage()]
    assert len(unmatched) == 1
    assert unmatched[0].levelno == level
    assert unmatched[0].getMessage().endswith(': ( 1| 1)')
    n_warnings = sum(r.levelno >= logging.WARNING for r in caplog.records)
    assert n_warnings == (level >= logging.WARNING)


def _smooth_csv(energies, phase):
    """Return the content of a CSV file with one smooth beam."""
    rows = (f'{energy:.2f}, {2 + np.sin((energy + phase) / 4):.6f}\n'
            for energy in energies)
    return 'E, ( 1| 0)\n' + ''.join(rows)


@pytest.mark.parametrize('which_r', ['pendry', 'r2'])
def test_rfactor_from_csv_large_v0r_range_is_not_nan(which_r):
    """Check that large shifts, with no overlap, do not give NaN."""
    # With a range of +-25 eV on 20 eV of data, the most extreme shifts
    # leave no overlap between experiment and theory.
    contents = (
        _smooth_csv(np.arange(40, 60.1, 1.0), 0.0),
        _smooth_csv(np.arange(40, 60.1, 0.5), 1.0),
    )
    results = [
        from_ase.rfactor_from_csv(
            contents,
            4.0,
            beams_file_is_content=(True, True),
            which_r=which_r,
            v0r_shift_range=shift_range,
        )
        for shift_range in ((-10, 10), (-25, 25))
    ]
    (r_small, v0r_small, *_), (r_large, v0r_large, *_) = results
    assert np.isfinite(r_large)
    assert r_large == pytest.approx(r_small)
    assert v0r_large == pytest.approx(v0r_small)


def test_rfactor_from_csv_rejects_only_single_point_beams():
    """Check that a file without any interpolable beam is refused."""
    rows = [f'{energy:.2f}, 1.0,\n' if energy == _FIRST_ENERGY
            else f'{energy:.2f}, , 2.0\n'
            for energy in range(_FIRST_ENERGY, 40)]
    only_single = _CSV_HEADER + ''.join(rows[:1] + rows[-1:])
    with pytest.raises(RfactorError, match='No experimental beam has'):
        from_ase.rfactor_from_csv((only_single, _GOOD_CSV),
                                  1.0,
                                  beams_file_is_content=(True, True))


def test_rfactor_from_csv_rejects_disjoint_energies():
    """Check that beams with no common energies are refused."""
    contents = _csv_beams(range(100, 110)), _GOOD_CSV
    with pytest.raises(RfactorError, match='No common energies'):
        from_ase.rfactor_from_csv(contents,
                                  1.0,
                                  beams_file_is_content=(True, True))


def test_rfactor_from_csv_auto_step():
    """Check that intpol_step=None uses the finer step of the data."""
    contents = _csv_beams(range(30, 50, 2)), _GOOD_CSV
    kwargs = {'beams_file_is_content': (True, True)}
    auto = from_ase.rfactor_from_csv(contents, 1.0, intpol_step=None,
                                     **kwargs)
    explicit = from_ase.rfactor_from_csv(contents, 1.0, intpol_step=1.0,
                                         **kwargs)
    assert auto[:2] == explicit[:2]
    np.testing.assert_array_equal(auto[2], explicit[2])


@pytest.mark.parametrize(('argument', 'value'),
                         [('intpol_deg', 5), ('return_beam_arrays', True)])
def test_rfactor_from_csv_deprecated_warns(argument, value, caplog):
    """Check the warning for an ignored argument names the function."""
    caplog.set_level(logging.WARNING)
    from_ase.rfactor_from_csv((_GOOD_CSV, _GOOD_CSV),
                              1.0,
                              beams_file_is_content=(True, True),
                              **{argument: value})
    expected = f'{argument} is deprecated in rfactor_from_csv()'
    assert any(expected in record.getMessage()
               for record in caplog.records)
