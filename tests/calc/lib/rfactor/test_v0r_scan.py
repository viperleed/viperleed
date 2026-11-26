"""Tests for module v0r_scan of viperleed.calc.lib.rfactor."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-20'
__license__ = 'GPLv3+'

import numpy as np
import pytest

from viperleed.calc.classes.beam import Beam
from viperleed.calc.classes.rparams.special import energy_range
from viperleed.calc.files import iorfactor
from viperleed.calc.lib.rfactor import v0r_scan
from viperleed.calc.lib.rfactor.selection import select_rfactor_by_name

_STEP = 1.0
_SHIFT = 5.0   # Larger than one theory step
_EDGE_CUTOFF = 2.5   # Experiment altered this far inside the overlap


def _make_beam(energies, intensities):
    beam = Beam((1, 0))
    beam.intens = dict(zip(energies, intensities))
    return beam


def _intensity(energies):
    return 2 + np.sin(energies / 4) + 0.5 * np.cos(energies / 7)


def _scan_fixed_shift(sign, tail_factor, r_name):
    """Scan a single shift; the experiment differs at the overlap edge.

    The experiment is the theory shifted by sign * _SHIFT, but its
    energies within _EDGE_CUTOFF of one end of the overlap (the top
    for sign > 0, the bottom otherwise) are scaled by tail_factor.
    These are the energies that are lost if the energy grid is cut at
    that end.
    """
    shift = sign * _SHIFT
    theo_energies = np.arange(100, 131) * _STEP
    exp_energies = np.arange(90, 141) * _STEP
    exp_intensity = _intensity(exp_energies - shift)
    if sign > 0:
        edge = theo_energies[-1] + shift - _EDGE_CUTOFF
        exp_intensity[exp_energies > edge] *= tail_factor
    else:
        edge = theo_energies[0] + shift + _EDGE_CUTOFF
        exp_intensity[exp_energies < edge] *= tail_factor
    theo_beams = [_make_beam(theo_energies, _intensity(theo_energies))]
    exp_beams = [_make_beam(exp_energies, exp_intensity)]

    iv_shift = energy_range.IVShiftRange(shift, shift, _STEP)
    theo_range = iorfactor.find_theory_range(
        energy_range.TheoEnergies.from_sorted_grid(theo_energies),
        iorfactor.experiment_energy_range(exp_beams),
        iv_shift,
        n_expand=1,
    )
    return v0r_scan.scan_v0r(
        select_rfactor_by_name(r_name),
        4.0,
        exp_beams,
        theo_beams,
        np.array([0]),
        iv_shift,
        theo_range,
    )


@pytest.mark.parametrize('r_name', ['r2', 'pendry'])
@pytest.mark.parametrize('sign', [1, -1], ids=['top', 'bottom'])
def test_scan_covers_overlap_at_large_shift(sign, r_name):
    """Check that the edges of the overlap at a large shift count."""
    r_plain = _scan_fixed_shift(sign, 1.0, r_name).best_r
    r_modified = _scan_fixed_shift(sign, 1.5, r_name).best_r
    assert r_modified > r_plain + 0.01


def test_find_best_shift_ignores_nan():
    """Check that shifts with an undefined R are not selected."""
    iv_shift = energy_range.IVShiftRange(-2, 2, 1)
    r_values = {-2: np.nan, -1: np.nan, 0: 0.5, 1: 0.2, 2: np.nan}
    best_shift, best_r = v0r_scan._find_best_shift(  # pylint: disable=protected-access
        lambda shift_2nd_spline: r_values[shift_2nd_spline],
        iv_shift,
    )
    assert best_shift == pytest.approx(1)
    assert best_r == pytest.approx(0.2)


def test_find_best_shift_rejects_all_nan():
    """Check that an error is raised if R is NaN at every shift."""
    iv_shift = energy_range.IVShiftRange(-2, 2, 1)
    with pytest.raises(iorfactor.RfactorError, match='undefined'):
        v0r_scan._find_best_shift(  # pylint: disable=protected-access
            lambda **_: np.nan,
            iv_shift,
        )
