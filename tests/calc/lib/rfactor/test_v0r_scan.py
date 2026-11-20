"""Tests for module v0r_scan of viperleed.calc.lib.rfactor."""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-20'
__license__ = 'GPLv3+'

import numpy as np

from viperleed.calc.classes.beam import Beam
from viperleed.calc.classes.rparams.special.energy_range import (
    IVShiftRange,
    TheoEnergies,
)
from viperleed.calc.files import iorfactor
from viperleed.calc.lib.rfactor import v0r_scan
from viperleed.calc.lib.rfactor.selection import select_rfactor_by_name

_STEP = 1.0
_SHIFT = 5.0   # Larger than one theory step


def _make_beam(energies, intensities):
    beam = Beam((1, 0))
    beam.intens = dict(zip(energies, intensities))
    return beam


def _intensity(energies):
    return 2 + np.sin(energies / 4) + 0.5 * np.cos(energies / 7)


def _scan_fixed_shift(exp_tail_factor):
    """Scan a single shift; the experiment differs at the overlap top.

    The experiment is the theory shifted by _SHIFT, but its top
    energies (those at which the theory is cut at the base) are
    scaled by exp_tail_factor.
    """
    theo_energies = np.arange(100, 131) * _STEP
    exp_energies = np.arange(100, 141) * _STEP
    exp_intensity = _intensity(exp_energies - _SHIFT)
    exp_intensity[exp_energies > 131 + 0.5] *= exp_tail_factor
    theo_beams = [_make_beam(theo_energies, _intensity(theo_energies))]
    exp_beams = [_make_beam(exp_energies, exp_intensity)]

    iv_shift = IVShiftRange(_SHIFT, _SHIFT, _STEP)
    theo_range = iorfactor.find_theory_range(
        TheoEnergies.from_sorted_grid(theo_energies),
        iorfactor.experiment_energy_range(exp_beams),
        iv_shift,
        n_expand=1,
    )
    return v0r_scan.scan_v0r(
        select_rfactor_by_name('r_1'),
        4.0,
        exp_beams,
        theo_beams,
        np.array([0]),
        iv_shift,
        theo_range,
    )


def test_scan_covers_overlap_at_large_shift():
    """Check that the top of the overlap at a large shift counts."""
    r_plain = _scan_fixed_shift(1.0).best_r
    r_modified = _scan_fixed_shift(1.5).best_r
    assert r_modified > r_plain + 0.01
