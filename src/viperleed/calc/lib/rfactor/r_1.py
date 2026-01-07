"""viperleed.calc.lib.rfactor.r_2."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2026-01-07'
__license__ = 'GPLv3+'

from .groups import group_rfactors
from .utils import nansum_trapezoid


def R_2(
    theo_spline, v0_imag, energy_step, energy_grid, exp_spline, groups=None
):
    """
    Notes
    -----
    Does not implement shifting intensity values abouve zero."""

    # experiment
    exp_intensity = exp_spline(energy_grid)

    # theory
    theo_intensity = theo_spline(energy_grid)

    # calculate normalization for each beam
    beam_normalization = nansum_trapezoid(
        exp_intensity, energy_step, axis=0
    ) / nansum_trapezoid(theo_intensity, energy_step, axis=0)

    numerators = nansum_trapezoid(
        abs(exp_intensity - beam_normalization * theo_intensity),
        energy_step,
        axis=0,
    )
    denominators = nansum_trapezoid(exp_intensity, energy_step, axis=0)

    # calculate R-factor with requested grouping
    return group_rfactors(numerators, denominators, groups=groups)
