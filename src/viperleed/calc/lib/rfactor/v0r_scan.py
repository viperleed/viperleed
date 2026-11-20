"""Module v0r_scan of viperleed.calc.lib.rfactor.

Host-side driver of an R-factor calculation: builds splines from
beams, scans the relative energy shift (V0r) between experiment and
theory, and exposes the R-factor at the best shift. This is NumPy
code operating on ragged data. It is not meant to be jit-compiled;
only the R-factor functions it calls are JAX-compatible.
"""

__authors__ = (
    'Florian Kraushofer (@fkraushofer)',
)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-11-19'
__license__ = 'GPLv3+'

from dataclasses import dataclass
import functools
import logging
from typing import Any
from typing import Callable

import numpy as np

from viperleed.calc.files import iorfactor
from viperleed.calc.lib import spline_interpolation
from viperleed.calc.lib.rfactor.utils import average_beam_array

logger = logging.getLogger(__name__)

# Fewest energies a beam needs for a cubic spline
_MIN_POINTS = 2


@dataclass(frozen=True)
# pylint: disable-next=too-many-instance-attributes  # Result record
class V0rScanResult:
    """Splines and best shift found by a scan of V0r.

    Attributes
    ----------
    r_func : callable
        The R-factor function used for the scan.
    v0_imag : float
        Imaginary part of the inner potential (eV).
    intpol_step : float
        Energy step (eV) of out_grid and of the V0r scan.
    out_grid : numpy.ndarray
        Energies on which the splines are compared.
    exp_spline : CachedSpline
        Spline of the experimental beams.
    theo_spline : CachedSpline
        Spline of the (averaged) theoretical beams, already
        shifted by best_v0r.
    best_v0r : float
        The shift (eV) with the lowest overall R factor.
    best_r : float
        The overall R factor at best_v0r, as found during the scan.
    """

    r_func: Callable
    v0_imag: float
    intpol_step: float
    out_grid: Any
    exp_spline: Any
    theo_spline: Any
    best_v0r: float
    best_r: float

    def r_factor(self, groups=None, **kwargs):
        """Return the R factor between the stored splines.

        Parameters
        ----------
        groups : str, Sequence or None, optional
            How to group beams. Passed to r_func. Default is None
            (a single overall R factor).
        **kwargs
            Further keyword arguments for r_func (e.g., num_groups,
            shift_2nd_spline).

        Returns
        -------
        r_factor : float or numpy.ndarray
            What r_func returns.
        """
        return self.r_func(
            self.v0_imag,
            self.intpol_step,
            self.out_grid,
            data_spline_1=self.exp_spline,
            data_spline_2=self.theo_spline,
            groups=groups,
            **kwargs,
        )


def _make_spline(grid, intensities, what):
    """Return a cached not-a-knot spline of ragged beam data.

    Beams with fewer than two points cannot be interpolated. They
    are turned into all-NaN beams, i.e., beams without overlap.

    Raises
    ------
    RfactorError
        If no beam has at least two points.
    """
    n_points = np.count_nonzero(~np.isnan(intensities), axis=0)
    if not np.any(n_points >= _MIN_POINTS):
        msg = (
            f'No {what} beam has at least {_MIN_POINTS} energies. '
            'Cannot interpolate.'
        )
        raise iorfactor.RfactorError(msg)
    intensities = np.where(n_points >= _MIN_POINTS, intensities, np.nan)
    spline = spline_interpolation.interpolate_ragged_array(
        x=grid,
        y=intensities,
        bc_type='not-a-knot',  # Cubic only
    )
    return spline_interpolation.CachedSpline(spline)


def _experimental_axis(exp_energies, theo_range, iv_shift):
    """Return the energies at which experiment and theory are compared.

    The theory is evaluated at these energies minus the shift. The grid
    thus covers theo_range moved by every shift in iv_shift, i.e., the
    whole overlap with experiment for each shift, restricted to the
    energies of the experiment.
    """
    step = iv_shift.step
    axis = np.arange(
        theo_range.min + iv_shift.start,
        theo_range.max + iv_shift.stop + 0.1 * step,
        step,
    )
    tolerance = 0.1 * step
    return axis[(axis >= exp_energies[0] - tolerance)
                & (axis <= exp_energies[-1] + tolerance)]


def _find_best_shift(r_at_shift, iv_shift):
    """Return the shift in iv_shift with the lowest R, and that R.

    Shifts at which R is NaN (e.g., because experiment and theory do
    not overlap there) are ignored.

    Raises
    ------
    RfactorError
        If R is NaN at all shifts.
    """
    logger.debug('Sweeping V0r...')
    # sweep V0r shifts - TODO: continuous shifts?
    shifts = np.arange(
        iv_shift.start,
        iv_shift.stop + iv_shift.step,
        iv_shift.step,
    )
    # Shifts without overlap give 0/0, i.e., NaN: this is expected
    with np.errstate(invalid='ignore'):
        r_values = np.array([
            r_at_shift(shift_2nd_spline=shift)
            for shift in shifts
        ])
    if np.all(np.isnan(r_values)):
        msg = (
            'The R factor is undefined (NaN) at all V0r shifts in '
            f'[{iv_shift.start}, {iv_shift.stop}] eV. Experiment and '
            'theory may have no energies in common at any of them.'
        )
        raise iorfactor.RfactorError(msg)
    best_index = np.nanargmin(r_values)
    best_shift = shifts[best_index]
    best_r = r_values[best_index]
    logger.debug('Best inner potential shift: %.2f eV with R = %.4f',
                 best_shift, best_r)
    return best_shift, best_r


# Plain inputs by design: no Rparams
# pylint: disable-next=too-many-arguments,too-many-positional-arguments,too-many-locals
def scan_v0r(r_func, v0_imag, exp_beams, theo_beams,  # noqa: PLR0913, PLR0917
             beam_correspondence, iv_shift, theo_range):
    """Find the V0r shift that minimizes the R factor.

    Parameters
    ----------
    r_func : callable
        The R-factor function (see lib.rfactor.selection).
    v0_imag : float
        Imaginary part of the inner potential (eV).
    exp_beams : Sequence
        Experimental Beam objects.
    theo_beams : Sequence
        Theoretical Beam objects.
    beam_correspondence : Sequence
        For each experimental beam, which theoretical beams are
        averaged into it, in average_beam_array format.
    iv_shift : IVShiftRange
        The shifts (eV) to scan. Its step is also used as the step
        of the output grid. Must have a defined step.
    theo_range : TheoEnergies
        The range of theory energies on which R is evaluated.

    Returns
    -------
    result : V0rScanResult
        Holds the best shift and R factor, and the splines at the
        best shift, from which further R factors can be evaluated.
    """
    intpol_step = iv_shift.step

    logger.debug('Constructing splines from data...')
    exp_energies, _, _, exp_intensities = iorfactor.beamlist_to_array(
        exp_beams
    )
    exp_spline = _make_spline(exp_energies, exp_intensities, 'experimental')
    theo_grid, _, _, theo_intensities = iorfactor.beamlist_to_array(theo_beams)
    theo_intensities = average_beam_array(theo_intensities,
                                          beam_correspondence)
    theo_spline = _make_spline(theo_grid, theo_intensities, 'theoretical')

    out_grid = _experimental_axis(exp_energies, theo_range, iv_shift)
    best_shift, best_r = _find_best_shift(
        functools.partial(
            r_func,
            v0_imag,
            intpol_step,
            out_grid,
            data_spline_1=exp_spline,
            data_spline_2=theo_spline,
            groups=None,
        ),
        iv_shift,
    )
    if best_shift:
        theo_spline = _make_spline(
            theo_grid + best_shift,
            theo_intensities,
            'theoretical',
        )
    return V0rScanResult(
        r_func,
        v0_imag,
        intpol_step,
        out_grid,
        exp_spline,
        theo_spline,
        best_shift,
        best_r,
    )
