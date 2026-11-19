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
from dataclasses import replace
import logging
from typing import Any
from typing import Callable

import numpy as np

from viperleed.calc.files.iorfactor import beamlist_to_array
from viperleed.calc.lib import spline_interpolation
from viperleed.calc.lib.rfactor.utils import average_beam_array

logger = logging.getLogger(__name__)


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
    best_v0r: float = None
    best_r: float = None

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


def _make_spline(grid, intensities):
    """Return a cached not-a-knot spline of ragged beam data."""
    spline = spline_interpolation.interpolate_ragged_array(
        x=grid,
        y=intensities,
        bc_type='not-a-knot',  # TODO: from Rfactor parameter?
    )
    return spline_interpolation.CachedSpline(spline)


def _find_best_shift(result, iv_shift):
    """Return the shift in iv_shift with the lowest R, and that R."""
    logger.debug('Sweeping V0r...')
    # sweep V0r shifts - TODO: continuous shifts?
    shifts = np.arange(
        iv_shift.start,
        iv_shift.stop + iv_shift.step,
        iv_shift.step,
    )
    r_values = np.array([
        result.r_factor(shift_2nd_spline=shift)
        for shift in shifts
    ])
    best_index = np.argmin(r_values)
    best_shift = shifts[best_index]
    best_r = r_values[best_index]
    logger.info(f'Best inner potential shift: {best_shift:.2f} eV '
                f'with R = {best_r:.4f}')
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
    out_grid = np.arange(
        theo_range.min, theo_range.max + 0.1 * intpol_step, intpol_step
    )

    logger.debug('Constructing splines from data...')
    exp_energies, _, _, exp_intensities = beamlist_to_array(exp_beams)
    exp_spline = _make_spline(exp_energies, exp_intensities)

    theo_grid, _, _, theo_intensities = beamlist_to_array(theo_beams)
    theo_intensities = average_beam_array(theo_intensities,
                                          beam_correspondence)
    result = V0rScanResult(
        r_func,
        v0_imag,
        intpol_step,
        out_grid,
        exp_spline,
        _make_spline(theo_grid, theo_intensities),
    )

    best_shift, best_r = _find_best_shift(result, iv_shift)
    theo_spline = result.theo_spline
    if best_shift != 0:
        # re-build theoretical splines one more time on the shifted grid
        theo_spline = _make_spline(theo_grid + best_shift, theo_intensities)
    return replace(result, theo_spline=theo_spline,
                   best_v0r=best_shift, best_r=best_r)
