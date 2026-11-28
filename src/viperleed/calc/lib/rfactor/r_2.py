"""viperleed.calc.lib.rfactor.r_2."""

__authors__ = ('Alexandra Mia Imre (@alexmiame)',
               'Florian Kraushofer (@fkraushofer)',)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2026-01-07'
__license__ = 'GPLv3+'

from . import utils
from .groups import group_beam_rfactors


def r_2(
    v0_imag,
    energy_step,
    energy_grid,
    data_spline_1=None,
    data_and_derivatives_1=None,
    data_spline_2=None,
    data_and_derivatives_2=None,
    shift_2nd_spline=0.0,  # only available if passed as spline
    **kwargs,
):
    """
    Calculate the R2 R factor.

    For each beam, with experimental intensity I_e and theoretical
    intensity I_t,

        R2 = int (I_e - c * I_t)^2 dE / int I_e^2 dE ,

    where c = int I_e dE / int I_t dE is the scale factor of the beam.
    All integrals (trapezoidal rule) run over the energy range in which
    both datasets are defined. This is the definition of Van Hove,
    Weinberg and Chan, "Low-Energy Electron Diffraction" (Springer,
    1986), https://doi.org/10.1007/978-3-642-82721-1, also used in
    Imre et al., arXiv:2511.05448 (Eq. 1). The R2 of the TensErLEED
    R-factor program (rfactor.f, used with R_FACTOR_LEGACY) includes
    an additional factor 1/2, i.e., it is half of the value calculated
    here. The R2 of the TensErLEED structure search (lib.search.f)
    has no such factor.

    The overall R2 value, and the value of each group of beams, is the
    mean of the per-beam R2 weighted with the length (n - 1) * step of
    the overlap of the beam, where n is the number of overlapping
    points. Beams without overlap, and beams whose R2 is undefined
    (zero denominator, e.g., I_e = 0 on the overlap), are excluded
    from the overall and group values.

    Uses two sets of beam data, either passed as splines or as
    already calculated arrays of intensities.

    For each dataset, pass either data_spline or data_and_derivatives, but not
    both. Using splines for one and pre-computed data for the other is allowed.
    Note that data_and_derivatives are expected to be 3-tuples for uniform call
    signatures with R-factors that also use the derivatives; however, the
    derivatives are not used and can be left empty.

    Parameters
    ----------
    v0_imag : float
        The imaginary part of the inner potential.
    energy_step : float
        The energy step of the data grid.
    energy_grid : array
        The grid on which data should be evaluated. All evaluations will ignore
        regions in which either of the datasets is nan.
    data_spline_1, data_spline_2 : arrays of splines
        Splines of the data, which will be evaluated on energy_grid. Ignored if
        the respective data_and_derivatives are passed.
    data_and_derivatives_1, data_and_derivatives_2 : 3-tuples of float arrays
        Pre-evaluated tuples (intensity, 1st derivative, 2nd derivative) at the
        energy_grid points. This avoids re-calculation of splines and is more
        JAX-friendly. The derivatives are not actually used, so passing
        (intensity, None, None) would be acceptable.
    shift_2nd_spline : float, optional
        Evaluate the 2nd spline on a shifted energy grid. This is meant for
        testing V0r variations. Note that this is NOT available when the 2nd
        dataset is passed as data_and_derivatives.

    Raises
    ------
    TypeError
        If neither data_spline nor data_and_derivatives is passed for
        either dataset, or if both are passed for the same dataset.
    """
    # Get data either as splines or as pre-computed arrays (mainly for JAX)
    data_1_intensity, *_ = utils.coerce_to_energy_grid(
        data_and_derivatives_1, data_spline_1, energy_grid, derivs=0
    )
    # 2nd data uses shifted grid if spline
    data_2_intensity, *_ = utils.coerce_to_energy_grid(
        data_and_derivatives_2,
        data_spline_2,
        energy_grid - shift_2nd_spline,
        derivs=0,
    )

    # all integrals are restricted to the overlap of the two datasets
    overlap = utils.overlap_mask(data_1_intensity, data_2_intensity)

    # calculate normalization for each beam (nan without overlap)
    beam_normalization = utils.overlap_trapezoid(
        data_1_intensity, energy_step, overlap, axis=0
    ) / utils.nonzero_or_nan(
        utils.overlap_trapezoid(data_2_intensity, energy_step, overlap, axis=0)
    )

    numerators = utils.overlap_trapezoid(
        (data_1_intensity - beam_normalization * data_2_intensity) ** 2,
        energy_step,
        overlap,
        axis=0,
    )

    # Groups and the overall value are the means over beams weighted
    # by the length of the overlap.
    beam_rfactors = numerators / utils.nonzero_or_nan(
        utils.overlap_trapezoid(
            data_1_intensity**2, energy_step, overlap, axis=0
        )
    )
    return group_beam_rfactors(
        beam_rfactors, overlap, energy_step, **kwargs
    )
