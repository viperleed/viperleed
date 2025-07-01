"""Module interpolation.

This module is a reworking of scipy's and my Bspline interpolation methods.
It can interpolate functions efficiently and in a JAX-compatible way.
"""

__authors__ = ('Alexander M. Imre (@amimre)',)
__created__ = '2024-02-19'

from abc import ABC, abstractmethod
from functools import partial

import interpax
import jax
import numpy as np
from jax import numpy as jnp
from jax.tree_util import register_pytree_node_class
from scipy import interpolate

# The two functions below are the only ones that are currently used
# Everything else is done by the interpax module.


def make_1d_ragged_cubic_spline(
    x, y, axis=0, bc_type='not-a-knot', extrapolate=False
):
    """Construct a piecewise cubic spline interpolator with ragged edges.

    The interpolator uses a cubic spline to interpolate data.
    """
    if x.ndim > 1 or y.ndim > 1:
        raise ValueError('x and y must be 1-dimensional arrays.')
    y_mask = jnp.isnan(y)
    x_subarray, y_subarray = x[~y_mask], y[~y_mask]
    start_index = jnp.where(~y_mask)[0][0]
    subarray_spline = interpax.CubicSpline(
        x_subarray, y_subarray, axis, bc_type, extrapolate, check=False
    )

    return subarray_spline, start_index


def interpolate_ragged_array(
    x, y, axis=0, bc_type='not-a-knot', extrapolate=False
):
    all_coeffs = np.full((4, y.shape[0], y.shape[1]), fill_value=jnp.nan)
    for dim in range(y.shape[1]):
        _y = y[:, dim]
        all_nans = jnp.all(jnp.isnan(_y))
        if all_nans:
            continue
        spline, start_id = make_1d_ragged_cubic_spline(
            x, y[:, dim], axis=0, bc_type=bc_type, extrapolate=None
        )
        all_coeffs[:, start_id : start_id + spline.c.shape[1], dim] = spline.c
    return interpax.PPoly.construct_fast(all_coeffs, x)


##############################################################################
# Currently unused code
##############################################################################


# Basis transformation from cradinal B-spline to cardinal piecewise polynomial
# in the form as given by de Boor, A practical guide to splines 1978, p.324
B_TO_PP_SPLINE_BASIS_TRANSFORMATION = np.array(
    [
        [1 / 6, 2 / 3, 1 / 6, 0.0],
        [-1 / 2, 0.0, 1 / 2, 0.0],
        [1.0, -2.0, 1.0, 0.0],
        [-1.0, 3.0, -3.0, 1.0],
    ]
)


def translate_cubic_pp_spline_coeffs(s):
    """Return a transformation matrix that translates spline coeffiecients.

    The return transformation can be applied to the coefficients (a,b,c,d) of a
    cubic (spline) in the piecewise-polynomial basis, i.e.:
    f_i(x) = a*x**3 + b*x**2 + c*x + d
    The resulting set of coefficients yields the coefficents for g(x) = f(x-s).
    Note that for s->0, the transformation approaches unity.
    Note also that for a piecewise polynomial splines, this shift in only valid
    while x+s is in the same knot-interval as the x.

    Parameters
    ----------
    d : float
        The amount to translate the coefficients.

    Returns
    -------
    ndarray
        The transformation in the form of a (4x4) matrix.
    """
    return np.array(
        [
            [1.0, 0.0, 0.0, 0.0],
            [3 * s, 1.0, 0.0, 0.0],
            [3 * s**2, 2 * s, 1.0, 0.0],
            [s**3, s**2, s, 1.0],
        ]
    )


@register_pytree_node_class
class CardinalSplineInterpolator(ABC):
    expected_kws = (
        'origin_grid',
        'target_grid',
        'intpol_deg',
        'knots',
        'intervals',
        'de_boor_coeffs',
        'full_colloc_matrix',
        'inv_colloc_matrix',
    )

    def __init__(self, origin_grid, target_grid, intpol_deg):
        self.origin_grid = origin_grid
        self.intpol_step = target_grid[1] - target_grid[0]
        self.orig_step = origin_grid[1] - origin_grid[0]
        self.target_grid = target_grid
        self.intpol_deg = intpol_deg
        self.knots = self._get_knots()
        self.intervals = self._find_intervals()
        # TODO: If we find a suitable solver, we could also use the banded
        # and pre-factorized collocation matrix

        # calculate De Boor coefficients
        de_boor_coeffs = [
            self._calc_de_boor_coeffs(deriv_order)
            for deriv_order in range(intpol_deg)
        ]
        self.de_boor_coeffs = jnp.array(de_boor_coeffs)

        # calcualte collocation matrix
        self.full_colloc_matrix = self.calculate_colloc_matrix()

        # Invert collocation matrix – this reduces the B-spline interpolation
        # to a simple matrix- vector multiplication
        self.inv_colloc_matrix = jnp.linalg.inv(self.full_colloc_matrix)

    def is_compatible(self, other):
        if not isinstance(other, type(self)):
            return False
        return (
            jnp.all(self.target_grid == other.target_grid)
            and self.intpol_deg == other.intpol_deg
        )

    @property
    def x_diffs(self):
        x_diffs = self.knots[self.intervals + 1] - self.target_grid
        return jnp.array([x_diffs**i for i in range(1, self.intpol_deg + 1)])

    @abstractmethod
    def _get_knots(self):
        raise NotImplementedError

    def _calc_de_boor_coeffs(self, deriv_order):
        """Calculate the De Boor coeffs for the given knots and target grid"""
        de_boor_coeffs = np.zeros((self.intpol_deg + 1, self.target_grid.size))
        for i, (interval, new_x) in enumerate(
            zip(self.intervals, self.target_grid)
        ):
            beta_coeffs = interpolate._bspl.evaluate_all_bspl(
                self.knots, self.intpol_deg, new_x, interval, deriv_order
            )
            de_boor_coeffs[:, i] = beta_coeffs
        return de_boor_coeffs

    def _find_intervals(self):
        """Return index of interval in knots that contains x_val"""
        # raise if knots are not sorted
        if not jnp.all(self.knots[:-1] <= self.knots[1:]):
            raise ValueError('knots must be sorted')
        return (
            jnp.clip(
                jnp.searchsorted(self.knots, self.target_grid, side='left'),
                a_min=self.intpol_deg + 1,
                a_max=self.knots.size - self.intpol_deg - 1,
            )
            - 1
        )

    def _banded_colloc_matrix_to_full(self, banded_colloc_matrix):
        kl = self.intpol_deg
        center_row = 2 * kl
        full_matrix = np.diag(banded_colloc_matrix[center_row, :])
        for k in range(1, kl + 1):
            full_matrix += np.diag(banded_colloc_matrix[center_row - k, k:], k)
            full_matrix += np.diag(
                banded_colloc_matrix[center_row + k, :-k], -k
            )
        return full_matrix

    def tree_flatten(self):
        # build kw_dict to rebuild the object
        kw_dict = {}
        for kw in self.expected_kws:
            kw_dict[kw] = getattr(self, kw)
        return None, kw_dict

    @classmethod
    def tree_unflatten(cls, aux_data, children):
        del children  # unused in this class
        new_object = cls.__new__(cls)
        new_object._rebuild_from_dict(aux_data[1])
        return new_object


class CardinalNotAKnotSplineInterpolator(CardinalSplineInterpolator):
    def _get_knots(self):
        return interpolate._bsplines._not_a_knot(
            self.origin_grid, self.intpol_deg
        )

    def calculate_colloc_matrix(self):
        # number of derivatives at boundaries (not-a-knot bc)
        deriv_l = None
        deriv_l_ords, deriv_l_vals = interpolate._bsplines._process_deriv_spec(
            deriv_l
        )
        deriv_r = None
        deriv_r_ords, deriv_r_vals = interpolate._bsplines._process_deriv_spec(
            deriv_r
        )

        nleft = deriv_l_ords.shape[0]
        nright = deriv_r_ords.shape[0]

        # basic collocation matrix
        kl = ku = self.intpol_deg
        nt = self.knots.size - self.intpol_deg - 1
        banded_colloc_matrix = np.zeros(
            (2 * kl + ku + 1, nt), dtype=np.float64, order='F'
        )
        interpolate._bspl._colloc(
            np.array(self.origin_grid, dtype=np.float64),
            np.array(self.knots, dtype=np.float64),
            self.intpol_deg,
            banded_colloc_matrix,
        )

        # derivatives at boundaries
        if nleft > 0:
            interpolate._bspl._handle_lhs_derivatives(
                self.knots,
                self.intpol_deg,
                self.origin_grid[0],
                banded_colloc_matrix,
                kl,
                ku,
                deriv_l_ords.astype(np.dtype('long')),
            )
        if nright > 0:
            interpolate._bspl._handle_lhs_derivatives(
                self.knots,
                self.intpol_deg,
                self.origin_grid[-1],
                banded_colloc_matrix,
                kl,
                ku,
                deriv_r_ords.astype(np.dtype('long')),
                offset=nt - nright,
            )

        # get full collocation matrix
        return self._banded_colloc_matrix_to_full(banded_colloc_matrix)

    # TODO: The below functions could (and probably should be) interpolator class
    #       methods. However, we need to figure out how to make this work with JAX.

    def get_bspline_coeffs(self, rhs):
        """Return the coefficients of the B-spline interpolant.

        Solves the linear system lhs * coeffs = rhs for coeffs.
        """
        # TODO: we could do this more efficiently. One easy improvement would be to
        # pre-factorize lhs by splitting .solve() into .lu_factor() and .lu_solve()
        # parts. Only the solve part depends on the right hand side.

        rhs_nan_mask = jnp.isnan(rhs)  # get mask for NaN values
        # solve the linear system
        raw_spline_coeffs = self.inv_colloc_matrix @ jnp.nan_to_num(rhs, 0)
        # mask NaN values in the result such that they are not used in the interpolation
        return jnp.where(
            rhs_nan_mask, jnp.nan, raw_spline_coeffs
        )

    @partial(jax.vmap, in_axes=(None, 1))
    def convert_b_to_pp_spline_coeffs(self, bspline_coeffs):
        """Converts B-spline to piecewise polynomial coefficents."""
        # pad with knot values
        _bspline_coeffs = jnp.pad(bspline_coeffs[1:], (3, 3), 'edge')

        a = (
            jnp.convolve(
                _bspline_coeffs,
                B_TO_PP_SPLINE_BASIS_TRANSFORMATION[3, :],
                'same',
            )
            / self.orig_step**3
        )
        a = a / 6

        b = (
            jnp.convolve(
                _bspline_coeffs,
                B_TO_PP_SPLINE_BASIS_TRANSFORMATION[2, :],
                'same',
            )
            / self.orig_step**2
        )
        b = (b - 6 * a) / 2

        c = (
            jnp.convolve(
                _bspline_coeffs,
                B_TO_PP_SPLINE_BASIS_TRANSFORMATION[1, :],
                'same',
            )
            / self.orig_step
        )
        c = c - 3 * a - 2 * b

        d = jnp.convolve(
            _bspline_coeffs, B_TO_PP_SPLINE_BASIS_TRANSFORMATION[0, :], 'same'
        )
        d = d - a - b - c
        return jnp.array([a, b, c, d])[
            :, :
        ]  # cut off last three coeffs from convolution
        # pp_spline_coeffs = jnp.pad(pp_spline_coeffs, ((0,0), (3,3)), 'edge')

    from functools import partial

    @partial(jax.vmap, in_axes=(None, 0, None, None))
    def evaluate_pp_spline_coeffs(self, pp_spline_coeffs, deriv=0, shift=0.0):
        knot_shift, frac_shift = divmod(shift, self.intpol_step)
        knot_shift = int(knot_shift)

        # a,b,c,d = pp_spline_coeffs
        shift_matrix = translate_cubic_pp_spline_coeffs(-frac_shift)
        a, b, c, d = shift_matrix @ pp_spline_coeffs
        # multiply wrapped coefficents by NaNs to remove them.
        n_intervals = len(self.intervals)
        invalid_ids = jnp.arange(n_intervals)
        invalidation_mask = jnp.logical_or(
            invalid_ids < -knot_shift, invalid_ids >= n_intervals - knot_shift
        )
        invalidation_mask = jnp.where(invalidation_mask, jnp.nan, 1.0)

        _intervals = jnp.roll(self.intervals + 1, -knot_shift)
        _x_diffs = jnp.roll(self.x_diffs, -knot_shift)
        _x_diffs = invalidation_mask * _x_diffs
        x, x2, x3 = _x_diffs

        if deriv == 0:
            return (
                a[_intervals] * x3
                + b[_intervals] * x2
                + c[_intervals] * x
                + d[_intervals] * jnp.ones_like(x)
            )
        if deriv == 1:
            return (
                3 * a[_intervals] * x2 + 2 * b[_intervals] * x + c[_intervals]
            )
        if deriv == 2:
            return 6 * a[_intervals] * x + 2 * b[_intervals]
        return None

    def translate_bspline_coeffs(self, bspline_coeffs, shift):
        """Somehow this doesn't work."""
        piecewise_translator = translate_cubic_pp_spline_coeffs(shift)

        transformation = (
            jnp.linalg.inv(B_TO_PP_SPLINE_BASIS_TRANSFORMATION)
            @ piecewise_translator
            @ B_TO_PP_SPLINE_BASIS_TRANSFORMATION
        )

        # trafo_eigen_vec = jnp.linalg.eig(transformation)[1][0,:]
        trafo_eigen_vec = transformation[1, :]
        trafo_eigen_vec = trafo_eigen_vec

        translated_bspline_coeffs = np.array(
            [
                np.convolve(bspline_coeffs[:, beam], trafo_eigen_vec, 'full')
                for beam in range(bspline_coeffs.shape[1])
            ]
        ).swapaxes(0, 1)
        # remove added dummy coeffs from convolution
        return translated_bspline_coeffs[1:-2]

    def evaluate_bspline_coeffs(self, bspline_coeffs, deriv_order=0):
        """Evaluate spline using the De Boor and the B-spline coefficients"""
        # Extract the relevant coefficients for each interval
        lower_indices = self.intervals - self.intpol_deg
        coeff_indices = lower_indices.reshape(-1, 1) + jnp.arange(
            self.intpol_deg + 1
        )
        coeff_subarrays = bspline_coeffs[coeff_indices]

        # Element-wise multiplication between coefficients and de_boor values
        # then sum over basis functions
        return jnp.einsum(
            'ijb,ji->ib', coeff_subarrays, self.de_boor_coeffs[deriv_order]
        )


def not_a_knot_rhs(values):
    return jnp.asarray(values)


# TODO: implement natural knot interpolator
## Natural knots spline boundary condition – currently unused –
def get_natural_knots(x, deg):
    return np.concatenate(
        [
            np.full(shape=(deg,), fill_value=x[0]),
            x,
            np.full(shape=(deg,), fill_value=x[-1]),
        ]
    )


def natural_derivative_bc(deg):
    _supported_degrees = (3, 5)

    if deg == 3:
        derivs_l_ord = np.array([2])
        derivs_l_val = np.array([0])

        derivs_r_ord = np.array([2])
        derivs_r_val = np.array([0])

    elif deg == 5:
        derivs_l_ord = np.array([3, 4])
        derivs_l_val = np.array([0, 0])

        derivs_r_ord = np.array([3, 4])
        derivs_r_val = np.array([0, 0])

    else:
        msg = f'unsupported degree {deg}'
        raise ValueError(msg)
    return derivs_l_ord, derivs_l_val, derivs_r_ord, derivs_r_val
