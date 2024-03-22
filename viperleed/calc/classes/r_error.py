"""Module r_error of viperleed.calc.classes.

Defines data structures and functions for storing and manipulating
errors curves.
"""

__authors__ = (
    'Alexander M. Imre (@amimre)',
    'Florian Kraushofer (@fkraushofer)',
    )
__copyright__ = 'Copyright (c) 2019-2024 ViPErLEED developers'
__created__ = '2023-01-12'
__license__ = 'GPLv3+'

import numpy as np
import copy
import logging

from viperleed.calc.lib.base import range_to_str


logger = logging.getLogger(__name__)


class R_Error():
    """Data structure for storing errors after they have been calculated.

    Stores the displacements that led to the R-factors even if they are
    modified later. Displacements are stored for one atom only, the other are
    linked.
    """

    def __init__(self, r_type, atoms, mode, rfacs, disp_label, lin_disp,
                 v0i=None, energy_range=None):
        self.r_type = r_type # Type of R-factor (int). Only Pendry R-factor (r_type=1) can give statistical error estimates.
        self.atoms = atoms  # atoms that have been varied together
        self.mode = mode    # vib, geo, or occ
        self.rfacs = rfacs  # the r-factors from the variations
        self.var_r = None
        self.displacements = []   # displacements of atoms[0]
        self.disp_label = disp_label
        self.lin_disp = lin_disp  # linearized displacement
        self.main_element = ""    # element occupation displayed in output
        self.elem_occ = None
        d = {}
        at = atoms[0]   # store displacements for this one; main element
        if mode == "occ":
            if atoms[0].el in atoms[0].disp_occ:
                self.main_element = atoms[0].el  # main site element
            else:
                # element with highest occupation in refcalc
                self.main_element = max(atoms[0].site.occ,
                                        key=lambda k: atoms[0].site.occ[k])
            self.displacements = copy.deepcopy(atoms[0].disp_occ[
                                                        self.main_element])
            self.lin_disp = copy.deepcopy(self.displacements)
            self.elem_occ = copy.deepcopy(atoms[0].disp_occ)
            vac_occ = np.ones_like(self.lin_disp)
            for el_occ in self.elem_occ.values():
                vac_occ -= el_occ
            if np.sum(vac_occ) > 1e-4:  # If vacancies, put them into output CSV
                self.elem_occ["Vac"] = vac_occ
        else:
            if mode == "geo":
                d = atoms[0].disp_geo
            elif mode == "vib":
                d = atoms[0].disp_vib
        if len(d) > 0:
            if at.el in d:
                k = at.el
            else:
                k = "all"
            self.displacements = copy.deepcopy(d[k])

        if v0i and energy_range and r_type==1:
            self.calc_var_r(v0i, energy_range)


    def calc_var_r(self, v0i, energy_range):
        r"""Calculates the variance of the Pendry R-factor.

        The variance of the Pendry R-factor can be estimated from the
        imaginary part of the inner potential and the measurement energy
        range as:

        .. math::
            \mathrm{var}(R_p) = \sqrt{ \frac{8V_{0\mathrm{i}}}{E_{\mathrm{range}}}}R_{\mathrm{min}}

        Parameters
        ----------
        v0i : float
            Imaginary part of the inner potential (eV).
        energy_range : float
            Energy range in eV.

        Raises
        ------
        ValueError:
            If v0i or energy_range are not float values.
        ValueError:
            If v0i or energy_range are negative. This would be
            unphysical.
        """
        if self.r_type != 1:
            raise ValueError("Var(R) can only be calculated when using the "
                             "Pendry R-factor.")
        if not isinstance(v0i, float) and isinstance(energy_range, float):
            raise ValueError("v0i and energy_range must be float values.")
        if v0i<= 0 or energy_range <= 0:
            raise ValueError("v0i and energy range must be positive ")
        if not self.rfacs:
            raise ValueError("Cannot calculate variance of R-factor without "
                             "stored R-factors.")
        self.var_r = np.sqrt(8*np.abs(v0i) / energy_range) * self.get_r_min


    @property
    def get_error_estimates(self):
        """Decide if statistical error estimates can be given and
        calculates them.

        Returns
        -------
        tuple of float or None
            A tuple containing upper and lower error bounds, es
            estimated from the variance of R_P, if errors can be
            estimated. If either bound can not be estimates (e.g. no
            intersection with line R_min + var(R)), the element is None.
            If neither bound can be estimated, both values are None.

        Raises
        ------
        ValueError
            If self.lin_disp is not initialized.
        """

        # error estimates are only available for Pendry R-factor
        if self.r_type != 1:  #TODO: when we introduce a new R-factor, make sure to update this.
            return (None, None)
        r_min = self.get_r_min
        if self.var_r is None:
            logger.warning("Cannot calculate statistical errors for "
                           f'atoms {range_to_str([at.num for at in self.atoms])}')
            return (None, None)

        if self.lin_disp is None:
            raise ValueError("Linear displacements not initialized.")

        # more than 4 points
        if len(self.rfacs) < 3 or not self.var_r:
            return (None, None)

        # minimum R factor of this error is not necessarily the same as
        # the overall minimum R factor
        err_min_R = min(self.rfacs)

        # TODO change!!!!!
        if err_min_R > r_min + self.var_r:
            return (None, None)

        # treat upper and lower half of error curve
        err_min_R_id = self.rfacs.index(err_min_R)
        best_param = self.lin_disp[self.rfacs.index(err_min_R)]
        left_rfacs = self.rfacs[:err_min_R_id+1]
        left_lin_disp = self.lin_disp[:err_min_R_id+1]
        right_rfacs = self.rfacs[err_min_R_id:]
        right_lin_disp = self.lin_disp[err_min_R_id:]

        error_estimates = []
        for rfacs, lin_disp in ((left_rfacs, left_lin_disp),
                                (right_rfacs, right_lin_disp)):
            # find zero crossing
            offset_rfacs = rfacs - (r_min + self.var_r)
            if (len(rfacs) < 3 or                                         # too few points
                get_n_zero_crossings(offset_rfacs) != 1): # unclear crossing of line
                error_estimates.append(None)
            elif(max(rfacs) <= r_min + self.var_r):       # maximum R below r_min + var(R) line
                error_estimates.append(None)
            else:
                crossing_point = get_zero_crossing(lin_disp, offset_rfacs)
                delta = abs(best_param - crossing_point)
                error_estimates.append(delta)
        return tuple(error_estimates)

    @property
    def get_r_min(self):
        """Returns minimum R-factor of error.

        Returns
        -------
        float
            Minimum R-factor stored in self.rfacs.
        """
        return min(self.rfacs)

    @property
    def get_p_min(self):
        """Returns parameter value at minimum R-factor.

        Returns
        -------
        float
            parameter value at minium R-factor (self.rfacs).
        """
        return self.lin_disp[self.rfacs.index(self.get_r_min)]


def get_n_zero_crossings(x_arr):
    """Counts and returns the number of zero crossings in input.

    x_arr must be discrete numerical data. This function counts the
    number of times the zero line is crossed when going through the
    input values in order, e.g. x_arr=[-3, 1 , 1, -5, -2] would return
    2.

    Parameters
    ----------
    x_arr : list or ndarray of numbers
        Discrete numerical input data.

    Returns
    -------
    int
        Number of times the discrete input data crossed the zero line.
    """
    arr = np.array(x_arr)
    n_zero_crossings = int((1 - np.sign((arr[:-1]*arr[1:]))).sum()/2)
    return n_zero_crossings


def get_zero_crossing(x_arr, y_arr, eps=1e-6):
    """Computes and returns the x_value where the provided dataset crosses zero.

    x_arr and y_arr are sets of x and y values respectively. This
    function first checks in between which points the zero line is
    crossed. It then performs a linear interpolation between those two
    points and gives the interpolated x values, where the zero line is
    crossed. If two subsequent points are zero, the middle is returned.
    If the dataset crosses the zero line multiple times, only the first
    x value will be returned. See also get_n_zero_crossings().

    Parameters
    ----------
    x_arr : list or array of numbers.
        x_values of the data to be checked.
        Must in ascending order.
    y_arr : _type_
        _description_
    eps: float
        Numerical accuracy of comparison.

    Returns
    -------
    float or None
        The linearly interpolated x value of the zero crossing.
        None if no zero crossing is found.

    Raises
    ------
    ValueError
        if x_arr and y_arr have incompatible shapes.
    ValueError
        if the values of x_arr are not monotonically increasing.
    """

    if len(x_arr) != len(y_arr):
        raise ValueError("Different sized arrays supplied.")
    # make sure x_arr is monotonically increasing
    if not np.all(x_arr[1:] >= x_arr[:-1]):
        raise ValueError("Values in x_arr are not monotonically increasing")

    # find location of zero crossing
    for id in range(len(x_arr)-1):
        if y_arr[id]*y_arr[id+1] <= 0:
            R_jump = abs(y_arr[id+1]-y_arr[id])
            if R_jump < eps:  # avoid division by zero
                return (x_arr[id]+x_arr[id+1])/2
            return x_arr[id]-y_arr[id]*(x_arr[id+1]-x_arr[id])/(y_arr[id+1]-y_arr[id])
    # return None if no zero crossing found
    return None
