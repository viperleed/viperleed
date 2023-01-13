# -*- coding: utf-8 -*-
"""
Created on Jan 12 2023

@author: Alexander M. Imre, Florian Kraushofer

Data structure for storing and manipulating errors curves for the error
calc after they have been calculated.
"""

import numpy as np
import copy
import logging

from viperleed.tleedmlib.base import range_to_str

logger = logging.getLogger("tleedm.classes.r_error")


class R_Error():
    """
    Data structure for storing errors after they have been calculated.
    Stores the displacements that led to the R-factors even if they are
    modified later. Displacements are stored for one atom only, the other are
    linked.
    """

    def __init__(self, atoms, mode, rfacs, disp_label, lin_disp, 
                 v0i=None, energy_range=None):
        self.atoms = atoms  # atoms that have been varied together
        self.mode = mode    # vib, geo, or occ
        self.rfacs = rfacs  # the r-factors from the variations
        self.displacements = []   # displacements of atoms[0]#
        self.disp_label = disp_label
        self.lin_disp = lin_disp  # linearized displacement
        self.main_element = ""    # element occupation displayed in output
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
            
        if v0i and energy_range:
            self.calc_var_r(v0i, energy_range)


    def calc_var_r(self, v0i, energy_range):
        self.var_r = np.sqrt(8*np.abs(v0i) / energy_range) * self.get_r_min


    @property
    def get_error_estimates(self):
    # decide if estimate errors

        r_min = self.get_r_min
        if self.var_r is None:
            logger.warning("Cannot calculate statistical errors for "
                           f'atoms {range_to_str([at.oriN for at in self.atoms])}')
            return (None, None)

        if self.lin_disp is None:
            raise ValueError("Linear displacements not initialized.")

        # more than 4 points
        if len(self.rfacs) < 3 or not self.var_r:
            return (None, None)

        # minimum R factor of this error is not necessarily the same as
        # the overall minimum R factor
        err_min_R = min(self.rfacs)

        if (err_min_R > r_min + self.var_r or              # minimum R is below line
            self.rfacs[0] < r_min + self.var_r or        # left max R is below line
            self.rfacs[-1] < r_min + self.var_r):        # right max R is below line
            return (None, None)

        # treat upper and lower half of error curve
        err_min_R_id = self.rfacs.index(err_min_R)
        left_rfacs = self.rfacs[:err_min_R_id+1]
        left_lin_disp = self.lin_disp[:err_min_R_id+1]
        right_rfacs = self.rfacs[err_min_R_id:]
        right_lin_disp = self.lin_disp[err_min_R_id:]

        error_estimates = []
        for lin_disp, rfacs in (zip(left_rfacs, left_lin_disp),
                                zip(right_rfacs, right_lin_disp)):
            # find zero crossing 
            if (len(rfacs) < 3 or                                    # too few points
                get_n_zero_crossings(rfacs - (r_min + self.var_r)) != 1): # unclear crossing of line
                error_estimates.append(None)
                continue
            error_estimates.append(get_zero_crossing(lin_disp, rfacs))
        return tuple(error_estimates)

    @property
    def get_r_min(self):
        return min(self.rfacs)

    @property
    def get_p_min(self):
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
    n_zero_crossings = int((1 - np.sign((arr[:-1]*arr[1:]))).sum())
    return n_zero_crossings


def get_zero_crossing(x_arr, y_arr):
    """Computes and returns the x_value where the provided dataset
    crosses the zero line.
    
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

    eps = 1e-6
    if len(x_arr) != len(y_arr):
        raise ValueError("Different sized arrays supplied.")
    # make sure x_arr is monotonically increasing
    if not np.all(x_arr[:, 1:] >= x_arr[:, :-1], axis=1):
        raise ValueError("Values in x_arr are not monotonically increasing")

    # find location of zero crossing
    for id in range(len(x_arr)):
        if y_arr[id]*y_arr[id+1] <= 0:
            R_jump = y_arr[id+1]-y_arr[id]
            if R_jump < eps:  # avoid division by zero
                return (x_arr[id]+x_arr[id+1])/2
            return -y_arr[id]*(x_arr[id+1]-x_arr[id])/(y_arr[id+1]-y_arr[id])
    # return None if no zero crossing found
    return None
