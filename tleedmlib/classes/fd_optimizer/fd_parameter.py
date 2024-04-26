# -*- coding: utf-8 -*-
"""
Module fd_parameter of viperleed.tleedmlib.classes.fd_optimizer.

Created on 2023-12-19

@author: Alexander M. Imre (@amimre)

Defines the FDParameter class, which is used to store information about
parameters that are optimized in the full-dynamic calculation.
"""
import numpy as np


class FDParameter():
    """Base class for parameters accessible via full dynamic calculation."""

    def __init__(self, name, start_value, bounds, transform_func):
        self.name = name
        self.start_value = start_value
        # make sure start value is within bounds
        if start_value < bounds[0] or start_value > bounds[1]:
            raise ValueError(f"Start value for {start_value} is not within "
                             "bounds {bounds}")
        self.bounds = bounds

        self.transform_func = transform_func

    def apply(self, rparams, slab, val):
        """Apply the parameter to the given rparams and slab."""
        self.transform_func(rparams, slab, val)


# parameters accessible in the full dynamic optimization
# must specify bounds, a function for altering the parameter and an initial 
# guess x0
FD_PARAMETERS = {
    'v0i': {
        'max_bounds': (0, 15),
        'eval': lambda r, s, v: setattr(r, "V0_IMAG", v),
        'x0': lambda rp: rp.V0_IMAG,
    },
    'theta': {
        'max_bounds': (-20, 20),
        'eval': lambda r, s, v: setattr(r, "THETA", v),
        'x0': lambda rp: rp.THETA,
    },
    'phi': {
        'max_bounds': (0, 360),
        'eval': lambda r, s, v: setattr(r, "PHI", v),
        'x0': lambda rp: rp.PHI,
    }
}

for scaling in ('a', 'b', 'c', 'ab', 'bc', 'abc'): # scaling of lattice vectors
    FD_PARAMETERS[scaling] = {
        'max_bounds': (0.1, 10),
        'eval': lambda r, s, v: apply_scaling(s, r, scaling, v),
        'x0': lambda rp: 1,
    }


# parameters that can be optimized in FD optimization
AVAILABLE_FD_PARAMETERS = tuple(FD_PARAMETERS.keys())

def apply_scaling(sl, rp, which, scale):
    m = np.eye(3)
    if "a" in which:
        m[0, 0] *= scale
    if "b" in which:
        m[1, 1] *= scale
    if "c" in which:
        m[2, 2] *= scale
    sl.update_fractional_from_cartesian()
    sl.ucell = np.dot(sl.ucell, m)
    sl.update_cartesian_from_fractional(updateOrigin=True)
    sl.bulkslab.update_fractional_from_cartesian()
    sl.bulkslab.ucell = np.dot(sl.bulkslab.ucell, m)
    sl.bulkslab.update_cartesian_from_fractional()
    if isinstance(rp.BULK_REPEAT, (np.floating, float)):
        rp.BULK_REPEAT *= scale
    elif rp.BULK_REPEAT is not None:
        rp.BULK_REPEAT = np.dot(rp.BULK_REPEAT, m)
