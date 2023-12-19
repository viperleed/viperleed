# -*- coding: utf-8 -*-
"""
Module fd_parameter of viperleed.tleedmlib.classes.fd_optimizer.

Created on 2023-12-19

@author: Alexander M. Imre (@amimre)

Defines the FDParameter class, which is used to store information about
parameters that are optimized in the full-dynamic calculation.
"""
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
