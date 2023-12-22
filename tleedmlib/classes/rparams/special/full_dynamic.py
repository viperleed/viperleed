# -*- coding: utf-8 -*-
"""Module fd_optimizers of viperleed.tleedmlib.classes.rparams.special.

Created on 2023-12-19

@author: Alexander M. Imre (@amimre)
"""
from dataclasses import dataclass

from ._base import SpecialParameter
from .._defaults import NO_VALUE

from viperleed.tleedmlib.classes.fd_optimizer import FDParameter
from viperleed.tleedmlib.classes.fd_optimizer.fd_parameter import FD_PARAMETERS, AVAILABLE_FD_PARAMETERS
from viperleed.tleedmlib.classes.fd_optimizer.errors import (
    FDParameterOutOfBoundsError, FDInvalidParameterError)


class FDSettings(SpecialParameter, param='FD'):

    def __init__(self, method=NO_VALUE):
        super().__init__()
        self.parameters = set()
        self.method = method


    def __post_init__(self):
        pass

    def __bool__(self):
        return bool(self.parameters)

    @property
    def n_params(self):
        return len(self.what)

    @property
    def what(self):
        if self._what is NO_VALUE:
            return set()
        return self._what

    def add_fd_parameter(self, rp, parameter, bounds):
        """Add a parameter to the list of parameters to be optimized."""

        if not isinstance(parameter, str):
            raise ValueError(f'parameter must be a string.')
        # check that parameter is a valid parameter
        if parameter.lower() not in AVAILABLE_FD_PARAMETERS:
            raise FDInvalidParameterError
        (f'Parameter {parameter} is not a valid parameter '
         'for full-dynamic optimization.')
        if len(bounds) == 2:
            _min, _max, _step = float(bounds[0]), float(bounds[1]), None
        elif len(bounds) == 3:
            _min, _max, _step = [float(val) for val in bounds]
        else:
            raise ValueError('Bounds must be a tuple of length 2 or 3.')
        if _min > _max:
            _min, _max = _max, _min
        # check user bounds vs. maximum possible bounds
        if _min < FD_PARAMETERS[parameter.lower()]['max_bounds'][0]:
            raise FDParameterOutOfBoundsError(
                f'Lower bound {_min} for {parameter} is too small.')
        if _max > FD_PARAMETERS[parameter.lower()]['max_bounds'][1]:
                raise FDParameterOutOfBoundsError(
                f'Upper bound {_max} for {parameter} is too large.')
        fd_param = FDParameter(
            name=parameter.lower(),
            start_value=FD_PARAMETERS[parameter.lower()]['x0'](rp),
            bounds=(_min, _max),
            transform_func=FD_PARAMETERS[parameter.lower()]['eval'],
        )
        self.parameters.add(fd_param)


    @property
    def user_defined_method(self):
        return self.method is not NO_VALUE

    def set_method(self, method, settings):
        pass
