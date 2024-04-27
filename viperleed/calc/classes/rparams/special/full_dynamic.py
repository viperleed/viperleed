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
from viperleed.tleedmlib.classes.fd_optimizer.fd_optimizers import (
    PARABOLA_SYNONYMS, ERROR_SYNONYMS, SCIPY_SYNONYMS)
from viperleed.tleedmlib.classes.fd_optimizer.fd_optimizers import (
    SingleParameterParabolaFit, SingleParameterBruteForceOptimizer,
    SingleParameterMinimizer, ParameterMinimizer
)


class FDSettings(SpecialParameter, param='FD'):

    def __init__(self, method=NO_VALUE):
        super().__init__()
        self.parameters = [] # list of FDParameter objects                    # TODO: Maybe this should use NO_VALUE?
        self.method = method

    def __bool__(self):
        return bool(self.parameters)

    @property
    def start_values(self):
        """Starting values for the parameters to be optimized."""
        return [param.start_value for param in self.parameters]

    @property
    def n_params(self):
        return len(self.parameters)

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
        self.parameters.append(fd_param)

    @property
    def user_defined_method(self):
        return self.method is not NO_VALUE

    def set_method(self, user_method, settings):
        """Called during parameter interpretation to store user defined method.

        We only store the method name and settings here. The actual
        interpretation happens later in the apply_method_settings method, 
        because it will depend of the number of parameters to be optimized."""
        self._method = user_method
        self.settings = settings

    def apply_method_settings(self):
        """Apply the user defined method and settings and """
        if self._method.lower() in PARABOLA_SYNONYMS:
            self._set_parabola_settings()
        elif self._method.lower() in ERROR_SYNONYMS:
            self._set_error_settings()
        elif self._method.lower() in SCIPY_SYNONYMS:
            self._set_scipy_minimizer_settings()
        else:
            raise NotImplementedError(
                f'Unknown optimization method for full-dynamic optimization: '
                f'{self._method}')

    def _set_parabola_settings(self):
        """Set settings for parabola method."""
        self.method = 'parabola'
        self.optimizer_class = SingleParameterParabolaFit
        if len(self.parameters) > 1:
            raise ValueError(
                'Parabola optimization can only be used for single parameter '
                'optimization. Parameters selected for optimization: '
                f'{self.parameters}.')

    def _set_error_settings(self):
            """Set settings for error method."""
            self.method = 'brute force'
            self.optimizer_class = SingleParameterBruteForceOptimizer
            # TODO: currently only single parameter optimization is supported
            if len(self.parameters) > 1:
                raise ValueError(
                    'Brute Force optimization and error estimation '
                    'can only be used for single parameter '
                    'optimization. Parameters selected for optimization: '
                    f'{self.parameters}.')

    def _set_scipy_minimizer_settings(self):
        """Set settings for the general one or multi parameter minimizer."""
        self.method = self._method.lower()  # which ever Scipy minimizer was used
        if len(self.parameters) == 1:
            self.optimizer_class = SingleParameterMinimizer
        else:
            self.optimizer_class = ParameterMinimizer

    def validate_optimizer_options(self, optimizer):
        """Check that the user defined optimizer options are valid.

        Must be run after the optimizer class is defined, i.e. after 
        apply_method_settings is called by parameters._checker."""
        # check that all required settings are defined
        for option in optimizer.required_settings:
            if option not in self.settings.keys():
                raise ValueError(
                    f'Minimizer type {self.method} requires the setting '
                    f'{option}.')
        # check that all settings are valid
        available_settings = {**optimizer.required_settings,
                              **optimizer.optional_settings}
        processed_settings = {}
        for setting, value in self.settings.items():
            if setting not in available_settings.keys():
                raise ValueError(
                    f'Unknown setting {setting} for method {self.method}.')
            # Type cast all the settings
            try:
                if available_settings[setting] is dict:  # needed for scipy options
                    processed_settings[setting] = eval(value)
                else:
                    processed_settings[setting] = available_settings[setting](value)
            except ValueError:
                raise ValueError(
                    f'Invalid value {value} for setting {setting} of FD_METHOD.')
        # set the settings
        self.settings = processed_settings
