# -*- coding: utf-8 -*-
"""
Module errors of viperleed.tleedmlib.classes.fd_optimizer.

Created on 2023-12-19

@author: Alexander M. Imre (@amimre)

Defines the exceptions raised by the full dynamic calculation.
"""
class FullDynamicCalculationError(Exception):
    """Base class for exceptions raised by the full dynamic calculation."""
    def __init__(self, message):
        super().__init__(message)


class FullDynamicOptimizationOutOfBoundsError(FullDynamicCalculationError):
    """Raised when the optimization is out of bounds."""
    def __init__(self, message):
        super().__init__(message)


class FDParameterError(Exception):
    """Base class for exceptions raised by FDParameter."""
    def __init__(self, message):
        super().__init__(message)


class FDInvalidParameterError(Exception):
    """Raised when the user-requested parameter is invalid."""
    def __init__(self, message):
        super().__init__(message)


class FDParameterOutOfBoundsError(FDParameterError):
    """Raised when the user-requested bounds are impossible."""
    def __init__(self, message):
        super().__init__(message)

