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

#TODO: add more exceptions
