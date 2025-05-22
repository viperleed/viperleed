"""Module errors of guilib.leedsim.classes.woods.

======================================
  ViPErLEED Graphical User Interface
======================================
 *** package guilib.leedsim.classes.woods ***

Author: Michele Riva (@michele-riva)
Created: 2024-03-01

Defines exceptions raised by the Woods class. Refactored from the
former woods.py module.
"""

from viperleed.guilib.helpers import array_to_string


class WoodsError(Exception):
    """Base exception for Wood's-related error."""


class MatrixIncommensurateError(WoodsError):
    """Matrix is incommensurate."""

    def __init__(self, matrix, message=''):
        self.matrix = array_to_string(matrix)
        if not message:
            message = (f'Matrix {self.matrix} gives an incommensurate lattice,'
                       ' i.e., it is singular or has non-integer elements')
        super().__init__(message)


class WoodsInvalidForBasisError(WoodsError):
    """A string woods is not appropriate for the current bulk basis."""


class WoodsNotRepresentableError(WoodsError):
    """Matrix is not Wood's-representable."""

    def __init__(self, matrix, message=''):
        self.matrix = array_to_string(matrix)
        if not message:
            message = f'Matrix {self.matrix} is not Woods-representable'
        super().__init__(message)


class WoodsSyntaxError(WoodsError):
    """A Wood's string has invalid syntax."""
