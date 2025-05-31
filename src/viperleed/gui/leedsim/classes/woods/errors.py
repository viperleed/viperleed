"""Module errors of viperleed.gui.leedsim.classes.woods.

Defines exceptions raised by the Woods class. Refactored from the
former woods.py module.
"""

__authors__ = (
    'Michele Riva (@michele-riva)',
    )
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2024-03-01'
__license__ = 'GPLv3+'

from viperleed.gui.helpers import array_to_string


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
