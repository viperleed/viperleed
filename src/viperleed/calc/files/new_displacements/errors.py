"""Module errors of viperleed.calc.files.new_displacements."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2024-10-04'
__license__ = 'GPLv3+'


class DisplacementsFileError(Exception):
    """Base class for all displacement-related errors."""

    def __init__(self, message):
        super().__init__(message)
        self.message = message


class AlreadyReadError(DisplacementsFileError):
    """Error raised when the file has already been read."""


class DisplacementsSyntaxError(DisplacementsFileError):
    """Base class for all errors related to the DISPLACEMENTS file."""


class OffsetsNotAtBeginningError(DisplacementsFileError):
    """Error raised when the OFFSETS block is not at the start of the file."""


class InvalidSearchBlocksError(DisplacementsFileError):
    """Error raised when the search blocks are not valid."""


class InvalidSearchLoopError(DisplacementsFileError):
    """Error raised when the search loop is not valid."""


class SymmetryViolationError(DisplacementsFileError):
    """Error raised when the requested displacements violate symmetry."""


class UnknownDisplacementsSegmentError(DisplacementsFileError):
    """Exception raised when a line does not match any known segment header."""
