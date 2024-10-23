"""Module errors."""

__authors__ = ("Alexander M. Imre (@amimre)",)
__created__ = "2024-10-04"

class InvalidDisplacementsSyntaxError(ValueError):
    pass

class OffsetsNotAtBeginningError(InvalidDisplacementsSyntaxError):
    pass


class InvalidSearchBlocksError(InvalidDisplacementsSyntaxError):
    pass

class InvalidSearchLoopError(InvalidDisplacementsSyntaxError):
    pass


class SymmetryViolationError(ValueError):
    pass
