"""Module errors."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__created__ = '2024-10-04'


class DisplacementsError(ValueError):
    """Base class for all displacement-related errors."""

    def __init__(self, message):
        super().__init__(message)
        self.message = message

class InvalidDisplacementsFileError(DisplacementsError):
    """Base class for all errors related to the DISPLACEMENTS file."""

    def __init__(self, message):
        super().__init__(message)
        self.message = message

class InvalidDisplacementsSyntaxError(DisplacementsError):
    """Base class for all errors related to the DISPLACEMENTS file."""

    def __init__(self, message):
        super().__init__(message)
        self.message = message

class IncompatibleBackendError(DisplacementsError):
    """Error raised when the chosen backend is incompatible with the file."""

    def __init__(self, message):
        super().__init__(message)
        self.message = message

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
