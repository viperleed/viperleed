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
