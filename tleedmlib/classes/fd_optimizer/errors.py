class FullDynamicCalculationError(Exception):
    """Base class for exceptions raised by the full dynamic calculation."""

    def __init__(self, message):
        super().__init__(message)


class FullDynamicOptimizationOutOfBoundsError(FullDynamicCalculationError):
    """Raised when the optimization is out of bounds."""
    def __init__(self, message):
        super().__init__(message)
