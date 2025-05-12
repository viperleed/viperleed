"""Module range."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__created__ = '2024-10-15'


class RangeToken:
    """Class to parse and represent displacement ranges.

    Ranges are specified in the form of strings like:
        <start> <stop> [<step>]
    where the step is optional (backwards compatibility with TensErLEED). The
    class can also be initialized from numeric values directly using the
    `from_floats` class method.

    Parameters
    ----------
    range_str : str
        The range string to parse.
    """

    _EPS = 1e-6

    def __init__(self, range_str: str):
        """Construct a DisplacementsRange from a string."""
        parts = range_str.strip().split()
        if len(parts) < 2 or len(parts) > 3:
            msg = (
                f'Invalid range format: "{range_str}". Expected format: '
                '"<start> <stop> [<step>]".'
            )
            raise ValueError(msg)

        try:
            start = float(parts[0])
            stop = float(parts[1])
            step = float(parts[2]) if len(parts) == 3 else None
        except ValueError as err:
            msg = f'Non-numeric value in range: "{range_str}"'
            raise ValueError(msg) from err

        # check that step is valid
        _check_step(step)

        self.start = start
        self.stop = stop
        self.step = step
        self.has_step = step is not None

    @classmethod
    def from_floats(
        cls, start: float, stop: float, step=None
    ) -> 'RangeToken':
        """Alternate constructor using numeric values directly."""
        _check_step(step)
        inst = cls.__new__(cls)
        inst.start = start
        inst.stop = stop
        inst.step = step
        inst.has_step = step is not None
        return inst

    def __eq__(self, other):
        """Compare two DisplacementsRange objects for equality."""
        if not isinstance(other, RangeToken):
            return False
        if self.has_step != other.has_step:
            return False
        if self.has_step:
            return (
                abs(self.start - other.start) < self._EPS
                and abs(self.stop - other.stop) < self._EPS
                and abs(self.step - other.step) < self._EPS
            )
        if other.step is not None:
            return False
        return (
            abs(self.start - other.start) < self._EPS
            and abs(self.stop - other.stop) < self._EPS
        )

    def __repr__(self):
        """Return a string representation of the DisplacementsRange object."""
        if self.has_step:
            return (f'DisplacementsRange(start={self.start}, stop={self.stop}, '
                    f'step={self.step})')
        return f'DisplacementsRange(start={self.start}, stop={self.stop})'

def _check_step(step):
    """Check if the step is positive."""
    if step is not None and step <= 0:
        msg = f'Step must be positive: "{step}"'
        raise ValueError(msg)
