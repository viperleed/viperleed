"""Module for the <range> token in the DISPLACEMENTS file."""

__authors__ = ("Alexander M. Imre (@amimre)",)
__copyright__ = "Copyright (c) 2019-2025 ViPErLEED developers"
__created__ = "2024-10-15"
__license__ = "GPLv3+"

from .base import DisplacementsFileToken, TokenParserError
from viperleed.calc.files.new_displacements import DISPLACEMENTS_FILE_EPS


class RangeTokenParserError(TokenParserError):
    """Class for parsing Errors in the RangeToken."""


class RangeToken(DisplacementsFileToken):
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

    def __init__(self, range_str: str):
        """Construct a RangeToken from a string."""
        # TODO: TensErLEED compatibility: backend requires the step to be specified
        parts = range_str.strip().split()
        if len(parts) not in (2, 3):
            msg = (
                f'Invalid range format: "{range_str}". Expected format: '
                '"<start> <stop> [<step>]".'
            )
            raise RangeTokenParserError(msg)

        self.start, sef.stop, self.step = None
        attrs = ('start', 'stop', 'step')
        for part, attr in zip(parts, attrs):
            try:
                value = float(part)
            except ValueError as err:
                msg = f'Non-numeric value in range: "{range_str}"'
                raise RangeTokenParserError(msg) from err
            setattr(self, attr, value)

    @property
    def has_step(self):
        """Check if the range has a step defined."""
        return self.step is not None
    @classmethod
    def from_floats(
        cls, start: float, stop: float, step=None
    ) -> 'RangeToken':
        """Alternate constructor using numeric values directly."""
        inst = cls.__new__(cls)
        inst.start = start
        inst.stop = stop
        inst.step = step
        return inst

    def __eq__(self, other):
        """Compare two RangeToken objects for equality."""
        if not isinstance(other, RangeToken):
            return NotImplemented
        if self.has_step != other.has_step:
            return False
        if self.has_step:
            return (
                abs(self.start - other.start) < DISPLACEMENTS_FILE_EPS
                and abs(self.stop - other.stop) < DISPLACEMENTS_FILE_EPS
                and abs(self.step - other.step) < DISPLACEMENTS_FILE_EPS
            )

        return (
            abs(self.start - other.start) < DISPLACEMENTS_FILE_EPS
            and abs(self.stop - other.stop) < DISPLACEMENTS_FILE_EPS
        )

    def __str__(self):
        """Return a string representation of the RangeToken object."""
        if self.has_step:
            return (f'RangeToken(start={self.start}, stop={self.stop}, '
                    f'step={self.step})')
        return f'RangeToken(start={self.start}, stop={self.stop})'
