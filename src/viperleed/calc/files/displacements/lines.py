"""Module lines."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__created__ = '2024-10-04'

from abc import ABC, abstractmethod
from collections import namedtuple

from viperleed_jax.perturbation_type import PerturbationType

from .direction import Direction
from .range import DisplacementsRange
from .targeting import BSTarget

LoopMarkerLine = namedtuple('LoopMarkerLine', ['type'])
SearchHeaderLine = namedtuple('SearchHeaderLine', ['label'])
SectionHeaderLine = namedtuple('SectionHeaderLine', ['section'])

_BELOW_DEBUG = 2


class ParsedLine(ABC):
    """Base class for parsing of non-header lines in the DISPLACEMENTS file.

    The class is used to parse lines in the displacements file. It can be
    initialized with a line string and will parse it into its components.
    """

    def __init__(self, line):
        """Perform basic parsing and checking of the line format."""
        # strip surplus whitespace
        self._raw_line = ' '.join(line.strip().split())

        # if not exactly one '=' is present, raise an error
        if self._raw_line.count('=') != 1:
            msg = (
                f'Invalid DISPLACEMENTS line format: "{self._raw_line}". '
                'Expected format: "<labels> = <values>".'
            )
            raise ValueError(msg)

        # split the line into raw left and right hand sides
        self._lhs, self._rhs = self._raw_line.split('=')


    @abstractmethod
    def __repr__(self):
        """Return the string representation of the line."""

class GeoDeltaLine(ParsedLine):
    """Class to parse lines in the GEO_DELTA block of DISPLACEMENTS.

    Lines in the GEO_DELTA block are of the form:
        <target> [, <target>] <direction> = <range>
    where <target>, <direction>, and <range> are tokes that are parsed by the
    `BSTarget`, `Direction`, and `DisplacementsRange` classes, respectively.
    """

    def __init__(self, line):
        super().__init__(line)

        # Left hand side
        # split off last element of the left hand side
        lhs_parts = self._lhs.split()
        # check if the last part is a direction
        try:
            self.direction = Direction(lhs_parts[-1])
        except ValueError as err:
            msg = ('Unable to parse direction information from line in '
                   f'GEO_DELTA block: {self.line}')
            raise ValueError(msg) from err

        # parse the rest of the left hand side into targets
        # TODO

        # parse right hand side to a range
        self.range = DisplacementsRange(self._rhs)


    def __repr__(self):
        """Return the string representation of the line."""
        if self._line is None:
            line = f'{self.label} {self.which}'
            if self.direction is not None:
                line += f' {self.direction}'
            line += f' = {self.range.start}'
            if self.range.stop is not None:
                line += f' {self.range.stop}'
            if self.range.step is not None:
                line += f' {self.range.step}'
        else:
            line = self._line
        return line


def _get_target(label, which):
    if which is None:
        return BSTarget(label)
    return BSTarget(f'{label} {which}')


class GeoDeltaLine:
    def __init__(self, label, which, direction, start, stop, step,
                 line=None):
        self._line= line
        self.label = label
        self.which = which
        self.targets = _get_target(label, which)
        self.direction = Direction(direction)
        self.range = DisplacementsRange(start, stop, step)

    def __eq__(self, other):
        if isinstance(other, GeoDeltaLine):
            return (
                self.targets == other.targets
                and self.direction == other.direction
                and self.range == other.range
            )
        return False

    def __repr__(self):
        """Return the string representation of the line."""
        if self._line is None:
            line = f'{self.label} {self.which}'
            if self.direction is not None:
                line += f' {self.direction}'
            line += f' = {self.range.start}'
            if self.range.stop is not None:
                line += f' {self.range.stop}'
            if self.range.step is not None:
                line += f' {self.range.step}'
        else:
            line = self._line
        return line


class VibDeltaLine:
    def __init__(self, label, which, start, stop, step, line=None):
        self._line = line
        self.label = label
        self.which = which
        self.targets = _get_target(label, which)
        self.range = DisplacementsRange(start, stop, step)

    def __eq__(self, other):
        if isinstance(other, VibDeltaLine):
            return (
                self.targets == other.targets
                and self.range == other.range
            )
        return False

    def __repr__(self):
        """Return the string representation of the line."""
        if self._line is None:
            line = f'{self.label} {self.which}'
            line += f' = {self.range.start}'
            if self.range.stop is not None:
                line += f' {self.range.stop}'
            if self.range.step is not None:
                line += f' {self.range.step}'
        else:
            line = self._line
        return line


class OccDeltaLine:
    def __init__(self, label, which, chem_blocks, line=None):
        self._line = line
        self.label = label
        self.which = which
        self.targets = _get_target(label, which)
        self.chem_blocks = chem_blocks

    def __eq__(self, other):
        if isinstance(other, OccDeltaLine):
            return (
                self.targets == other.targets
                and self.chem_blocks == other.chem_blocks
            )
        return False

    def __repr__(self):
        """Return the string representation of the line."""
        if self._line is None:
            line = f'{self.label} {self.which}'
            line += f' = {self.chem_blocks}'
        else:
            line = self._line
        return line


class ConstraintLine:
    def __init__(self, constraint_type, targets, direction, value, line=None):
        self._line = line
        self.constraint_type = PerturbationType.from_string(constraint_type)
        self.direction = None
        if self.constraint_type == PerturbationType.GEO and direction is not None:
            self.direction = Direction(direction)
        elif self.constraint_type != PerturbationType.GEO and direction is not None:
            msg = (
                f'A direction is not allowed for {self.constraint_type} '
                f'constraints.'
            )
            raise ValueError(msg)
        self.targets = BSTarget(targets)
        self.value = value

    def __eq__(self, other):
        if isinstance(other, ConstraintLine):
            return (
                self.constraint_type == other.constraint_type
                and self.targets == other.targets
                and self.value == other.value
                and (self.direction == other.direction
                     if self.direction is not None else other.direction is None)
            )
        return False

    def __repr__(self):
        """Return the string representation of the line."""
        if self._line is None:
            line = f'{self.constraint_type} {self.targets}'
            if self.direction is not None:
                line += f' {self.direction}'
            line += f' = {self.value}'
        else:
            line = self._line
        return line


class OffsetsLine:
    def __init__(self, offset_type, targets, direction, value, line=None):
        self._line= line
        self.offset_type = offset_type
        self.targets = BSTarget(targets)
        self.value = value
        if self.offset_type == 'geo':
            if direction is not None:
                self.direction = Direction(direction)
            else:
                self.direction = None
        elif self.offset_type != 'geo' and direction is not None:
            msg = (
                f'A direction is not allowed for {self.offset_type} '
                f'offsets.'
            )
            raise ValueError(msg)

    def __eq__(self, other):
        if isinstance(other, ConstraintLine):
            return (
                self.offset_type == other.constraint_type
                and self.parameters == other.parameters
                and self.value == other.value
            )
        return False

    def __repr__(self):
        """Return the string representation of the line."""
        if self._line is None:
            line = f'{self.offset_type} {self.targets}'
            if self.direction is not None:
                line += f' {self.direction}'
            line += f' = {self.value}'
        else:
            line = self._line
        return line
