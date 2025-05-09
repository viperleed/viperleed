"""Module lines."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__created__ = '2024-10-04'

import re
from abc import ABC, abstractmethod
from collections import namedtuple

from viperleed_jax.files.displacements.regex import DIRECTION_PATTERN
from viperleed_jax.perturbation_type import PerturbationType

from .direction import Direction
from .errors import InvalidDisplacementsSyntaxError
from .range import DisplacementsRange
from .targeting import Targets

LoopMarkerLine = namedtuple('LoopMarkerLine', ['type'])
SearchHeaderLine = namedtuple('SearchHeaderLine', ['label'])
SectionHeaderLine = namedtuple('SectionHeaderLine', ['section'])

_BELOW_DEBUG = 2

# precompiled direction-at-end regex for separating out the direction token
_DIR_AT_END = re.compile(
    rf'(?P<dir>{DIRECTION_PATTERN})\s*$'
)

class ParsedLine(ABC):
    """Base class for parsing of non-header lines in the DISPLACEMENTS file.

    The class is used to parse lines in the displacements file. It can be
    initialized with a line string and will parse it into its components.
    """

    def __init__(self, line: str):
        """Perform basic parsing and checking of the line format."""
        # strip surplus whitespace
        self.raw_line = ' '.join(line.strip().split())

        # if not exactly one '=' is present, raise an error
        if self.raw_line.count('=') != 1:
            msg = (
                f'Invalid DISPLACEMENTS line format: "{self.raw_line}". '
                'Expected format: "<labels> = <values>".'
            )
            raise InvalidDisplacementsSyntaxError(msg)

        # split the line into raw left and right hand sides
        self._lhs, self._rhs = self.raw_line.split('=')


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

    def __init__(self, line: str):
        super().__init__(line)

        # Left hand side
        # split off last element of the left hand side
        targets_str, dir_str = separate_direction_from_targets(self._lhs)
        if not targets_str or not dir_str:
            msg = (
                f'Invalid GEO_DELTA line format: "{self._lhs}". '
                'Expected format: "<targets> <direction> = <range>".'
            )
            raise InvalidDisplacementsSyntaxError(msg)
        # check if the last part is a direction
        try:
            self.direction = Direction(dir_str)
        except InvalidDisplacementsSyntaxError as err:
            msg = ('Unable to parse direction information from line in '
                   f'GEO_DELTA block: {self.raw_line}')
            raise InvalidDisplacementsSyntaxError(msg) from err

        # parse the rest of the left hand side into targets
        self.targets = Targets(targets_str)

        # Right hand side
        _check_moire_tag(self._rhs)
        # parse to a range
        try:
            self.range = DisplacementsRange(self._rhs)
        except ValueError as err:
            msg = ('Unable to parse range information from line in '
                   f'GEO_DELTA block: {self.raw_line}')
            raise InvalidDisplacementsSyntaxError(msg) from err


    def __repr__(self):
        """Return the string representation of the line."""
        return f'{self.targets} {self.direction} = {self.range}'


def _get_target(label, which):
    if which is None:
        return BSTarget(label)
    return BSTarget(f'{label} {which}')


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


def separate_direction_from_targets(targets_and_direction: str):
    """Separate a string into targets and direction.

    Uses a regex pattern to identify and extract an optional direction token at
    the end of the string. The rest of the string is considered the targets.

    Parameters
    ----------
    targets_and_direction : str
        The string to be separated.

    Returns
    -------
    tuple
        A tuple containing the targets string and the direction string. Either
        string may be empty if the corresponding part was not found.
    """
    match = _DIR_AT_END.search(targets_and_direction)
    if not match:
        # no direction found
        return targets_and_direction.strip(), ''

    dir_str = match.group('dir')
    # everything before the direction is the targets
    return targets_and_direction[: match.start()].strip(), dir_str.strip()


def _check_moire_tag(line_rhs):
    """Check if the right hand side of the line contains a moire tag."""
    if 'moire' in line_rhs.lower():
        msg = ('Moiré structures are not yet supported.')
        raise NotImplementedError(msg)
