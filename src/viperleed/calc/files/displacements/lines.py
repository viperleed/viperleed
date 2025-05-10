"""Module lines."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__created__ = '2024-10-04'

import re
from abc import ABC, abstractmethod
from collections import namedtuple

from viperleed_jax.files.displacements.regex import DIRECTION_PATTERN
from viperleed_jax.perturbation_type import PerturbationType, PerturbationTypeError

from .direction import Direction
from .errors import InvalidDisplacementsSyntaxError
from .range import DisplacementsRange
from .targeting import Targets, TargetingError

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

    @abstractmethod
    def block_name(self):
        """Name of the Block in the DISPLACEMENTS file."""

    def _parse_targets(self, targets_str):
        try:
            return  Targets(targets_str)
        except TargetingError as err:
            msg = (
                'Unable to parse target information from line in '
                f'{self.block_name} block: {self.raw_line}.'
            )
            raise InvalidDisplacementsSyntaxError(msg) from err

    def _parse_direction(self, dir_str):
        try:
            return Direction(dir_str)
        except InvalidDisplacementsSyntaxError as err:
            msg = (
                'Unable to parse direction information from line in '
                f'{self.block_name} block: {self.raw_line}.'
            )
            raise InvalidDisplacementsSyntaxError(msg) from err

    def _parse_range(self, range_str):
        try:
            return DisplacementsRange(range_str)
        except ValueError as err:
            msg = (
                'Unable to parse range information from line in '
                f'{self.block_name} block: {self.raw_line}.'
            )
            raise InvalidDisplacementsSyntaxError(msg) from err

class GeoDeltaLine(ParsedLine):
    """Class to parse lines in the GEO_DELTA block of DISPLACEMENTS.

    Lines in the GEO_DELTA block are of the form:
        <target> [, <target>] <direction> = <range>
    where <target>, <direction>, and <range> are tokes that are parsed by the
    `Targets`, `Direction`, and `DisplacementsRange` classes, respectively.
    """

    block_name = 'GEO_DELTA'

    def __init__(self, line: str):
        super().__init__(line)

        # Left hand side
        # check if the last part is a direction
        targets_str, dir_str = separate_direction_from_targets(self._lhs)
        if not targets_str or not dir_str:
            msg = (
                f'Invalid GEO_DELTA line format: "{self._lhs}". '
                'Expected format: "<targets> [, <target>] <direction> = '
                '<range>".'
            )
            raise InvalidDisplacementsSyntaxError(msg)

        # parse the into targets and direction
        self.targets = self._parse_targets(targets_str)
        self.direction = self._parse_direction(dir_str)

        # Right hand side
        _check_moire_tag(self._rhs)
        # parse to a range
        self.range = self._parse_range(self._rhs)

    def __repr__(self):
        """Return the string representation of the line."""
        return f'{self.targets} {self.direction} = {self.range}'


class VibDeltaLine(ParsedLine):
    """Class to parse lines in the VIB_DELTA block of DISPLACEMENTS.

    Lines in the VIB_DELTA block are of the form:
        <target> [, <target>] = <range>
    where <target>, and <range> are tokes that are parsed by the
    `Targets`, `Direction`, and `DisplacementsRange` classes, respectively.
    """

    block_name = 'VIB_DELTA'

    def __init__(self, line: str):
        super().__init__(line)

        # Left hand side
        # check if the last part is a direction
        targets_str, dir_str = separate_direction_from_targets(self._lhs)
        if dir_str:
            msg = (
                f'Invalid VIB_DELTA line format: "{self._raw_line}". '
                'Expected format: "<targets> [, <target>] = <range>".'
            )
            raise InvalidDisplacementsSyntaxError(msg)

        # parse the into targets and direction
        self.targets = self._parse_targets(targets_str)

        # Right hand side
        _check_moire_tag(self._rhs)
        # parse to a range
        self.range = self._parse_range(self._rhs)

    def __repr__(self):
        """Return the string representation of the line."""
        return f'{self.targets} = {self.range}'


def _get_target(label, which):
    if which is None:
        return BSTarget(label)
    return BSTarget(f'{label} {which}')



class OccDeltaLine:
    """Class to parse lines in the OCC_DELTA block of DISPLACEMENTS.

    Lines in the OCC_DELTA block are of the form:
        <target> [, <target>] = <element> <range> [, <element> <range> ...]
    where <target>, <element> and <range> are tokes that are parsed by the
    `Targets`, `Element` and `DisplacementsRange` classes, respectively.
    """ # TODO: Element class?

    block_type = 'OCC_DELTA'

    def __init__(self, line: str):
        super().__init__(line)

        # Left hand side
        # check if the last part is a direction
        targets_str, dir_str = separate_direction_from_targets(self._lhs)
        if dir_str:
            msg = (
                f'Invalid OCC_DELTA line format: "{self._raw_line}". '
                'Expected format: "<targets> [, <target>] = <range>".'
            )
            raise InvalidDisplacementsSyntaxError(msg)

        # parse the into targets and direction
        self.targets = self._parse_targets(targets_str)

        # TODO RHS parsing


    def __repr__(self):
        """Return the string representation of the line."""
        # TODO

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
