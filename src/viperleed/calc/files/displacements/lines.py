"""Module lines."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__created__ = '2024-10-04'

import re
from abc import ABC, abstractmethod
from collections import namedtuple

from viperleed_jax.files.displacements.regex import DIRECTION_PATTERN
from viperleed_jax.perturbation_type import PerturbationType, PerturbationTypeError

from .direction import DirectionToken
from .errors import InvalidDisplacementsSyntaxError
from .range import RangeToken
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
            return DirectionToken(dir_str)
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
    `Targets`, `Direction`, and `RangeToken` classes, respectively.
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
    `Targets`, `Direction`, and `RangeToken` classes, respectively.
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
    """Class to parse lines in the DISPLACEMENTS block of DISPLACEMENTS.

    Lines in the OCC_DELTA block are either of the form:
        <type> <target> [, <target> ...] = [<operation>] <target>
    where <target>, ... are tokes that are parsed by the # TODO

    or:
        <type> <target>, <target> [, <target> ...] = linked
    The latter syntax is a shorthand for direct linking of all targets on the
    left hand side. The 'linked' assignment will raise an
    InvalidDisplacementsSyntaxError if only one target is specified on the left
    hand side.

    block_type = 'CONSTRAIN'

    """
    def __init__(self, line: str):
        super().__init__(line)

        # check for deprecated 'offset' tag
        if 'offset' in self._rhs.lower:
            msg = ('Offset assignment in the CONSTRAIN block is deprecated. '
                   'Use the OFFSETS block instead.')
            raise InvalidDisplacementsSyntaxError(msg)



class OffsetsLine:
    """Class to parse lines in the OFFSETS block of DISPLACEMENTS.

    Lines in the OFFSETS block are of the form:
        <type> <target> [, <target> ...] = <offset>
    where <target>, ... are tokes that are parsed by the # TODO
    """

    block_type = 'OFFSET'

    def __init__(self, line):
        super().__init__(line)
        # TODO



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
