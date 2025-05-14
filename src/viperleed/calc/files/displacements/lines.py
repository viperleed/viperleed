"""Module lines."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__created__ = '2024-10-04'

import logging
import re
from abc import ABC, abstractmethod
from collections import namedtuple

import numpy as np

from viperleed_jax.files.displacements.perturbation_type import (
    PerturbationType,
)
from viperleed_jax.files.displacements.regex import DIRECTION_PATTERN

from .errors import InvalidDisplacementsSyntaxError
from .tokens import (
    DirectionToken,
    ElementToken,
    LinearOperationToken,
    OffsetToken,
    RangeToken,
    TargetToken,
    TokenParserError,
    TypeToken,
)

LoopMarkerLine = namedtuple('LoopMarkerLine', ['type'])
SearchHeaderLine = namedtuple('SearchHeaderLine', ['label'])
SectionHeaderLine = namedtuple('SectionHeaderLine', ['section'])

_BELOW_DEBUG = 2

logger = logging.getLogger(__name__)

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
        logger.log(_BELOW_DEBUG, f'Parsing DISPLACEMENTS: {line}')
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

    @abstractmethod
    def expected_format(self):
        """Name of the Block in the DISPLACEMENTS file."""

    def _parse_targets(self, targets_str):
        target_parts = targets_str.split(',')
        try:
            return tuple(TargetToken(part) for part in target_parts)
        except TokenParserError as err:
            msg = (
                'Unable to parse <target> tokens from line in '
                f'{self.block_name} block: {self.raw_line}.'
            )
            raise InvalidDisplacementsSyntaxError(msg) from err

    def _parse_direction(self, dir_str):
        try:
            return DirectionToken(dir_str)
        except TokenParserError as err:
            msg = (
                'Unable to parse <direction> token from line in '
                f'{self.block_name} block: {self.raw_line}.'
            )
            raise InvalidDisplacementsSyntaxError(msg) from err

    def _parse_range(self, range_str):
        try:
            return RangeToken(range_str)
        except TokenParserError as err:
            msg = (
                'Unable to parse <range> token from line in '
                f'{self.block_name} block: {self.raw_line}.'
            )
            raise InvalidDisplacementsSyntaxError(msg) from err

    def _parse_type(self, type_str):
        try:
            return TypeToken(type_str)
        except TokenParserError as err:
            msg = (
                'Unable to parse <type> information from line in '
                f'{self.block_name} block: {self.raw_line}.'
            )
            raise InvalidDisplacementsSyntaxError(msg) from err

    def _parse_element(self, element_str):
        try:
            return ElementToken(element_str)
        except TokenParserError as err:
            msg = (
                'Unable to parse <element> information from line in '
                f'{self.block_name} block: {self.raw_line}.'
            )
            raise InvalidDisplacementsSyntaxError(msg) from err

    def _parse_linear_operation(self, operation_str):
        try:
            self.linear_operation = LinearOperationToken(operation_str)
        except TokenParserError as err:
            msg = (
                'Unable to parse <linear_operation> token from line in '
                f'{self.block_name} block: {self.raw_line}.'
            )
            raise InvalidDisplacementsSyntaxError(msg) from err

    @property
    def invalid_format_msg(self):
        """Return a string with a general invalid format error message."""
        return (
            f'Invalid {self.block_name} line format: "{self.raw_line}". '
            f'Expected format: "{self.expected_format}".'
        )


class GeoDeltaLine(ParsedLine):
    """Class to parse lines in the GEO_DELTA block of DISPLACEMENTS.

    Lines in the GEO_DELTA block are of the form:
        <target> [, <target>] <direction> = <range>
    where <target>, <direction>, and <range> are tokes that are parsed by the
    `Targets`, `Direction`, and `RangeToken` classes, respectively.
    """

    block_name = 'GEO_DELTA'
    expected_format = '<target> [, <target>] <direction> = <range>'

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
    expected_format = '<targets> [, <target>] = <range>'

    def __init__(self, line: str):
        super().__init__(line)

        # Left hand side
        # check if the last part is a direction
        targets_str, dir_str = separate_direction_from_targets(self._lhs)
        if dir_str:
            raise InvalidDisplacementsSyntaxError(self.invalid_format_msg)

        # parse the into targets and direction
        self.targets = self._parse_targets(targets_str)

        # Right hand side
        _check_moire_tag(self._rhs)
        # parse to a range
        self.range = self._parse_range(self._rhs)

    def __repr__(self):
        """Return the string representation of the line."""
        return f'{self.targets} = {self.range}'



class OccDeltaLine(ParsedLine):
    """Class to parse lines in the OCC_DELTA block of DISPLACEMENTS.

    Lines in the OCC_DELTA block are of the form:
        <target> [, <target>] = <element> <range> [, <element> <range> ...]
    where <target>, <element> and <range> are tokes that are parsed by the
    `Targets`, `Element` and `RangeToken` classes, respectively.
    """

    block_name = 'OCC_DELTA'
    expected_format = (
        '<target> [, <target>] = <element> <range> [, <element> <range> ...]'
    )

    def __init__(self, line: str):
        super().__init__(line)

        # Left hand side
        # check if the last part is a direction
        targets_str, dir_str = separate_direction_from_targets(self._lhs)
        if dir_str:
            raise InvalidDisplacementsSyntaxError(self.invalid_format_msg)

        # parse the into targets and direction
        self.targets = self._parse_targets(targets_str)

        # Right hand side
        # consists of one or more <element> <range> combinations
        elem_ranges_strs = self._rhs.split(',')
        # iterate over all
        element_ranges = []
        for er in elem_ranges_strs:
            try:
                # may raise ValueError if Element or range is missing
                elem_str, range_str = er.strip().split(' ', maxsplit=1)
            except ValueError as err:
                msg = (
                    'Unable to parse chemical ranges from line in '
                    f'{self.block_name} block: {self.raw_line}.'
                )
                raise InvalidDisplacementsSyntaxError(msg) from err
            element_ranges.append(
                (
                    self._parse_element(elem_str),
                    self._parse_range(range_str),
                )
            )
        if len(element_ranges) < 1:
            # must contain at least one pair
            raise InvalidDisplacementsSyntaxError(self.invalid_format_msg)
        self.element_ranges = tuple(element_ranges)

    def __repr__(self):
        """Return the string representation of the line."""
        txt = f'{self.targets[0]}'
        for target in self.targets[1:]:
            txt += f', {target}'
        txt += f' = {self.element_ranges[0]}'
        for er in self.element_ranges[1:]:
            txt += f', {er}'
        return txt


class ConstraintLine(ParsedLine):
    """Class to parse lines in the CONSTRAIN block of DISPLACEMENTS.

    Lines in the OCC_DELTA block are given in the form of :
        <type> <target> [, <target> ...] = [<linear_operation>] <target>
    Alternatively, for geometric, vibrational or occupational parameters the
    syntax
        <type> <target_1> [, <target_2> ...], <target_n> = linked
    an allowed shorthand notation for direct links and requires at least two
    targets to be specified on the left hand side. It will be treated as
    equivalently to
        <target_1> [, <target_2> ...], <target_{n-1}> = <target_n>
    .
    """

    block_name = 'CONSTRAIN'
    expected_format = '<type> <target> [, <target>] = [<linear_operation>] <target>'

    def __init__(self, line: str):
        super().__init__(line)

        # Left hand side
        lhs_parts = self._lhs.split()
        if len(lhs_parts) < 2:  # at least type and one target
            raise InvalidDisplacementsSyntaxError(self.invalid_format_msg)

        self.type = self._parse_type(lhs_parts[0])

        # TODO: Implement Domains here
        if self.type.type is PerturbationType.DOMAIN:
            raise NotImplementedError(
                'Domain constraints are not yet supported.'
            )

        targets_str, dir_str = separate_direction_from_targets(
            ''.join(lhs_parts[1:])
        )

        # parse targets
        self.targets = self._parse_targets(targets_str)

        # Right hand side

        # check for deprecated 'offset' tag
        if 'offset' in self._rhs.lower():
            msg = (
                'Offset assignment in the CONSTRAIN block is deprecated. '
                'Use the OFFSETS block instead.'
            )
            raise InvalidDisplacementsSyntaxError(msg)

        # check for 'linked' tag
        if self._rhs.lower().strip() == 'linked':
            logger.log(_BELOW_DEBUG, 'Detected "linked" tag.')
            if len(self.targets) < 2:
                raise InvalidDisplacementsSyntaxError(
                    'Direct link assignment using the "linked" tag requires '
                    'specifying at least two targets.'
                )
            # "move" the last target to the rhs
            self.link_target = self.targets[-1]
            self.targets = self.targets[:-1]
            # treat as if array is identity
            self.linear_operation = LinearOperationToken.from_array(np.eye(1))
            return

        # The default case is to treat it as <linear_operation> and <target>.
        # It's not immediately obvious where to split the tokens since both
        # token types may contain spaces. Instead we iterate over white-space
        # split parts from the right (since the <linear_operation> is optional)
        # and try to parse the target.
        rhs_parts = self._rhs.split('')
        for i in range(1, len(rhs_parts) + 1):
            # Try treating the last i parts as the target
            op_part = ' '.join(rhs_parts[:-i]) if i < len(rhs_parts) else ''
            target_part = ' '.join(rhs_parts[-i:])
            try:
                link_target = self._parse_targets(target_part)
                break
            except InvalidDisplacementsSyntaxError:
                continue
        else:
            raise InvalidDisplacementsSyntaxError(self.invalid_format_msg)

        # If we reach here, we have a valid target part
        self.link_target = link_target

        # if the operation part is empty, default to identity
        if not op_part:
            self.linear_operation = LinearOperationToken.from_array(np.eye(1))
            return

        self.linear_operation = self._parse_linear_operation(op_part)


class OffsetsLine(ParsedLine):
    """Class to parse lines in the OFFSETS block of DISPLACEMENTS.

    Lines in the OFFSETS block are of the form:
        <type> <target> [, <target> ...] [<direction>] = <offset>
    where <target>, ... are tokes that are parsed by the OffsetToken class.
    A direction token must be specified if and only if the type of the offset
    is geometric.
    """

    block_name = 'OFFSET'
    expected_format = ('<type> <target> [, <target> ...] [<direction>] '
                       '= <offset>')

    def __init__(self, line):
        super().__init__(line)
        self.direction = None

        # parse LHS

        parts = self._lhs.split()
        if len(parts) < 2:  # at least type and one target
            raise InvalidDisplacementsSyntaxError(self.invalid_format_msg)

        self.type = self._parse_type(parts[0])
        targets_str, dir_str = separate_direction_from_targets(
            ''.join(parts[1:]))

        # parse targets
        self.targets = self._parse_targets(targets_str)

        if self.type.type is PerturbationType.GEO:
            # expect and parse direction specifier
            # will raise if no direction is given
            self.direction = self._parse_direction(dir_str)

        if self.type.type is not PerturbationType.GEO and dir_str:
            raise InvalidDisplacementsSyntaxError(
                'Direction tokens in the OFFSETS block are only allowed for '
                'geometric offsets.'
            )

        # parse RHS into offset
        try:
            self.offset = OffsetToken(self._rhs)
        except TokenParserError as err:
            raise InvalidDisplacementsSyntaxError(
                self.invalid_format_msg) from err

    def __repr__(self):
        """Return the string representation of the line."""
        txt = f'{self.type} {self.targets[0]}'
        for target in self.targets[1:]:
            txt += f', {target}'
        if self.type.type is PerturbationType.GEO:
            txt += f', {self.direction}'
        txt += f' = {self.offset}'
        return txt


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
