"""Module lines of viperleed.calc.files.new_displacements."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2024-10-04'
__license__ = 'GPLv3+'

import logging
import re
from abc import ABC, abstractmethod

import numpy as np

from viperleed.calc.classes.perturbation_mode import PerturbationMode

from .errors import DisplacementsSyntaxError
from .tokens import (
    CartesianDirectionToken,
    ElementToken,
    LinearOperationToken,
    OffsetToken,
    RangeToken,
    TargetToken,
    TokenParserError,
    ModeToken,
)


SEARCH_HEADER_PATTERN = re.compile(r"^==\s+(?i:search)\s+(?P<label>.*)$")
SECTION_HEADER_PATTERN = re.compile(
    r"^=\s*(?P<section>GEO_DELTA|VIB_DELTA|OCC_DELTA|CONSTRAIN)$"
)
OFFSETS_HEADER_PATTERN = re.compile(r"^=+\s*(?i:OFFSET[S]?)$")
LOOP_START_PATTERN = re.compile(r'<loop>')
LOOP_END_PATTERN = re.compile(r'<\\loop>|</loop>')


_BELOW_DEBUG = 2

logger = logging.getLogger(__name__)

# precompiled direction-at-end regex for separating out the direction token
DIRECTION_PATTERN = (
    r'(?P<direction>('
    # longer, more specific first:
    r'azi\((?:ab|xy)?\[\s*-?\d+(?:\.\d+)?(?:\s+-?\d+(?:\.\d+)?)*\s*\]\)'
    r'|r\((?:ab|xy)?\[\s*-?\d+(?:\.\d+)?(?:\s+-?\d+(?:\.\d+)?)*\s*\]\)'
    # 3-letter shorthands *with* optional bracket
    r'|xyz(?:\[\s*-?\d+(?:\.\d+)?(?:\s+-?\d+(?:\.\d+)?)*\s*\])?'
    r'|abc(?:\[\s*-?\d+(?:\.\d+)?(?:\s+-?\d+(?:\.\d+)?)*\s*\])?'
    # 2-letter shorthands *with* optional bracket
    r'|xy(?:\[\s*-?\d+(?:\.\d+)?(?:\s+-?\d+(?:\.\d+)?)*\s*\])?'
    r'|ab(?:\[\s*-?\d+(?:\.\d+)?(?:\s+-?\d+(?:\.\d+)?)*\s*\])?'
    # bare bracketed vector
    r'|\[\s*-?\d+(?:\.\d+)?(?:\s+-?\d+(?:\.\d+)?)*\s*\]'
    # single-letter shorthands
    r'|[xyzabc]'
    r'))'
)
_DIR_AT_END = re.compile(
    rf'(?P<dir>{DIRECTION_PATTERN})\s*$'
)


class HeaderLine(ABC):
    """Base class for header lines in the DISPLACEMENTS file.

    This class is used to parse header and loop marker lines that may contain
    a label or section name.
    """

    @abstractmethod
    def __init__(self, line: str):
        """Initialize the header line with a line string."""

    @abstractmethod
    def __str__(self):
        """Return the string representation of the header line."""
        pass

class OffsetsHeaderLine(HeaderLine):
    """Class to parse the offsets header line in the DISPLACEMENTS file.

    The offsets header line is of the form:
        = OFFSETS
    """

    def __init__(self, line: str):
        """Initialize the OffsetsHeaderLine with a line string."""
        if not OFFSETS_HEADER_PATTERN.match(line.strip()):
            raise DisplacementsSyntaxError(
                f'Invalid offsets header line: "{line}".')
        self.line = line.strip()
        self.section = 'OFFSETS'

    def __str__(self):
        """Return the string representation of the offsets header."""
        return self.line

class SearchHeaderLine(HeaderLine):
    """Class to parse the search header line in the DISPLACEMENTS file.

    The search header line is of the form:
        == SEARCH <label>
    where <label> is a string that identifies the search block.
    """

    def __init__(self, line: str):
        """Initialize the SearchHeaderLine with a line string."""
        match = SEARCH_HEADER_PATTERN.match(line.strip())
        if not match:
            raise DisplacementsSyntaxError(
                f'Invalid search header line: "{line}".')
        self.label = match.group('label').strip()

    def __str__(self):
        """Return the string representation of the search header."""
        return f'== SEARCH {self.label}'


class SectionHeaderLine(HeaderLine):
    """Class to parse section header lines in the DISPLACEMENTS file.

    Section header lines are of the form:
        = <section>
    where <section> is one of the defined sections in the DISPLACEMENTS file.
    """

    def __init__(self, section: str):
        """Initialize the SectionHeaderLine with a section name."""
        match = SECTION_HEADER_PATTERN.match(section.strip().upper())
        if not match:
            raise DisplacementsSyntaxError(
                f'Invalid section header line: "{section}".')
        # extract the section name from the match
        self.section = match.group('section').strip()

    def __str__(self):
        """Return the string representation of the section header."""
        return f'= {self.section}'


class LoopMarkerLine(HeaderLine):
    """Class to parse loop marker lines in the DISPLACEMENTS file.

    Loop marker lines are of the form:
        <loop> or </loop>
    They indicate the start or end of a loop in the DISPLACEMENTS file.
    """

    def __init__(self, line: str):
        """Initialize the LoopMarkerLine with a line string."""
        stripped_line = line.strip().lower()
        if LOOP_START_PATTERN.match(stripped_line):
            self.kind = 'start'
        elif LOOP_END_PATTERN.match(stripped_line):
            self.kind = 'end'
        else:
            raise DisplacementsSyntaxError(
                f'Invalid loop marker line: "{line}".'
            )

    def __str__(self):
        """Return the string representation of the loop marker."""
        if self.kind == 'start':
            return '<loop>'
        else:
            return '</loop>'


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
            raise DisplacementsSyntaxError(msg)

        # split the line into raw left and right hand sides
        self._lhs, self._rhs = self.raw_line.split('=')


    @abstractmethod
    def __str__(self):
        """Return the string representation of the line."""

    @property
    @abstractmethod
    def block_name(self):
        """Name of the block in the DISPLACEMENTS file."""

    @property
    @abstractmethod
    def expected_format(self):
        """Name of the block in the DISPLACEMENTS file."""

    def _parse_targets(self, targets_str):
        target_parts = targets_str.split(',')
        try:
            return tuple(TargetToken(part) for part in target_parts)
        except TokenParserError as err:
            msg = (
                'Unable to parse <target> tokens from line in '
                f'{self.block_name} block: {self.raw_line}.'
            )
            raise DisplacementsSyntaxError(msg) from err

    def _parse_direction(self, dir_str):
        try:
            return CartesianDirectionToken(dir_str)
        except TokenParserError as err:
            msg = (
                'Unable to parse <direction> token from line in '
                f'{self.block_name} block: {self.raw_line}.'
            )
            raise DisplacementsSyntaxError(msg) from err

    def _parse_range(self, range_str):
        try:
            return RangeToken(range_str)
        except TokenParserError as err:
            msg = (
                'Unable to parse <range> token from line in '
                f'{self.block_name} block: {self.raw_line}.'
            )
            raise DisplacementsSyntaxError(msg) from err

    def _parse_type(self, type_str):
        try:
            return ModeToken(type_str)
        except TokenParserError as err:
            msg = (
                'Unable to parse <type> information from line in '
                f'{self.block_name} block: {self.raw_line}.'
            )
            raise DisplacementsSyntaxError(msg) from err

    def _parse_element(self, element_str):
        try:
            return ElementToken(element_str)
        except TokenParserError as err:
            msg = (
                'Unable to parse <element> information from line in '
                f'{self.block_name} block: {self.raw_line}.'
            )
            raise DisplacementsSyntaxError(msg) from err

    def _parse_linear_operation(self, operation_str):
        try:
            return LinearOperationToken(operation_str)
        except TokenParserError as err:
            msg = (
                'Unable to parse <linear_operation> token from line in '
                f'{self.block_name} block: {self.raw_line}.'
            )
            raise DisplacementsSyntaxError(msg) from err

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
    where <target>, <direction>, and <range> are tokens that are parsed by the
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
            raise DisplacementsSyntaxError(msg)

        # parse the into targets and direction
        self.targets = self._parse_targets(targets_str)
        self.direction = self._parse_direction(dir_str)

        # Right hand side
        _check_moire_tag(self._rhs)
        # parse to a range
        self.range = self._parse_range(self._rhs)

    def __str__(self):
        """Return the string representation of the line."""
        return f'{self.targets} {self.direction} = {self.range}'


class VibDeltaLine(ParsedLine):
    """Class to parse lines in the VIB_DELTA block of DISPLACEMENTS.

    Lines in the VIB_DELTA block are of the form:
        <target> [, <target>] = <range>
    where <target>, and <range> are tokens that are parsed by the
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
            raise DisplacementsSyntaxError(self.invalid_format_msg)

        # parse the into targets and direction
        self.targets = self._parse_targets(targets_str)

        # Right hand side
        _check_moire_tag(self._rhs)
        # parse to a range
        self.range = self._parse_range(self._rhs)

    def __str__(self):
        """Return the string representation of the line."""
        return f'{self.targets} = {self.range}'



class OccDeltaLine(ParsedLine):
    """Class to parse lines in the OCC_DELTA block of DISPLACEMENTS.

    Lines in the OCC_DELTA block are of the form:
        <target> [, <target>] = <element> <range> [, <element> <range> ...]
    where <target>, <element> and <range> are tokens that are parsed by the
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
            raise AttributeError(self.invalid_format_msg)

        # check for deprecated combined element ranges
        _check_combined_element_ranges(self._rhs)

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
                raise DisplacementsSyntaxError(msg) from err
            element_ranges.append(
                (
                    self._parse_element(elem_str),
                    self._parse_range(range_str),
                )
            )
        if len(element_ranges) < 1:
            # must contain at least one pair
            raise DisplacementsSyntaxError(self.invalid_format_msg)
        self.element_ranges = tuple(element_ranges)

    def __str__(self):
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

    Lines in the CONSTRAIN block are given in the form of :
        <type> <target> [, <target> ...] = [<linear_operation>] <target>
    Alternatively, for geometric, vibrational or occupational parameters the
    syntax
        <type> <target_1> [, <target_2> ...] = linked
    is an allowed shorthand notation to easily create direct link.
    It will be treated as equivalent to
        <type> <target_1> [, <target_2> ...] = <target_1>
    where the linear operation implicitly is the identity matrix.
    """

    block_name = 'CONSTRAIN'
    expected_format = '<type> <target> [, <target>] = [<linear_operation>] <target>'

    def __init__(self, line: str):
        super().__init__(line)

        # Left hand side
        lhs_parts = self._lhs.split()
        if len(lhs_parts) < 2:  # at least type and one target
            raise DisplacementsSyntaxError(self.invalid_format_msg)

        self.type = self._parse_type(lhs_parts[0])

        # TODO: Implement Domains here
        if self.type.mode is PerturbationMode.DOM:
            raise NotImplementedError(
                'Domain constraints are not yet supported.'
            )

        targets_str, dir_str = separate_direction_from_targets(
            ' '.join(lhs_parts[1:])
        )
        if dir_str:
            raise DisplacementsSyntaxError(
                "Direction tokens are not allowed in CONSTRAIN block."
            )
            # TODO, may be implemented in a future version: Do we want this?
            # It would be generally possible (e.g. link only z, but not xy for
            # some atoms, reducing DOF from 6 to 5), but I'm not sure if one
            # would realistically use this.

        # parse targets
        self.targets = self._parse_targets(targets_str)

        # Right hand side

        # check for deprecated 'offset' tag
        if 'offset' in self._rhs.lower():
            msg = (
                'Offset assignment in the CONSTRAIN block is deprecated. '
                'Use the OFFSETS block instead.'
            )
            raise DisplacementsSyntaxError(msg)

        # check for 'linked' tag
        if self._rhs.lower().strip() == 'linked':
            logger.log(_BELOW_DEBUG, 'Detected "linked" tag.')
            # "copy" the first target to the rhs
            self.link_target = self.targets[0]
            # treat as if array is identity
            self.linear_operation = LinearOperationToken.from_array(np.eye(1))
            return

        # The default case is to treat it as [<linear_operation>] <target>
        # It's not immediately obvious where to split the tokens since both
        # token types may contain spaces.
        # Split by whitespace, parse from the right: find the last complete
        # target token and treat everything before it as the operation
        rhs_parts = self._rhs.split(' ')
        for i in range(1, len(rhs_parts) + 1):
            # Try treating the last i parts as the target
            op_part = ' '.join(rhs_parts[:-i]) if i < len(rhs_parts) else ''
            target_part = ' '.join(rhs_parts[-i:])
            try:
                link_targets = self._parse_targets(target_part)
                break
            except DisplacementsSyntaxError:
                continue
        else:
            raise DisplacementsSyntaxError(self.invalid_format_msg)

        if len(link_targets) != 1:
            raise DisplacementsSyntaxError(
                "The target part of the CONSTRAIN line must contain exactly "
                "one target."
            )

        # If we reach here, we have a valid target part
        self.link_target = link_targets[0]

        # if the operation part is empty, default to identity
        if not op_part:
            self.linear_operation = LinearOperationToken.from_array(np.eye(1))
            return

        self.linear_operation = self._parse_linear_operation(op_part)

    def __str__(self):
        """Return the string representation of the line."""
        txt = f'ConstraintLine({self.targets[0]}'
        for target in self.targets[1:]:
            txt += f', {target}'
        txt += f' = {self.linear_operation}'
        txt += f'{self.link_target})'
        return txt


class OffsetsLine(ParsedLine):
    """Class to parse lines in the OFFSETS block of DISPLACEMENTS.

    Lines in the OFFSETS block are of the form:
        <type> <target> [, <target> ...] [<direction>] = <offset>
    where <target>, ... are tokens that are parsed by the OffsetToken class.
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
            raise DisplacementsSyntaxError(self.invalid_format_msg)

        self.type = self._parse_type(parts[0])
        targets_str, dir_str = separate_direction_from_targets(
            ' '.join(parts[1:]))

        # parse targets
        self.targets = self._parse_targets(targets_str)

        if self.type.mode is PerturbationMode.GEO:
            # expect and parse direction specifier
            # will raise if no direction is given
            self.direction = self._parse_direction(dir_str)

        if self.type.mode is not PerturbationMode.GEO and dir_str:
            raise DisplacementsSyntaxError(
                "Direction tokens in the OFFSETS block are only allowed for "
                "geometric offsets."
            )

        # parse RHS into offset
        try:
            self.offset = OffsetToken(self._rhs)
        except TokenParserError as err:
            raise DisplacementsSyntaxError(self.invalid_format_msg) from err

    def __str__(self):
        """Return the string representation of the line."""
        txt = f'{self.type} {self.targets[0]}'
        for target in self.targets[1:]:
            txt += f', {target}'
        if self.type.mode is PerturbationMode.GEO:
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
    matches = list(_DIR_AT_END.finditer(targets_and_direction))
    if len(matches) > 1:
        raise ValueError("Only one directional specification is allowed.")
    elif matches:
        match = matches[0]
    else:
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

def _check_combined_element_ranges(line_rhs):
    """Check if the rhs of the line contains combined element ranges."""
    if '+' in line_rhs:
        msg = (
            'Combined element ranges have been deprecated. Use individual '
            'element ranges for each element instead, combined with '
            'constraints in the CONSTRAIN block.'
        )
        raise DisplacementsSyntaxError(msg)

