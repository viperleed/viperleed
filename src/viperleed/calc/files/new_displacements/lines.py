"""Module lines of viperleed.calc.files.new_displacements."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2024-10-04'
__license__ = 'GPLv3+'

from abc import ABC, abstractmethod
import logging
import re

import numpy as np

from viperleed.calc.classes.perturbation_mode import PerturbationMode

from .errors import DisplacementsSyntaxError
from .tokens import (
    CartesianDirectionToken,
    ElementToken,
    LinearOperationToken,
    ModeToken,
    OffsetToken,
    RangeToken,
    TargetToken,
    TokenParserError,
    TotalOccupationToken,
)

SEARCH_HEADER_PATTERN = re.compile(
    r'^==\s+(?i:search)\s+(?P<label>.*)$', re.IGNORECASE
)
SECTION_HEADER_PATTERN = re.compile(
    r'^=\s*(?P<section>GEO_DELTA|VIB_DELTA|OCC_DELTA|CONSTRAIN)$',
    re.IGNORECASE,
)
OFFSETS_HEADER_PATTERN = re.compile(r'^=+\s*(?i:OFFSET[S]?)$', re.IGNORECASE)
LOOP_START_PATTERN = re.compile(r'<loop>', re.IGNORECASE)
LOOP_END_PATTERN = re.compile(r'<\\loop>|</loop>', re.IGNORECASE)


_BELOW_DEBUG = 2

logger = logging.getLogger(__name__)

# precompiled direction-at-end regex for separating out the direction token
_FLOAT_RE = r'-?\d+(?:\.\d+)?'
_DIRECTION_BRACKETS_RE = rf'\[\s*{_FLOAT_RE}(?:\s+{_FLOAT_RE})*\s*\]'
DIRECTION_PATTERN = (
    r'(?P<direction>('
    # longer, more specific first:
    rf'azi\((?:ab|xy)?{_DIRECTION_BRACKETS_RE}\)'
    rf'|r\((?:ab|xy)?{_DIRECTION_BRACKETS_RE}\)'
    # 3-letter shorthands *with* optional bracket
    rf'|xyz(?:{_DIRECTION_BRACKETS_RE})?'
    rf'|abc(?:{_DIRECTION_BRACKETS_RE})?'
    # 2-letter shorthands *with* optional bracket
    rf'|xy(?:{_DIRECTION_BRACKETS_RE})?'
    rf'|ab(?:{_DIRECTION_BRACKETS_RE})?'
    # bare bracketed vector
    rf'|{_DIRECTION_BRACKETS_RE}'
    # single-letter shorthands
    r'|[xyzabc]'
    r'))'
)
_DIR_AT_END = re.compile(rf'(?P<dir>{DIRECTION_PATTERN})\s*$')


class HeaderLine(ABC):
    """Base class for header lines in the DISPLACEMENTS file.

    This class is used to parse header and loop marker lines that may contain
    a label or section name.
    """

    @abstractmethod
    def __init__(self, line):
        """Initialize the header line with a line string."""

    @abstractmethod
    def __str__(self):
        """Return the string representation of the header line."""


class OffsetsHeaderLine(HeaderLine):
    """Class to parse the offsets header line in the DISPLACEMENTS file.

    The offsets header line is of the form:
        = OFFSETS
    """

    def __init__(self, line):
        """Initialize the OffsetsHeaderLine with a line string."""
        if not OFFSETS_HEADER_PATTERN.match(line.strip()):
            msg = f'Invalid offsets header line: {line!r}.'
            raise DisplacementsSyntaxError(msg)
        self.section = 'OFFSETS'

    def __str__(self):
        """Return the string representation of the offsets header."""
        return '== OFFSETS'


class SearchHeaderLine(HeaderLine):
    """Class to parse the search header line in the DISPLACEMENTS file.

    The search header line is of the form:
        == SEARCH <label>
    where <label> is a string that identifies the search block.
    """

    def __init__(self, line):
        """Initialize the SearchHeaderLine with a line string."""
        match = SEARCH_HEADER_PATTERN.match(line.strip())
        if not match:
            msg = f'Invalid search header line: "{line}".'
            raise DisplacementsSyntaxError(msg)
        self.label = match['label'].strip()

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
            msg = f'Invalid section header line: "{section}".'
            raise DisplacementsSyntaxError(msg)
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
            msg = f'Invalid loop marker line: "{line}".'
            raise DisplacementsSyntaxError(msg)

    def __str__(self):
        """Return the string representation of the loop marker."""
        if self.kind == 'start':
            return '<loop>'
        return '</loop>'


class ParsedLine(ABC):
    """Base class for parsing of non-header lines in DISPLACEMENTS.

    The class is used to parse lines in the displacements file. It can
    be initialized with a line string and will parse it into its
    components.
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

    def _parse_mode(self, mode_str):
        try:
            return ModeToken(mode_str)
        except TokenParserError as err:
            msg = (
                'Unable to parse <type> information from line in '
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

    def __str__(self):
        """Return the string representation of the line."""
        lhs = self._format_lhs_str()
        rhs = self._format_rhs_str()
        return f'{type(self).__name__}({lhs} = {rhs})'

    @abstractmethod
    def _format_lhs_str(self):
        """Return a string representation of the left hand side."""

    @abstractmethod
    def _format_rhs_str(self):
        """Return a string representation of the right hand side."""


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
            raise DisplacementsSyntaxError(self.invalid_format_msg)

        # parse into targets and direction
        self.targets = self._parse_targets(targets_str)
        self.direction = self._parse_direction(dir_str)

        # Right hand side
        _check_moire_tag(self._rhs)
        # parse to a range
        self.range = self._parse_range(self._rhs)

    def _format_lhs_str(self):
        lhs = ', '.join(str(t) for t in self.targets)
        lhs += f' {self.direction}'
        return lhs

    def _format_rhs_str(self):
        return f'{self.range}'


class VibDeltaLine(ParsedLine):
    """Class to parse lines in the VIB_DELTA block of DISPLACEMENTS.

    Lines in the VIB_DELTA block are of the form:
        <target> [, <target>] = <range>
    where <target>, and <range> are tokens that are parsed by the
    `Targets` and `RangeToken` classes, respectively.
    """

    block_name = 'VIB_DELTA'
    expected_format = '<target> [, <target>] = <range>'

    def __init__(self, line: str):
        super().__init__(line)

        # Left hand side
        # check if the last part is a direction
        targets_str, dir_str = separate_direction_from_targets(self._lhs)
        if dir_str:
            raise DisplacementsSyntaxError(self.invalid_format_msg)

        # parse the targets
        self.targets = self._parse_targets(targets_str)

        # Right hand side
        _check_moire_tag(self._rhs)
        # parse to a range
        self.range = self._parse_range(self._rhs)

    def _format_lhs_str(self):
        return ', '.join(str(t) for t in self.targets)

    def _format_rhs_str(self):
        return f'{self.range}'


class OccDeltaLine(ParsedLine):
    """Class to parse lines in the OCC_DELTA block of DISPLACEMENTS.

    Lines in the OCC_DELTA block are of the form:
        <target> [, <target>] = <element> <range> [, <element> <range> ...]
    where <target>, <element> and <range> are tokens that are parsed by the
    `Targets`, `Element` and `RangeToken` classes, respectively.
    """

    block_name = 'OCC_DELTA'
    expected_format = (
        '<target> [, <target>] = <element> <range> [, <element> <range>]'
    )

    def __init__(self, line: str):
        super().__init__(line)

        # Left hand side
        # check if the last part is a direction
        targets_str, dir_str = separate_direction_from_targets(self._lhs)
        if dir_str:
            raise DisplacementsSyntaxError(self.invalid_format_msg)

        # check for deprecated combined element ranges
        _check_combined_element_ranges(self._rhs)

        # parse the targets
        self.targets = self._parse_targets(targets_str)

        # Right hand side
        # consists of one or more <element> <range> combinations
        element_ranges = tuple(
            self._parse_element_and_range(er) for er in self._rhs.split(',')
        )
        if len(element_ranges) < 1:
            # must contain at least one pair
            raise DisplacementsSyntaxError(self.invalid_format_msg)
        self.element_ranges = tuple(element_ranges)

    def _format_lhs_str(self):
        return ', '.join(str(t) for t in self.targets)

    def _format_rhs_str(self):
        return ', '.join(
            f'{elem_tok} {range_tok}'
            for elem_tok, range_tok in self.element_ranges
        )

    def _parse_element_and_range(self, el_range_str):
        """Parse a string containing an element and range token."""
        try:
            # may raise ValueError if Element or range is missing
            elem_str, range_str = el_range_str.strip().split(' ', maxsplit=1)
        except ValueError as err:
            msg = (
                'Unable to parse chemical ranges from line in '
                f'{self.block_name} block: {self.raw_line}.'
            )
            raise DisplacementsSyntaxError(msg) from err
        return (
            self._parse_element(elem_str),
            self._parse_range(range_str),
        )


class ConstraintLine(ParsedLine):
    """Class to parse lines in the CONSTRAIN block of DISPLACEMENTS.

    Lines in the CONSTRAIN block are given in the form of :
        <mode> <target> [, <target> ...] = [<linear_operation>] <target>
    Alternatively, for geometric, vibrational or occupational parameters the
    syntax
        <mode> <target_1> [, <target_2> ...] = linked
    is an allowed shorthand notation to easily create direct link.
    It will be treated as equivalent to
        <mode> <target_1> [, <target_2> ...] = <target_1>
    where the linear operation implicitly is the identity matrix.

    Additionally, for occupational parameters, the syntax
        <mode> <target_1> [, <target_2> ...] = total <total_occupation>
    is allowed, which will set a constant sum of the occupations, i.e. it holds
    the amount of vacancies constant, but allows the individual occupations to
    vary freely. This can, for example, be used to vary the occupation in an
    alloy, without introducing vacancies.
    """

    block_name = 'CONSTRAIN'
    expected_format = (
        '<type> <target> [, <target>] = [<linear_operation>] <target>'
    )

    def __init__(self, line: str):
        super().__init__(line)
        self.is_total_occupation = False

        # Left hand side
        self._treat_lhs()
        # Right hand side
        self._treat_rhs()

    def _treat_lhs(self):
        lhs_parts = self._lhs.split()
        if len(lhs_parts) < 2:  # at least mode and one target
            raise DisplacementsSyntaxError(self.invalid_format_msg)

        self.type = self._parse_mode(lhs_parts[0])

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
                'Direction tokens are not allowed in CONSTRAIN block.'
            )
            # TODO, may be implemented in a future version: Do we want this?
            # It would be generally possible (e.g. link only z, but not xy for
            # some atoms, reducing DOF from 6 to 5), but I'm not sure if one
            # would realistically use this.

        # parse targets
        self.targets = self._parse_targets(targets_str)

    def _treat_rhs(self):
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

        # check for 'total' tag
        if self._rhs.lower().strip().startswith('total '):
            logger.log(_BELOW_DEBUG, 'Detected "total" tag.')

            if self.mode.mode is not PerturbationMode.OCC:
                raise DisplacementsSyntaxError(
                    'The "total" tag is only allowed for occupational '
                    'constraints.'
                )

            # set flag for special treatment
            self.is_total_occupation = True

            # parse the total occupation value
            total_occupation_str = self._rhs[len('total ') :].strip()

            # this is a special case, where the linear operation is
            self.linear_operation = TotalOccupationToken(total_occupation_str)
            self.link_target = None
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
                'The target part of the CONSTRAIN line must contain exactly '
                'one target.'
            )

        # If we reach here, we have a valid target part
        self.link_target = link_targets[0]

        # if the operation part is empty, default to identity
        self.linear_operation = (
            self._parse_linear_operation(op_part)
            if op_part
            else LinearOperationToken.from_array(np.eye(1))
        )

    def _format_lhs_str(self):
        return ', '.join(str(t) for t in self.targets)

    def _format_rhs_str(self):
        return f'{self.linear_operation} {self.link_target}'


class OffsetsLine(ParsedLine):
    """Class to parse lines in the OFFSETS block of DISPLACEMENTS.

    Lines in the OFFSETS block are of the form:
        <mode> <target> [, <target> ...] [<direction>] = <offset>
    where <target>, ... are tokens that are parsed by the OffsetToken class.
    A direction token must be specified if and only if the mode of the offset
    is geometric.
    """

    block_name = 'OFFSET'
    expected_format = (
        '<type> <target> [, <target> ...] [<direction>] = <offset>'
    )

    def __init__(self, line):
        super().__init__(line)
        self.direction = None

        # expected DOF is 1 for vib and occ
        expected_dof = 1

        # parse LHS
        parts = self._lhs.split()
        if len(parts) < 2:  # at least mode and one target
            raise DisplacementsSyntaxError(self.invalid_format_msg)

        self.mode = self._parse_mode(parts[0])
        targets_str, dir_str = separate_direction_from_targets(
            ' '.join(parts[1:])
        )

        # parse targets
        self.targets = self._parse_targets(targets_str)

        if self.mode.mode is PerturbationMode.GEO:
            # expect and parse direction specifier
            # will raise if no direction is given
            self.direction = self._parse_direction(dir_str)
            # if geometric, expected DOF is given by direction
            expected_dof = self.direction.dof

        if self.mode.mode is not PerturbationMode.GEO and dir_str:
            raise DisplacementsSyntaxError(
                'Direction tokens in the OFFSETS block are only allowed for '
                'geometric offsets.'
            )

        # parse RHS into offset
        try:
            offset = OffsetToken(self._rhs)
        except TokenParserError as err:
            raise DisplacementsSyntaxError(self.invalid_format_msg) from err

        # check that the offset DOF matches the expected DOF
        if offset.dof != expected_dof:
            msg = (
                f'Offset DOF ({offset.dof}) does not match expected DOF '
                f'({expected_dof}) for line in {self.block_name} block: '
                f'"{self.raw_line}".'
            )
            raise DisplacementsSyntaxError(msg)
        self.offset = offset

    def _format_lhs_str(self):
        lhs = ', '.join(str(t) for t in self.targets)
        if self.mode.mode is PerturbationMode.GEO:
            lhs += f' {self.direction}'
        return f'{self.mode} {lhs}'

    def _format_rhs_str(self):
        return f'{self.offset}'


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
    matches = tuple(_DIR_AT_END.finditer(targets_and_direction))
    if len(matches) > 1:
        raise ValueError('Only one directional specification is allowed.')
    if matches:
        match = matches[0]
    else:
        # no direction found
        return targets_and_direction.strip(), ''

    # everything before the direction is the targets
    return targets_and_direction[: match.start()].strip(), match['dir'].strip()


def _check_moire_tag(line_rhs):
    """Check if the right hand side of the line contains a moire tag."""
    if 'moire' in line_rhs.lower():
        msg = 'Moiré structures are not yet supported.'
        raise NotImplementedError(msg)


def _check_combined_element_ranges(line_rhs):
    """Check if the rhs of the line contains combined element ranges."""
    if '+' in line_rhs:
        msg = (
            'Combined element ranges have been deprecated. Specify each '
            'element range individually.'
        )
        raise DisplacementsSyntaxError(msg)
