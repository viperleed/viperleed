"""Module lines of viperleed.calc.files.new_displacements."""

__authors__ = ('Alexandra Mia Imre (@alexmiame)',)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2024-10-04'
__license__ = 'GPLv3+'

from abc import ABC, abstractmethod
import logging
import re

import numpy as np

from viperleed.calc.classes.perturbation_mode import PerturbationMode
from viperleed.calc.constants import LOG_VERY_VERBOSE
from viperleed.calc.files.new_displacements.errors import (
    DisplacementsSyntaxError,
)
from viperleed.calc.files.new_displacements.tokens.base import (
    TokenParserError,
    parse_plain_float,
)
from viperleed.calc.files.new_displacements.tokens.direction import (
    make_direction_token,
)
from viperleed.calc.files.new_displacements.tokens.element import ElementToken
from viperleed.calc.files.new_displacements.tokens.linear_operation import (
    LinearOperationToken,
)
from viperleed.calc.files.new_displacements.tokens.mode import ModeToken
from viperleed.calc.files.new_displacements.tokens.offset import OffsetToken
from viperleed.calc.files.new_displacements.tokens.range_ import RangeToken
from viperleed.calc.files.new_displacements.tokens.target import TargetToken
from viperleed.calc.files.new_displacements.tokens.total_occupation import (
    TotalOccupationToken,
)

SEARCH_HEADER_PATTERN = re.compile(
    r'^==\s+(?i:search)(?:\s+(?P<label>.*))?$', re.IGNORECASE
)
MODE_HEADER_PATTERN = re.compile(
    r'^=\s*(?P<section>GEO_DELTA|VIB_DELTA|OCC_DELTA|CONSTRAIN)$',
    re.IGNORECASE,
)
OFFSETS_HEADER_PATTERN = re.compile(r'^=+\s*(?i:OFFSET[S]?)$', re.IGNORECASE)
LOOP_START_PATTERN = re.compile(r'<\s*loop\s*>', re.IGNORECASE)
LOOP_END_PATTERN = re.compile(r'<\s*[/\\]\s*loop\s*>', re.IGNORECASE)


logger = logging.getLogger(__name__)

_OFFSET_KEYWORD = 'offset'
_MIN_LHS_PARTS_WITH_OFFSET = 2  # targets and the keyword
_ZERO_OFFSET_WORDS = {'clear', 'original'}

# precompiled direction-at-end regex for separating out the direction token
_FLOAT_RE = r'-?(?:\d+(?:\.\d*)?|\.\d+)'
_ELEMENT_AND_VALUE = 2  # parts of an occupation rhs
_INTEGER_RE = re.compile(r'[+-]?\d+')
_DIRECTION_BRACKETS_RE = rf'\[\s*{_FLOAT_RE}(?:\s+{_FLOAT_RE})*\s*\]'
DIRECTION_PATTERN = (
    r'(?P<direction>('
    # longer, more specific first:
    rf'azi\((?:ab|xy)?{_DIRECTION_BRACKETS_RE}\)'
    rf'|r\((?:ab|xy)?{_DIRECTION_BRACKETS_RE}\)'
    # Cartesian shorthands: only the forms with x, y, z in order
    rf'|xyz(?:{_DIRECTION_BRACKETS_RE})?'
    rf'|xy(?:{_DIRECTION_BRACKETS_RE})?'
    r'|[xyz]'
    # Frame letters with at least one of a, b, c, and optional bracket.
    # Also catches invalid ones (e.g., 'c', 'az') so that the direction
    # token can complain.
    rf'|[xyzabc]*[abc][xyzabc]*(?:{_DIRECTION_BRACKETS_RE})?'
    # bare bracketed vector
    rf'|{_DIRECTION_BRACKETS_RE}'
    r'))'
)
_DIR_AT_END = re.compile(rf'\s+(?P<dir>{DIRECTION_PATTERN})\s*$')


class HeaderLine(ABC):
    """Base class for header lines in the DISPLACEMENTS file.

    This class is used to parse header and loop marker lines that may contain
    a label or mode name.
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

    section = 'OFFSETS'

    def __init__(self, line):
        """Initialize the OffsetsHeaderLine with a line string."""
        if not OFFSETS_HEADER_PATTERN.match(line.strip()):
            msg = f'Invalid offsets header line: {line!r}.'
            raise DisplacementsSyntaxError(msg)

    def __str__(self):
        """Return the string representation of the offsets header."""
        return f'== {self.section}'


class SearchHeaderLine(HeaderLine):
    """Class to parse the search header line in the DISPLACEMENTS file.

    The search header line is of the form:
        == SEARCH <label>
    where <label> is a string that identifies the search block. The
    label is optional: if missing, `label` is an empty string.
    """

    def __init__(self, line):
        """Initialize the SearchHeaderLine with a line string."""
        match = SEARCH_HEADER_PATTERN.match(line.strip())
        if not match:
            msg = f'Invalid search header line: "{line}".'
            raise DisplacementsSyntaxError(msg)
        self.label = (match['label'] or '').strip()

    def __str__(self):
        """Return the string representation of the search header."""
        return f'== SEARCH {self.label}'


class ModeHeaderLine(HeaderLine):
    """Class to parse mode header lines in the DISPLACEMENTS file.

    Mode header lines are of the form:
        = <mode>
    where <mode> is one of the perturbation modes that open a mode block
    in the DISPLACEMENTS file.
    """

    def __init__(self, section):
        """Initialize the ModeHeaderLine with a mode name.

        Parameters
        ----------
        section : str
            The text of the mode header line.
        """
        match = MODE_HEADER_PATTERN.match(section.strip().upper())
        if not match:
            msg = f'Invalid mode header line: "{section}".'
            raise DisplacementsSyntaxError(msg)
        # extract the mode name from the match
        self.section = match['section'].strip()

    def __str__(self):
        """Return the string representation of the mode header."""
        return f'= {self.section}'


class LoopMarkerLine(HeaderLine):
    """Class to parse loop marker lines in the DISPLACEMENTS file.

    Loop marker lines are of the form:
        <loop> or </loop>
    They indicate the start or end of a loop in the DISPLACEMENTS file.
    """

    def __init__(self, line):
        """Initialize the LoopMarkerLine with a line string.

        Parameters
        ----------
        line : str
            The text of the loop marker line.
        """
        stripped_line = line.strip().lower()
        if LOOP_START_PATTERN.fullmatch(stripped_line):
            self.kind = 'start'
        elif LOOP_END_PATTERN.fullmatch(stripped_line):
            self.kind = 'end'
        else:
            msg = f'Invalid loop marker line: "{line}".'
            raise DisplacementsSyntaxError(msg)

    @staticmethod
    def is_marker(line):
        """Return whether `line` is exactly a loop start or end tag."""
        stripped_line = line.strip().lower()
        return bool(
            LOOP_START_PATTERN.fullmatch(stripped_line)
            or LOOP_END_PATTERN.fullmatch(stripped_line)
        )

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

    def __init__(self, line):
        """Perform basic parsing and checking of the line format.

        Parameters
        ----------
        line : str
            The text of the line.
        """
        logger.log(LOG_VERY_VERBOSE, f'Parsing DISPLACEMENTS: {line}')
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

    def __init_subclass__(cls, **kwargs):
        """Check that subclasses define required class attributes."""
        super().__init_subclass__(**kwargs)
        required_attrs = (
            'block_name',  # Name of the block in the DISPLACEMENTS file
            'expected_format',  # Format of the block
        )
        for attr in required_attrs:
            if not hasattr(cls, attr):
                msg = (
                    f'Class {cls.__name__} must define '
                    f'class attribute {attr!r}'
                )
                raise TypeError(msg)

    def _reject_inline_offset(self, mode):
        """Raise if the line has an in-line 'offset' keyword.

        Parameters
        ----------
        mode : {'geo', 'vib'}
            The perturbation mode of the line, used to suggest the
            equivalent line for the OFFSETS block.

        Raises
        ------
        DisplacementsSyntaxError
            If the left-hand side ends with 'offset'.
        """
        lhs_parts = self._lhs.split()
        if len(lhs_parts) < _MIN_LHS_PARTS_WITH_OFFSET:
            return
        if lhs_parts[-1].lower() != _OFFSET_KEYWORD:
            return
        value = self._rhs.strip()
        if value.lower() in _ZERO_OFFSET_WORDS:
            value = '0'
        targets = ' '.join(lhs_parts[:-1])
        equivalent = f'{mode} {targets} = {value}'
        msg = (
            f'In-line "offset" is not supported: "{self.raw_line}". '
            'Offsets are set in the OFFSETS block. Write '
            f'"{equivalent}" there instead.'
        )
        raise DisplacementsSyntaxError(msg)

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
            return make_direction_token(dir_str)
        except (TokenParserError, ValueError) as err:
            msg = (
                'Unable to parse <direction> token from line in '
                f'{self.block_name} block: {self.raw_line}.'
            )
            raise DisplacementsSyntaxError(msg) from err

    def _parse_range(self, range_str, single_value=False):
        """Return a RangeToken for `range_str`.

        If `single_value` is True, a lone number is read as the fixed
        range <value> <value>.
        """
        if single_value and len(range_str.split()) == 1:
            range_str = f'{range_str.strip()} {range_str.strip()}'
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
                'Unable to parse <mode> information from line in '
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

    def __init__(self, line):
        super().__init__(line)
        self._reject_inline_offset('geo')

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
    `Targets` and `RangeToken` classes, respectively. A single value
    is read as a fixed range, i.e., <value> <value>.
    """

    block_name = 'VIB_DELTA'
    expected_format = '<target> [, <target>] = <range>'

    def __init__(self, line):
        super().__init__(line)
        self._reject_inline_offset('vib')

        # Left hand side
        # check if the last part is a direction
        targets_str, dir_str = separate_direction_from_targets(self._lhs)
        if dir_str:
            raise DisplacementsSyntaxError(self.invalid_format_msg)

        # parse the targets
        self.targets = self._parse_targets(targets_str)

        # Right hand side
        _check_moire_tag(self._rhs)
        # parse to a range (a single value is a fixed one)
        self.range = self._parse_range(self._rhs, single_value=True)

    def _format_lhs_str(self):
        return ', '.join(str(t) for t in self.targets)

    def _format_rhs_str(self):
        return f'{self.range}'


class OccDeltaLine(ParsedLine):
    """Class to parse lines in the OCC_DELTA block of DISPLACEMENTS.

    Lines in the OCC_DELTA block are of the form:
        <target> [, <target>] = <element> <range> [, <element> <range> ...]
    where <target>, <element> and <range> are tokens that are parsed by the
    `Targets`, `Element` and `RangeToken` classes, respectively. A single
    value for an element is read as a fixed range, i.e., <value> <value>.
    """

    block_name = 'OCC_DELTA'
    expected_format = (
        '<target> [, <target>] = <element> <range> [, <element> <range>]'
    )

    def __init__(self, line):
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
        self.element_ranges = element_ranges

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
        range_token = self._parse_range(range_str, single_value=True)
        if len(range_str.split()) == 1 and not 0 <= range_token.start <= 1:
            msg = (
                f'Fixed occupation in line "{self.raw_line}" must be in '
                f'the range [0, 1], got {range_token.start}.'
            )
            raise DisplacementsSyntaxError(msg)
        return self._parse_element(elem_str), range_token


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

    Finally, a parameter can be frozen at a fixed value with
        <mode> <target> [, <target> ...] = <value>
    for geometric and vibrational parameters, and with
        occ <target> [, <target> ...] = <element> <value>
    for occupational ones. For a fixed value, the `link_target` and
    `linear_operation` attributes are None, `fixed_value` is the float
    <value>, and `fixed_element` is the ElementToken of an occupation
    (None for geo and vib). Both are None for all other forms. For
    occupations, a rhs of the form <element> <number> is read as a link
    if it is a valid target (e.g., "O 2"), and as a fixed value
    otherwise (e.g., "O 0.8").

    A parameter can also be constrained to a given index of its
    displacement range with
        <mode> <target> [, <target> ...] = ind(<index>)
    (case-insensitive) for all modes. <index> is a positive integer,
    counted from 1 (i.e., ind(1) is the first value of the range). For
    this form, `link_target`, `linear_operation`, `fixed_value` and
    `fixed_element` are all None.

    Attributes
    ----------
    fixed_index : int or None
        The 1-based index of the displacement range to which the
        parameter is constrained. None for all other forms.
    fixed_value : float or None
        The value at which the parameter is frozen.
    fixed_element : ElementToken or None
        The element whose occupation is frozen.
    """

    block_name = 'CONSTRAIN'
    expected_format = (
        '<mode> <target> [, <target>] = [<linear_operation>] <target>'
    )

    def __init__(self, line):
        super().__init__(line)

        # Left hand side
        self._treat_lhs()
        # Right hand side
        self._treat_rhs()

    def _treat_lhs(self):
        lhs_parts = self._lhs.split()
        if len(lhs_parts) < 2:  # at least mode and one target
            raise DisplacementsSyntaxError(self.invalid_format_msg)

        self.mode_token = self._parse_mode(lhs_parts[0])

        # TODO: Implement Domains here
        if self.mode_token.mode is PerturbationMode.DOM:
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
        self.fixed_value = None
        self.fixed_element = None
        self.fixed_index = None

        # check for an 'offset' tag, which CONSTRAIN never supported
        if 'offset' in self._rhs.lower():
            msg = (
                'CONSTRAIN lines do not support "offset". Offsets are '
                'set in the OFFSETS block.'
            )
            raise DisplacementsSyntaxError(msg)

        # check for 'linked' tag
        if self._rhs.lower().strip() == 'linked':
            logger.log(LOG_VERY_VERBOSE, 'Detected "linked" tag.')
            # "copy" the first target to the rhs
            self.link_target = self.targets[0]
            # treat as if array is identity
            self.linear_operation = LinearOperationToken.from_array(np.eye(1))
            return

        # check for 'total' tag
        if self._rhs.lower().strip().startswith('total '):
            logger.log(LOG_VERY_VERBOSE, 'Detected "total" tag.')

            if self.mode_token.mode is not PerturbationMode.OCC:
                raise DisplacementsSyntaxError(
                    'The "total" tag is only allowed for occupational '
                    'constraints.'
                )

            # parse the total occupation value
            total_occupation_str = self._rhs[len('total ') :].strip()

            # This is a special case: the 'linear operation' is the
            # total occupation itself, and there is no link target.
            self.linear_operation = TotalOccupationToken(total_occupation_str)
            self.link_target = None
            return

        if self._rhs.lower().strip().startswith('ind('):
            self._treat_fixed_index()
            return

        # A number is never a target, but a '-' makes it look like a label
        if self._treat_fixed_value():
            return

        self._treat_link()

    def _treat_fixed_index(self):
        """Read the rhs as ind(<index>) with a positive integer index."""
        match = re.fullmatch(
            r'ind\(\s*(\d+)\s*\)', self._rhs.lower().strip()
        )
        if not match or int(match.group(1)) < 1:
            msg = (
                f'Invalid index constraint in line "{self.raw_line}". '
                'Expected ind(<n>) with n a positive integer.'
            )
            raise DisplacementsSyntaxError(msg)
        self.fixed_index = int(match.group(1))
        self.link_target = None
        self.linear_operation = None

    def _treat_link(self):
        """Read the rhs as [<linear_operation>] <target> or a fixed value."""
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
            except DisplacementsSyntaxError:
                continue
            break
        else:
            self._reject_integer_occupation()
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

    def _treat_fixed_value(self):
        """Read the rhs as <value> (geo, vib) or <element> <value> (occ).

        An occupation given as an integer is not a fixed value, but a
        link to an atom number (e.g., "O 2").

        Returns
        -------
        bool
            Whether the rhs is a fixed value. If so, the fixed-value
            attributes are set.

        Raises
        ------
        DisplacementsSyntaxError
            If the value looks like a number but is not a finite number
            in plain decimal notation, or if a fixed occupation is
            outside [0, 1].
        """
        parts = self._rhs.split()
        is_occ = self.mode_token.mode is PerturbationMode.OCC
        if len(parts) != 1 + is_occ:
            return False
        try:
            float(parts[-1])
        except ValueError:
            return False
        if is_occ and _INTEGER_RE.fullmatch(parts[-1]):
            return False
        value = self._parse_fixed_number(parts[-1])
        if is_occ:
            if not 0 <= value <= 1:
                msg = (
                    f'Fixed occupation in line "{self.raw_line}" must be '
                    f'in the range [0, 1], got {value}.'
                )
                raise DisplacementsSyntaxError(msg)
            self.fixed_element = self._parse_element(parts[0])
        self.fixed_value = value
        self.link_target = None
        self.linear_operation = None
        return True

    def _parse_fixed_number(self, value_str):
        """Return a float from `value_str`, or raise a syntax error."""
        try:
            return parse_plain_float(value_str)
        except ValueError as err:
            msg = (
                f'Invalid value "{value_str}" in line "{self.raw_line}". '
                'Fixed values must be finite numbers in plain decimal '
                'notation (e.g., 0.5).'
            )
            raise DisplacementsSyntaxError(msg) from err

    def _reject_integer_occupation(self):
        """Raise if the rhs is an <element> and an integer non-target."""
        parts = self._rhs.split()
        if (
            self.mode_token.mode is PerturbationMode.OCC
            and len(parts) == _ELEMENT_AND_VALUE
            and _INTEGER_RE.fullmatch(parts[1])
        ):
            msg = (
                f'Invalid line "{self.raw_line}": "{parts[1]}" is read as '
                'an atom number, which is not a valid target. To freeze '
                f'the occupation of {parts[0]}, write a non-integer value '
                f'(e.g., "{parts[0]} {parts[1]}.0").'
            )
            raise DisplacementsSyntaxError(msg)

    def _format_lhs_str(self):
        return ', '.join(str(t) for t in self.targets)

    def _format_rhs_str(self):
        if self.fixed_index is not None:
            return f'ind({self.fixed_index})'
        if self.fixed_value is not None:
            if self.fixed_element is None:
                return f'{self.fixed_value}'
            return f'{self.fixed_element} {self.fixed_value}'
        return f'{self.linear_operation} {self.link_target}'


class OffsetsLine(ParsedLine):
    """Class to parse lines in the OFFSETS block of DISPLACEMENTS.

    Lines in the OFFSETS block are of the form:
        <mode> <target> [, <target> ...] [<direction>] = <offset>
    where <target>, ... are tokens that are parsed by the TargetToken
    class, and <offset> is parsed by the OffsetToken class.
    A direction token must be specified if and only if the mode of the offset
    is geometric.
    """

    block_name = 'OFFSET'
    expected_format = '<mode> <target> [, <target>] [<direction>] = <offset>'

    def __init__(self, line):
        super().__init__(line)
        self.direction = None

        # expected DOF is 1 for vib and occ
        expected_dof = 1

        # parse LHS
        parts = self._lhs.split()
        if len(parts) < 2:  # at least mode and one target
            raise DisplacementsSyntaxError(self.invalid_format_msg)

        self.mode_token = self._parse_mode(parts[0])
        targets_str, dir_str = separate_direction_from_targets(
            ' '.join(parts[1:])
        )

        # parse targets
        self.targets = self._parse_targets(targets_str)

        if self.mode_token.mode is PerturbationMode.GEO:
            # expect and parse direction specifier
            # will raise if no direction is given
            self.direction = self._parse_direction(dir_str)
            # if geometric, expected DOF is given by direction
            expected_dof = self.direction.dof

        if self.mode_token.mode is not PerturbationMode.GEO and dir_str:
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
        if self.mode_token.mode is PerturbationMode.GEO:
            lhs += f' {self.direction}'
        return f'{self.mode_token} {lhs}'

    def _format_rhs_str(self):
        return f'{self.offset}'


def separate_direction_from_targets(targets_and_direction):
    """Separate a string into targets and direction.

    Uses a regex pattern to identify and extract an optional direction
    token at the end of the string. The rest of the string is considered
    the targets.

    Parameters
    ----------
    targets_and_direction : str
        The string to be separated.

    Returns
    -------
    tuple
        A tuple containing the targets string and the direction string.
        Either string may be empty if the corresponding part was not
        found.
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
