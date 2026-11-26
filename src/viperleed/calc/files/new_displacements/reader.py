"""Module reader of viperleed.calc.files.new_displacements."""

__authors__ = ('Alexandra Mia Imre (@alexmiame)',)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2024-10-04'
__license__ = 'GPLv3+'


from enum import Enum
import logging

from viperleed.calc.files.input_reader import (
    InputFileReader,
    ShouldSkipLineError,
)
from viperleed.calc.files.new_displacements.errors import (
    DisplacementsSyntaxError,
)
from viperleed.calc.files.new_displacements.lines import (
    LOOP_END_PATTERN,
    LOOP_START_PATTERN,
    ConstraintLine,
    GeoDeltaLine,
    LoopMarkerLine,
    ModeHeaderLine,
    OccDeltaLine,
    OffsetsHeaderLine,
    OffsetsLine,
    SearchHeaderLine,
    VibDeltaLine,
)
from viperleed.calc.lib.string_utils import strip_comments

DISPLACEMENTS_FILE_SECTION = {
    'OFFSETS': OffsetsLine,
    'GEO_DELTA': GeoDeltaLine,
    'VIB_DELTA': VibDeltaLine,
    'OCC_DELTA': OccDeltaLine,
    'CONSTRAIN': ConstraintLine,
}

_ABBREVIATED_HEADERS = {
    'g': '= GEO_DELTA',
    'v': '= VIB_DELTA',
    'o': '= OCC_DELTA',
    'c': '= CONSTRAIN',
}

LoopMarker = Enum('LoopMarker', ['LOOP_START', 'LOOP_END'])

logger = logging.getLogger(__name__)


class DisplacementsReader(InputFileReader):
    """Reader/Parser for the DISPLACEMENTS file based on InputFileReader.

    This class reads the DISPLACEMENTS file and parses its content into
    structured data. It handles the different blocks of the file,
    including headers, loop markers, constraints and delta mode blocks
    (ranges and implicit constraints). Individual lines are parsed into
    specific line objects, which are organized into mode blocks and
    search blocks.

    Note that this class only performs the reading and parsing of the
    file, but no interpretation. It does not perform any cross-checking
    and validation of the user input vs. the structure.
    """

    def __init__(self, filename, noisy=True):
        """Initialize instance."""
        super().__init__(filename, noisy)

        self._current_section = None
        self._follows_loop_start = False

    def _read_one_line(self, line):
        """Extract mode header or data line from a single line."""
        line = strip_comments(line)
        stripped_line = line.strip().lower()
        if not stripped_line:
            raise ShouldSkipLineError

        _check_loop_tags_isolated(line)
        new_header = _read_header(line)
        self._check_follows_loop_start(line, new_header)
        # headers have been validated by their own parsers: only
        # other lines can contain invalid or deprecated syntax
        if new_header is None:
            _check_dropped_header(line)
            _check_line_generally_valid(stripped_line)

        # update current mode if a new header is found
        if isinstance(new_header, (LoopMarkerLine, SearchHeaderLine)):
            # search headers and loop markers reset the current mode
            self._current_section = None
            return new_header
        if isinstance(new_header, (ModeHeaderLine, OffsetsHeaderLine)):
            # mode headers set the current mode
            self._current_section = new_header.section
            return new_header

        # parse mode-block lines
        if self._current_section is None:
            msg = f'Cannot parse line "{line}" without a block header.'
            raise DisplacementsSyntaxError(msg)

        try:
            # call line parsers
            return DISPLACEMENTS_FILE_SECTION[self._current_section](line)
        except (ValueError, KeyError) as err:
            msg = (
                f'Cannot parse line "{line}" in block '
                f'"{self._current_section}".'
            )
            raise DisplacementsSyntaxError(msg) from err


    def _check_follows_loop_start(self, line, header):
        """Complain if `<loop>` is not followed by a search header."""
        if self._follows_loop_start and not _may_follow_loop_start(header):
            msg = (
                f'Found line "{line.strip()}" right after "<loop>". '
                'A "<loop>" must be followed by "== SEARCH <label>" '
                '(or by another "<loop>").'
            )
            raise DisplacementsSyntaxError(msg)
        self._follows_loop_start = (
            isinstance(header, LoopMarkerLine) and header.kind == 'start'
        )


def _may_follow_loop_start(header):
    """Return whether `header` is acceptable right after a <loop>."""
    if isinstance(header, SearchHeaderLine):
        return True
    # An empty loop is reported when the file is validated
    return isinstance(header, LoopMarkerLine)


def _check_loop_tags_isolated(line):
    """Complain if `line` has <loop> tags next to other content."""
    stripped = line.strip()
    if LoopMarkerLine.is_marker(stripped):
        return
    for pattern in (LOOP_START_PATTERN, LOOP_END_PATTERN):
        tag = pattern.search(stripped)
        if tag:
            msg = (
                f'Found "{tag.group(0)}" together with other text in '
                f'line "{stripped}". Loop tags must be on a line of their '
                'own, followed by a "== SEARCH <label>" line.'
            )
            raise DisplacementsSyntaxError(msg)


def _check_dropped_header(line):
    """Complain if `line` is a header written in an unsupported form."""
    stripped = line.strip()
    if not stripped.startswith('='):
        return
    words = stripped.lstrip('=').split()
    if not words:
        return
    word = words[0].lower()
    if stripped.startswith('=='):
        if 'search'.startswith(word):
            msg = (
                f'Cannot parse header "{stripped}": headers cannot be '
                'abbreviated. Write "== SEARCH <label>".'
            )
            raise DisplacementsSyntaxError(msg)
        return
    if word == 'regex':
        msg = (
            f'Cannot parse header "{stripped}": "= REGEX" is not '
            'supported. Use "*" wildcards in the atom labels instead.'
        )
        raise DisplacementsSyntaxError(msg)
    full_name = _ABBREVIATED_HEADERS.get(word[0])
    if word.startswith('off'):
        full_name = '= OFFSETS'
    if full_name:
        msg = (
            f'Cannot parse header "{stripped}": headers cannot be '
            f'abbreviated. Write "{full_name}".'
        )
        raise DisplacementsSyntaxError(msg)


def _read_header(line):
    """Return a header or loop marker read from `line`, or None."""
    for header in (
        LoopMarkerLine,
        OffsetsHeaderLine,
        SearchHeaderLine,
        ModeHeaderLine,
    ):
        try:
            return header(line)
        except DisplacementsSyntaxError:
            continue
    return None


def _check_line_generally_valid(line):
    """Check the line for invalid or deprecated syntax."""
    # check if the line contains at least one '='
    if '=' not in line:
        msg = f"Could not parse line '{line}'."
        raise DisplacementsSyntaxError(msg)

    if 'sym_delta' in line.lower():
        msg = (
            'The SYM_DELTA Tag has been deprecated. Use the SYMMETRY_FIX '
            'parameter instead to manually lower the system symmetry.'
        )
        raise DisplacementsSyntaxError(msg)
