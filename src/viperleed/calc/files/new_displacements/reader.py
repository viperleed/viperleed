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
    LOOP_START,
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

# Full names of the headers that open a mode block. In this order,
# the first one a word abbreviates is the one suggested for it.
_MODE_HEADER_NAMES = (
    'GEO_DELTA',
    'VIB_DELTA',
    'OCC_DELTA',
    'CONSTRAIN',
    'OFFSETS',
)
# Spellings that, followed by extra text, are valid headers that are
# not abbreviated
_COMPLETE_MODE_WORDS = {name.lower() for name in _MODE_HEADER_NAMES} | {
    'constraint',
    'offset',
}
_SEARCH_WORD = 'search'
_REGEX_WORD = 'regex'
_SYM_DELTA_WORD = 'sym_delta'

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
        self._pending_markers = []

    def __next__(self):
        """Return the next line, including extra tags of a line."""
        if self._pending_markers:
            return self._pending_markers.pop(0)
        return super().__next__()

    def _read_one_line(self, line):
        """Extract mode header or data line from a single line."""
        line = strip_comments(line)
        stripped_line = line.strip().lower()
        if not stripped_line:
            raise ShouldSkipLineError

        _check_loop_tags_isolated(line)
        tags = LoopMarkerLine.split_tags(line)
        new_header = _read_header(tags[0] if tags else line)
        # headers have been validated by their own parsers: only
        # other lines can contain invalid or deprecated syntax
        if new_header is None:
            _check_dropped_header(line)
        self._check_follows_loop_start(line, new_header)
        for tag in tags[1:]:
            marker = LoopMarkerLine(tag)
            self._check_follows_loop_start(line, marker)
            self._pending_markers.append(marker)
        if new_header is None:
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
            isinstance(header, LoopMarkerLine) and header.kind == LOOP_START
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
    equals = stripped[: len(stripped) - len(stripped.lstrip('='))]
    rest = stripped[len(equals) :]
    words = rest.split()
    if not equals or not words:
        return
    word = words[0].lower()
    if len(equals) > 1:
        _check_search_header(stripped, rest, word)
        return
    if word == _SYM_DELTA_WORD:
        return  # Reported as deprecated
    if _REGEX_WORD.startswith(word) and word.startswith('r'):
        msg = (
            f'Cannot parse header "{stripped}": "= REGEX" is not '
            'supported. Use "*" wildcards in the atom labels instead.'
        )
        raise DisplacementsSyntaxError(msg)
    if word in _COMPLETE_MODE_WORDS:
        msg = (
            f'Cannot parse header "{stripped}": unexpected text after '
            f'header "= {words[0].upper()}".'
        )
        raise DisplacementsSyntaxError(msg)
    for name in _MODE_HEADER_NAMES:
        if name.lower().startswith(word):
            msg = (
                f'Cannot parse header "{stripped}": headers cannot be '
                f'abbreviated. Write "= {name}".'
            )
            raise DisplacementsSyntaxError(msg)
    valid = ', '.join(f'= {name}' for name in _MODE_HEADER_NAMES)
    msg = (
        f'Cannot parse header "{stripped}": unknown header "= '
        f'{words[0]}". Valid headers: {valid}.'
    )
    raise DisplacementsSyntaxError(msg)


def _check_search_header(stripped, rest, word):
    """Complain if `stripped` is a malformed "== SEARCH" header."""
    if _SEARCH_WORD.startswith(word):
        if word != _SEARCH_WORD:
            msg = (
                f'Cannot parse header "{stripped}": headers cannot be '
                'abbreviated. Write "== SEARCH <label>".'
            )
            raise DisplacementsSyntaxError(msg)
        if rest[0].isspace():
            return  # Valid, parsed already
        msg = (
            f'Cannot parse header "{stripped}": a space is needed '
            f'after "==". Write "== {rest.strip()}".'
        )
        raise DisplacementsSyntaxError(msg)
    if word.startswith(_SEARCH_WORD):
        label = rest.strip()[len(_SEARCH_WORD) :].strip()
        msg = (
            f'Cannot parse header "{stripped}": a space is needed '
            'between SEARCH and the label. Write '
            f'"== SEARCH {label}".'
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
