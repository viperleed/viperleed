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

    def _read_one_line(self, line):
        """Extract mode header or data line from a single line."""
        line = strip_comments(line)
        stripped_line = line.strip().lower()
        if not stripped_line:
            raise ShouldSkipLineError

        new_header = _read_header(line)
        # headers have been validated by their own parsers: only
        # other lines can contain invalid or deprecated syntax
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
