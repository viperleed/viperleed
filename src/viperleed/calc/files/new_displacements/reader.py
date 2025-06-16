"""Module reader of viperleed.files.displacements."""

__authors__ = ("Alexander M. Imre (@amimre)",)
__copyright__ = "Copyright (c) 2019-2025 ViPErLEED developers"
__created__ = "2024-10-04"
__license__ = "GPLv3+"


import logging
import re
from enum import Enum

from viperleed.calc.files.input_reader import InputFileReader
from viperleed.calc.lib.string_utils import strip_comments

from .errors import DisplacementsSyntaxError
from .lines import (
    ConstraintLine,
    GeoDeltaLine,
    LoopMarkerLine,
    OccDeltaLine,
    OffsetsLine,
    SearchHeaderLine,
    SectionHeaderLine,
    VibDeltaLine,
)

SECTION_HEADER_PATTERN = re.compile(
    r"^=\s*(?P<section>OFFSETS|GEO_DELTA|VIB_DELTA|OCC_DELTA|CONSTRAIN)$"
)
LOOP_START_PATTERN = re.compile(r'<loop>')
LOOP_END_PATTERN = re.compile(r'<\\loop>|</loop>')

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
    structured data. It handles different sections of the file, including
    headers, loop markers, constraints and delta sections (ranges and
    implicit constraints).
    Individual lines are parsed into specific line objects, which are organized
    into sections and search blocks.

    Note that this class only performs the reading and parsing of the file, but
    no interpretation. It does not perform any cross-checking and validation of
    the user input vs. the structure.
    """

    def __init__(self, filename, noisy=True):
        """Initialize instance."""
        super().__init__(filename, noisy)

        self._current_section = None

    def __next__(self):
        """Read the next line from the file."""
        for line in self._file_obj:
            self._current_line += 1
            parsed_line = self._read_one_line(line)
            if isinstance(parsed_line, (LoopMarkerLine, SearchHeaderLine)):
                # search headers and loop markers reset the current section
                self._current_section = None
            elif isinstance(parsed_line, SectionHeaderLine):
                # section headers set the current section
                self._current_section = parsed_line.section
            if parsed_line is not None:
                return parsed_line

        raise StopIteration

    def _read_one_line(self, line):
        """Extract section header or data line from a single line."""
        line = strip_comments(line)
        stripped_line = line.strip().lower()
        if not stripped_line:
            return None

        # check for invalid or deprecated syntax
        _check_line_generally_valid(stripped_line)

        # check for search headers, loop markers and section headers
        for header in (SearchHeaderLine, LoopMarkerLine, SectionHeaderLine):
            try:
                return header(line)
            except DisplacementsSyntaxError:
                continue

        # parse section lines
        if self._current_section is None:
            msg = f'Cannot parse line "{line}" without a section header.'
            raise DisplacementsSyntaxError(msg)

        try:
            # call line parsers
            return DISPLACEMENTS_FILE_SECTION[self._current_section](line)
        except (ValueError, IndexError) as err:
            msg = (
                f'Cannot parse line "{line}" in section '
                f'"{self._current_section}".'
            )
            raise DisplacementsSyntaxError(msg) from err


def _check_line_generally_valid(line):
    """Check the line for invalid or deprecated syntax."""
    # check if the line contains at least one '='
    if '=' not in line:
        msg = f"Could not parse line '{line}'."
        raise DisplacementsSyntaxError(msg)

    if 'sym_delta' in line.lower():
        msg = ('The SYM_DELTA Tag has been deprecated. Use the SYMMETRY_FIX '
               'parameter instead to manually lower the system symmetry.')
        raise DisplacementsSyntaxError(msg)
