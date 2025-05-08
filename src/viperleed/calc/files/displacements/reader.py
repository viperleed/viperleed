"""Module file."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__created__ = '2024-10-04'

import re
from enum import Enum

from .errors import InvalidDisplacementsSyntaxError
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
from .regex import (
    LOOP_END_PATTERN,
    LOOP_START_PATTERN,
    SEARCH_HEADER_PATTERN,
    SECTION_HEADER_PATTERN,
    match_constrain_line,
    match_geo_line,
    match_occ_line,
    match_offsets_line,
    match_vib_line,
)

DisplacementFileSections = Enum(
    'DisplacementFileSections',
    ['OFFSETS', 'GEO_DELTA', 'VIB_DELTA', 'OCC_DELTA', 'CONSTRAIN'],
)

LoopMarker = Enum('LoopMarker', ['LOOP_START', 'LOOP_END'])

from viperleed.calc.files.input_reader import InputFileReader
from viperleed.calc.lib.string_utils import strip_comments


class DisplacementsReader(InputFileReader):
    """Reader for the DISPLACEMENTS file based on InputFileReader."""

    def __init__(self, filename, noisy=True):
        """Initialize instance."""
        super().__init__(filename, noisy)

        self.current_section = None

    def __next__(self):
        """Read the next line from the file."""
        for line in self._file_obj:
            self._current_line += 1
            line_content = self._read_one_line(line)
            if line_content is not None:
                return line_content
        raise StopIteration

    def _read_one_line(self, line):
        """Extract section header or data line from a single line."""
        line = strip_comments(line)
        stripped_line = line.strip().lower()
        if not stripped_line:
            return None

        # check for invalid or deprecated syntax
        _check_line_generally_valid(stripped_line)

        # check for loop markers
        if LOOP_START_PATTERN.match(stripped_line):
            self.current_section = None
            return LoopMarkerLine(LoopMarker.LOOP_START)
        if LOOP_END_PATTERN.match(stripped_line):
            self.current_section = None
            return LoopMarkerLine(LoopMarker.LOOP_END)

        # check for search header
        search_header_match = re.match(SEARCH_HEADER_PATTERN, line)
        if search_header_match:
            label = search_header_match.group(1)
            self.current_section = None
            return SearchHeaderLine(label)

        # check for section headers
        section_header_match = re.match(SECTION_HEADER_PATTERN, line)
        if section_header_match:
            section = section_header_match.group(1)
            self.current_section = section
            return SectionHeaderLine(section)

        # parse section lines
        try:
            return self._parse_line(line)
        except (ValueError, IndexError) as err:
            msg = (
                f'Cannot parse line "{line}" in section '
                f'"{self.current_section}".'
            )
            raise InvalidDisplacementsSyntaxError(msg) from err

    def _parse_line(self, line):
        """Get content from line."""
        # Decide based on the current section
        section = DisplacementFileSections[self.current_section]
        if section is DisplacementFileSections.OFFSETS:
            return self._read_offsets_line(line)
        if section is DisplacementFileSections.GEO_DELTA:
            return self._read_geo_delta_line(line)
        if section is DisplacementFileSections.VIB_DELTA:
            return self._read_vib_delta_line(line)
        if section is DisplacementFileSections.OCC_DELTA:
            return self._read_occ_delta_line(line)
        if section is DisplacementFileSections.CONSTRAIN:
            return self._read_constraints_line(line)
        # If no section is set, raise an error
        msg = f'Cannot parse line "{line}" without a section header.'
        raise ValueError(msg)

    def _read_offsets_line(self, line):
        """Parse a line in the OFFSETS section."""
        match = match_offsets_line(line)
        if match is None:
            msg = f"Cannot parse line '{line}' in OFFSETS section."
            raise InvalidDisplacementsSyntaxError(msg)
        offset_type, targets, direction, value = match
        return OffsetsLine(offset_type, targets, direction, value, line=line)

    def _read_geo_delta_line(self, line):
        """Parse a line in the GEO_DELTA section."""
        match = match_geo_line(line)
        if match is None:
            msg = f"Cannot parse line '{line}' in GEO_DELTA section."
            raise InvalidDisplacementsSyntaxError(msg)

        label, which, direction, start, stop, step = match
        return GeoDeltaLine(label, which, direction, start, stop, step, line=line)

    def _read_vib_delta_line(self, line):
        """Parse a line in the VIB_DELTA section."""
        match = match_vib_line(line)
        if match is None:
            raise InvalidDisplacementsSyntaxError(
                f"Cannot parse line '{line}' in VIB_DELTA section."
            )

        label, which, start, stop, step = match
        return VibDeltaLine(label, which, start, stop, step, line=line)

    def _read_occ_delta_line(self, line):
        """Parse a line in the OCC_DELTA section."""
        match = match_occ_line(line)
        if match is None:
            raise InvalidDisplacementsSyntaxError(
                f"Cannot parse line '{line}' in OCC_DELTA section."
            )

        label, which, chem_blocks = match
        return OccDeltaLine(label, which, chem_blocks, line=line)

    def _read_constraints_line(self, line):
        """Parse a line in the CONSTRAIN section."""
        match = match_constrain_line(line)
        if match is None:
            raise InvalidDisplacementsSyntaxError(
                f"Cannot parse line '{line}' in CONSTRAIN section."
            )

        constraint_type, targets, direction, value = match
        return ConstraintLine(constraint_type, targets, direction, value, line=line)


def _check_line_generally_valid(line):
    """Check the line for invalid or deprecated syntax."""
    # check if the line contains at least one '='
    if '=' not in line:
        msg = f"Could not parse line '{line}'."
        raise InvalidDisplacementsSyntaxError(msg)

    if 'sym_delta' in line.lower():
        msg = ('The SYM_DELTA Tag has been deprecated. Use the SYMMETRY_FIX '
               'parameter instead to manually lower the system symmetry.')
        raise InvalidDisplacementsSyntaxError(msg)
