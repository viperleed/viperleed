from enum import Enum
from collections import namedtuple
import re

from .errors import InvalidDisplacementsSyntaxError
from .errors import SymmetryViolationError
from .lines import GeoDeltaLine, VibDeltaLine, OccDeltaLine
from .lines import ConstraintLine, OffsetsLine
from .lines import LoopMarkerLine, SearchHeaderLine, SectionHeaderLine
from .regex import match_constrain_line
from .regex import match_geo_line
from .regex import match_occ_line
from .regex import match_offsets_line
from .regex import match_vib_line
from .regex import SEARCH_HEADER_PATTERN
from .regex import SECTION_HEADER_PATTERN

DisplacementFileSections = Enum('DisplacementFileSections', [
    'OFFSETS',
    'GEO_DELTA',
    'VIB_DELTA',
    'OCC_DELTA',
    'CONSTRAIN'
])

LoopMarker = Enum('LoopMarker', [
    'LOOP_START',
    'LOOP_END'
])

from viperleed.calc.files.parameters.file_reader import SettingsFileReader
from viperleed.calc.lib.string_utils import strip_comments

class DisplacementsReader(SettingsFileReader):
    """Reader for the DISPLACEMENTS file based on SettingsFileReader."""

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

        # check for loop markers
        if stripped_line == r"<loop>":
            self.current_section = None
            return LoopMarkerLine(LoopMarker.LOOP_START)
        elif stripped_line == r"<\loop>":
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
        return self._parse_line(line)

    def _parse_line(self, line):
        """Get content from line."""
        # Decide based on the current section
        section = DisplacementFileSections[self.current_section]
        if section is DisplacementFileSections.OFFSETS
            return self._parse_offsets_line(line)
        elif section is DisplacementFileSections.GEO_DELTA:
            return self._parse_geo_delta_line(line)
        elif section is DisplacementFileSections.VIB_DELTA:
            return self._parse_vib_delta_line(line)
        elif section is DisplacementFileSections.OCC_DELTA:
            return self._parse_occ_delta_line(line)
        elif section is DisplacementFileSections.CONSTRAIN:
            return self._parse_constraints_line(line)
        else:
            raise ValueError(
                f"Cannot parse line '{line}' without a section header."
            )

    def _parse_offsets_line(self, line):
        """Parse a line in the OFFSETS section."""
        match = match_offsets_line(line)
        if match is None:
            raise InvalidDisplacementsSyntaxError(f"Cannot parse line '{line}' "
                                     "in OFFSETS section.")
        offset_type, parameters, value = match
        return OffsetsLine(offset_type, parameters, value)

    def _parse_geo_delta_line(self, line):
        """Parse a line in the GEO_DELTA section."""
        match = match_geo_line(line)
        if match is None:
            raise InvalidDisplacementsSyntaxError(f"Cannot parse line '{line}' "
                                     "in GEO_DELTA section.")

        label, which, direction, start, stop, step = match
        return GeoDeltaLine(label, which, direction, start, stop, step)

    def _parse_vib_delta_line(self, line):
        """Parse a line in the VIB_DELTA section."""
        match = match_vib_line(line)
        if match is None:
            raise InvalidDisplacementsSyntaxError(f"Cannot parse line '{line}' "
                                     "in VIB_DELTA section.")

        label, which, start, stop, step = match
        return VibDeltaLine(label, which, start, stop, step)

    def _parse_occ_delta_line(self, line):
        """Parse a line in the OCC_DELTA section."""
        match = match_occ_line(line)
        if match is None:
            raise InvalidDisplacementsSyntaxError(f"Cannot parse line '{line}' "
                                     "in OCC_DELTA section.")

        label, which, chem_blocks = match
        return OccDeltaLine(label, which, chem_blocks)

    def _parse_constraints_line(self, line):
        """Parse a line in the CONSTRAIN section."""
        match = match_constrain_line(line)
        if match is None:
            raise InvalidDisplacementsSyntaxError(f"Cannot parse line '{line}' "
                                     "in CONSTRAIN section.")

        constraint_type, parameters, value = match
        return ConstraintLine(constraint_type, parameters, value)
