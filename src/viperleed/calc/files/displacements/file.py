from .errors import InvalidSearchLoopError
from .reader import LOOP_START_MARKER, LOOP_END_MARKER
from .reader import DisplacementsReader
from .reader import DisplacementFileSections
from .lines import GeoDeltaLine, VibDeltaLine, OccDeltaLine, ConstraintLine
from .lines import LoopMarkerLine, SearchHeaderLine, SectionHeaderLine

class SearchBlock:
    """Class to hold all information for a single search block."""

    def __init__(self, label):
        """Initialize with a label and empty sections."""
        self.label = label
        self.sections = {
            DisplacementFileSections.GEO_DELTA: [],
            DisplacementFileSections.VIB_DELTA: [],
            DisplacementFileSections.OCC_DELTA: [],
            DisplacementFileSections.CONSTRAIN: []
        }

    def add_line(self, section, line):
        """Add a line to the corresponding section."""
        if section not in self.sections:
            raise ValueError(f"Invalid section: {section}")
        self.sections[section].append(line)

    def __repr__(self):
        return f"SearchBlock(label={self.label}, sections={self.sections})"

