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

    @property
    def geo_delta(self):
        return self.sections[DisplacementFileSections.GEO_DELTA]

    @property
    def vib_delta(self):
        return self.sections[DisplacementFileSections.VIB_DELTA]

    @property
    def occ_delta(self):
        return self.sections[DisplacementFileSections.OCC_DELTA]

    @property
    def constrain(self):
        return self.sections[DisplacementFileSections.CONSTRAIN]

    def __repr__(self):
        return f"SearchBlock(label={self.label}, sections={self.sections})"


class DisplacementsFile:

    def __init__(self):
        self.blocks = []
        self.current_search_block = None
        self.current_section_lines = []
        self.has_been_read = False

    @property
    def unclosed_loop(self):
        # Check if there is a loop that was started but not closed
        for block in reversed(self.blocks):
            if block == LOOP_START_MARKER:
                return True
            elif block == LOOP_END_MARKER:
                break
        return False

    def finish_block(self):
        """Finalize and store the current search block."""
        if self.current_search_block is not None:
            self.blocks.append(self.current_search_block)
            self.current_search_block = None


    def read(self, filename):
        """Read the file using the DisplacementsReader."""
        if self.has_been_read:
            raise ValueError("read() has already been called.")
        self.has_been_read = True

        with DisplacementsReader(filename, noisy=True) as reader:
            while True:
                try:
                    read = next(reader)
                except StopIteration:
                    self.finish_block()
                    break

                if isinstance(read, LoopMarkerLine):
                    if read.type == LOOP_START_MARKER and self.unclosed_loop:
                        raise InvalidSearchLoopError(
                            "Loop started before the previous loop was closed."
                        )
                    if read.type == LOOP_END_MARKER and not self.unclosed_loop:
                        raise InvalidSearchLoopError(
                            "Loop ended without a matching start."
                        )
                    self.finish_block()
                    self.blocks.append(read)

                elif isinstance(read, SearchHeaderLine):
                    self.finish_block()
                    # Start a new search block with the label from SearchHeaderLine
                    self.current_search_block = SearchBlock(read.label)

                elif isinstance(read, SectionHeaderLine):
                    # Update the current section in the active search block
                    self.current_section = DisplacementFileSections[read.section]

                elif isinstance(read, (GeoDeltaLine, VibDeltaLine, OccDeltaLine, ConstraintLine)):
                    # Add lines to the current section in the active search block
                    if not self.current_search_block:
                        raise ValueError("No active search block for section line.")
                    self.current_search_block.add_line(self.current_section, read)

                else:
                    raise ValueError("Unexpected line type")

    def parse(self):
        # After reading, process the blocks, e.g., validate the structure or check symmetry.
        pass
