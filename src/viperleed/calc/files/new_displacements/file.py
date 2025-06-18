"""Module file of viperleed.calc.files.displacements."""

__authors__ = ("Alexander M. Imre (@amimre)",)
__copyright__ = "Copyright (c) 2019-2025 ViPErLEED developers"
__created__ = "2024-10-04"
__license__ = "GPLv3+"

from .errors import (
    InvalidSearchBlocksError,
    InvalidSearchLoopError,
    OffsetsNotAtBeginningError,
)
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
from .reader import (
    DISPLACEMENTS_FILE_SECTION,
    DisplacementsReader,
    LoopMarker,
)


class SearchBlock:
    """Class to hold all information for a single search block."""

    def __init__(self, label):
        """Initialize with a label and empty sections."""
        self.label = label
        self.sections = {s: [] for s in DISPLACEMENTS_FILE_SECTION}

    def add_line(self, section, line):
        """Add a line to the corresponding section."""
        if section not in self.sections:
            msg = f'Invalid section: {section}'
            raise ValueError(msg)
        self.sections[section].append(line)

    @property
    def geo_delta(self):
        return self.sections['GEO_DELTA']

    @property
    def vib_delta(self):
        return self.sections['VIB_DELTA']

    @property
    def occ_delta(self):
        return self.sections['OCC_DELTA']

    @property
    def explicit_constraints(self):
        return self.sections['CONSTRAIN']

    def __str__(self):
        return f'SearchBlock(label={self.label}, sections={self.sections})'


class OffsetsBlock:
    """Class to hold all information for the (optional) OFFSETS block."""

    def __init__(self):
        self.lines = []

    def add_line(self, line):
        """Add a line to the corresponding section."""
        self.lines.append(line)


class DisplacementsFile:
    def __init__(self):
        self.blocks = []
        self.current_search_block = None
        self.current_section_lines = []
        self.has_been_read = False
        self.has_been_parsed = False
        # an OFFSETS block is only allowed at the very beginning of the file
        # we check this by setting this flag to False after the first block
        self.offsets_block_allowed = True

    def unclosed_loop(self):
        """Return True if there is a loop that was started but not closed."""
        for block in reversed(self.blocks):
            if not isinstance(block, LoopMarkerLine):
                continue
            if block.kind == 'start':
                return True
            if block.kind == 'end':
                break
        return False

    def finish_block(self):
        """Finalize and store the current search block."""
        if self.current_search_block is not None:
            self.blocks.append(self.current_search_block)
            self.current_search_block = None

    def check_valid(self):
        """Check that the read file is valid."""
        # Syntax and logic checks
        # check for unclosed loops
        if self.unclosed_loop():
            raise InvalidSearchLoopError
        if not self.blocks:
            raise InvalidSearchBlocksError
        if any(isinstance(block, OffsetsBlock) for block in self.blocks[1:]):
            raise OffsetsNotAtBeginningError
        return True

    def offsets_block(self):
        """Return the OFFSETS block if present, else None."""
        if isinstance(self.blocks[0], OffsetsBlock):
            return self.blocks[0]
        return None

    def first_block(self):
        """Return the first non loop-marker search block."""
        non_loop_blocks = [
            block for block in self.blocks
            if not isinstance(block, LoopMarkerLine)
        ]
        return non_loop_blocks[0]

    def read(self, filename):
        """Read the file using the DisplacementsReader."""
        if self.has_been_read:
            raise ValueError('read() has already been called.')
        self.has_been_read = True

        with DisplacementsReader(filename, noisy=True) as reader:
            while True:
                try:
                    read = next(reader)
                except StopIteration:
                    self.finish_block()
                    break
                # TODO: low level logging: print(f"Read: {read}")

                if isinstance(read, LoopMarkerLine):
                    if (
                        read.kind == 'start'
                        and self.unclosed_loop()
                    ):
                        raise InvalidSearchLoopError(
                            'Loop started before the previous loop was closed.'
                        )
                    if (
                        read.kind == LoopMarker.LOOP_END
                        and not self.unclosed_loop()
                    ):
                        raise InvalidSearchLoopError(
                            'Loop ended without a matching start.'
                        )
                    self.finish_block()
                    self.blocks.append(read)

                elif isinstance(read, SearchHeaderLine):
                    self.finish_block()
                    # Start a new search block with the label from SearchHeaderLine
                    self.current_search_block = SearchBlock(read.label)

                elif isinstance(read, SectionHeaderLine):
                    if read.section == 'OFFSETS':
                        if not self.offsets_block_allowed:
                            raise OffsetsNotAtBeginningError(
                                'The OFFSETS block is only allowed at the beginning of the file.'
                            )
                        self.finish_block()
                        self.current_search_block = OffsetsBlock()
                    # Update the current section in the active search block
                    self.current_section = read.section

                elif isinstance(
                    read,
                    (GeoDeltaLine, VibDeltaLine, OccDeltaLine, ConstraintLine),
                ):
                    # Add lines to the current section in the active search block
                    if not self.current_search_block:
                        raise ValueError(
                            'No active search block for section line.'
                        )
                    self.current_search_block.add_line(
                        self.current_section, read
                    )

                elif isinstance(read, OffsetsLine):
                    if self.current_section != 'OFFSETS':
                        raise ValueError(
                            'Offsets line found outside of an OFFSETS block.'
                        )
                    self.current_search_block.add_line(self.current_section,
                                                       read)

                else:
                    raise ValueError(f'Unexpected line type: {read}')

                # OFFSETS block is only allowed in the first search block, thus
                # we set the flag to False after the first block
                self.offsets_block_allowed = False

    def parse(self, tl_backend):
        if self.has_been_parsed:
            raise ValueError('parse() has already been called.')
        # After reading, process the blocks
        # There are two steps to this:
        # 1. Check that the general structure is valid, ie. line order, etc.
        # 2. check the requested displacements vs. the backend capabilities
        # Checks about matching the search blocks to the structure are performed
        # when the parameter space is generated.

        # check that the file is valid
        self.check_valid()

        # check against the backend capabilities
        for block in self.blocks:
            new_blocks = []
            if not isinstance(block, SearchBlock):
                new_blocks.append(block)
                continue
            can_handle, replacement = tl_backend.can_handle_search_block(
                self.offsets_block, block
            )
            if can_handle:
                new_blocks.append(block)
                continue
            if replacement is None:
                raise ValueError(
                    f'Tensor LEED backend {tl_backend.name } cannot handle '
                    f'search block {block.label}.'
                )
            # replace the block with the replacement blocks
            new_blocks.extend(replacement)

        # finally, replace the blocks with the new ones
        self.blocks = new_blocks
        self.has_been_parsed = True


class UnknownDisplacementsSegmentError(Exception):
    """Exception raised when a line does not match any known segment header."""
    pass



class DeltaBlock(DisplacementsSegmentABC):
    """Class to hold all information for the DELTA block."""

    required_class_attrs = ['HEADER', 'LINE_TYPE']

    def __init_subclass__(cls, **kwargs):
            super().__init_subclass__(**kwargs)
            for attr in cls.required_class_attrs:
                if not hasattr(cls, attr):
                    raise TypeError(f"Can't instantiate subclass {cls.__name__} without required class attribute '{attr}'")

    @staticmethod
    def is_my_header_line(line):
        """Check if the line is a header for this DELTA block."""
        return isinstance(line, LoopMarker) and line.label == 'DELTA'

    def _belongs_to_me(self, line):
        """Check if the line belongs to this DELTA block."""
        try:
            self.LINE_TYPE(line)
        except AttributeError:
            return False
        return True

    def _parse_block(self):
        """Parse the DELTA block."""
        # Implementation of parsing logic goes here
        self.header_line = SectionHeaderLine(self.HEADER)
        self.lines = [self.LINE_TYPE(line) for line in self._lines]


class GeoDeltaBlock(DeltaBlock):
    HEADER = 'GEO_DELTA'
    LINE_TYPE = GeoDeltaLine

class VibDeltaBlock(DeltaBlock):
    HEADER = 'VIB_DELTA'
    LINE_TYPE = VibDeltaLine

class OccDeltaBlock(DeltaBlock):
    HEADER = 'OCC_DELTA'
    LINE_TYPE = OccDeltaLine

class ConstraintBlock(DeltaBlock):
    HEADER = 'CONSTRAIN'
    LINE_TYPE = ConstraintLine


class UnknownDisplacementsSegmentError(Exception):
    """Exception raised when a line does not match any known segment header."""


class DisplacementsSegmentABC(ABC, NodeMixin):
    """Base class for logical segments in a DISPLACEMENTS file."""

    SUBSEGMENTS = ()  # override in subclasses
    _subclasses = set()

    def __init_subclass__(cls):
        """Register known search segments."""
        cls._subclasses.add(cls)

    def __init__(self, header_line):  # subclasses may use the varargs
        """Initialize instance."""
        self._header = header_line
        self._lines = []

    @classmethod
    def from_header_line(cls, line):
        """Return the a subclass instance for a given header `line`."""
        for subsegment in cls._subclasses:
            if subsegment.is_my_header_line(line):
                return subsegment(line)
        raise UnknownDisplacementsSegmentError(line)

    @staticmethod
    @abstractmethod
    def is_my_header_line(line):
        """Return whether `line` is a header for this segment."""

    @staticmethod
    @abstractmethod
    def _belongs_to_me(line):
        """Return whether `line` is a header for this segment."""

    def read_lines(self, lines):
        while True:
            try:
                line = next(lines)
            except StopIteration:
                return iter([])  # No more lines to read

            processed = False
            print("current line:", line, flush=True)
            if self._belongs_to_me(line):
                self._lines.append(line)
                continue
            # Check if this line starts one of *my* direct subsegments
            for subsegment in self.SUBSEGMENTS:
                if subsegment.is_my_header_line(line):
                    child = subsegment(line)
                    child.parent = self
                    lines = child.read_lines(lines)
                    processed = True
                    break  # after child consumes what it needs, continue outer loop
            if processed:
                print("done: ", line, flush=True)
                continue
            else:
                print("giving up on this line:", line, flush=True)
                return itertools.chain([line], lines)

    def __str__(self):
        return str(RenderTree(self, style=ContStyle()).by_attr("_render_name"))

    @property
    @abstractmethod
    def _render_name(self):
        pass