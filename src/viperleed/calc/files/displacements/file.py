"""Module file."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__created__ = '2024-10-04'

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
from .reader import DisplacementFileSections, DisplacementsReader, LoopMarker


class SearchBlock:
    """Class to hold all information for a single search block."""

    def __init__(self, label):
        """Initialize with a label and empty sections."""
        self.label = label
        self.sections = {
            DisplacementFileSections.GEO_DELTA: [],
            DisplacementFileSections.VIB_DELTA: [],
            DisplacementFileSections.OCC_DELTA: [],
            DisplacementFileSections.CONSTRAIN: [],
        }

    def add_line(self, section, line):
        """Add a line to the corresponding section."""
        if section not in self.sections:
            msg = f'Invalid section: {section}'
            raise ValueError(msg)
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
    def explicit_constraints(self):
        return self.sections[DisplacementFileSections.CONSTRAIN]

    def __repr__(self):
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
            if block.type == LoopMarker.LOOP_START:
                return True
            if block.type == LoopMarker.LOOP_END:
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
                        read.type == LoopMarker.LOOP_START
                        and self.unclosed_loop()
                    ):
                        raise InvalidSearchLoopError(
                            'Loop started before the previous loop was closed.'
                        )
                    if (
                        read.type == LoopMarker.LOOP_END
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
                    self.current_section = DisplacementFileSections[
                        read.section
                    ]

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
                    if (
                        self.current_section
                        is not DisplacementFileSections.OFFSETS
                    ):
                        raise ValueError(
                            'Offsets line found outside of an OFFSETS block.'
                        )
                    self.current_search_block.add_line(read)

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
                    f'Tensor LEED backend {tl_backend.name }cannot handle '
                    f'search block  {block.label}.'
                )
            # replace the block with the replacement blocks
            new_blocks.extend(replacement)

        # finally, replace the blocks with the new ones
        self.blocks = new_blocks
        self.has_been_parsed = True
