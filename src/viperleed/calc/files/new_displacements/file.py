"""Module file of viperleed.calc.files.new_displacements."""

__authors__ = ('Alexandra Mia Imre (@alexmiame)',)
__copyright__ = 'Copyright (c) 2019-2026 ViPErLEED developers'
__created__ = '2024-10-04'
__license__ = 'GPLv3+'

import itertools
import logging

from anytree import (
    NodeMixin,
    PreOrderIter,
    RenderTree,
)
from anytree.render import ContStyle

from .errors import (
    AlreadyReadError,
    DisplacementsSyntaxError,
    OffsetsNotAtBeginningError,
)
from .lines import (
    LOOP_END,
    LoopMarkerLine,
    ModeHeaderLine,
    SearchHeaderLine,
)
from .reader import DisplacementsReader
from .segments import LoopBlock, OffsetsBlock, SearchBlock

logger = logging.getLogger(__name__)

TOP_LEVEL_SEGMENTS = (
    OffsetsBlock,
    SearchBlock,
    LoopBlock,
)


class DisplacementsFile(NodeMixin):
    """Class representing the DISPLACEMENTS file.

    This class is responsible for reading and parsing the DISPLACEMENTS file,
    validating its structure, and providing an interface to access the parsed
    segments.
    Use the `read` method to read the file and parse its content.
    After reading, you can access the next search block using the `next`
    method, which will return the next search block or raise StopIteration if
    there are no more blocks to process.

    Note
    ----
    Since loop blocks need to check convergence criteria, the next() method
    requires the current R-factor to be passed as an argument. This means that
    this class cannot be used as a standard Python iterator. Instead, it
    provides a `next` method that takes the current R-factor as an argument
    and returns the next search block.
    """

    def __init__(self):
        super().__init__()
        self.name = self._render_name
        self._has_been_read = False

        # attributes related to iterating over the search blocks
        self._child_id = 0

    def _check_has_been_read(self):
        """Check if the file has been read."""
        if not self._has_been_read:
            raise ValueError('File has not been read yet. Call read() first.')

    def reset_iterator(self):
        """Reset to the beginning of the DISPLACEMENTS file.

        This is needed, e.g., when the search repeats after a new
        reference calculation.
        """
        self._child_id = 0
        for child in self.children:
            if isinstance(child, LoopBlock):
                child.reset_iterator()

    @property
    def offsets(self):
        """Return the OFFSETS block if present, else None."""
        self._check_has_been_read()
        if isinstance(self.children[0], OffsetsBlock):
            return self.children[0]
        return None

    def read(self, filename):
        """Read the file using the DisplacementsReader.

        Reads the DISPLACEMENTS file and parses its content into
        structured data using the DisplacementsReader.
        Since the DISPLACEMENTS file is not supposed to change during a
        run, calling this method more than once will raise an
        AlreadyReadError.

        Parameters
        ----------
        filename : path-like
            The path to the DISPLACEMENTS file to be read.

        Raises
        ------
        AlreadyReadError
            If the file has already been read.
        DisplacementsSyntaxError
            If the file cannot be parsed or has an invalid structure.
        OffsetsNotAtBeginningError
            If a OFFSETS block is found to not be at the beginning of the file.
        """
        if self._has_been_read:
            raise AlreadyReadError('read() has already been called.')
        self._has_been_read = True

        with DisplacementsReader(filename) as reader:
            lines = itertools.chain(reader)
            while True:
                try:
                    header_line = next(lines)
                except StopIteration:
                    break
                if isinstance(header_line, ModeHeaderLine):
                    # Mode blocks before the first search header form
                    # an implicit search block
                    if self._has_search_blocks():
                        msg = (
                            f'Line "{header_line}" is not inside a search '
                            'block. Start a new one with '
                            '"== SEARCH <label>" before it.'
                        )
                        raise DisplacementsSyntaxError(msg)
                    lines = itertools.chain([header_line], lines)
                    header_line = SearchHeaderLine('== SEARCH')
                for segment_class in TOP_LEVEL_SEGMENTS:
                    if not segment_class.is_my_header_line(header_line):
                        continue
                    new_segment = segment_class(header_line)
                    lines = new_segment.read_lines(lines)
                    new_segment.parent = self
                    break
                else:
                    raise DisplacementsSyntaxError(
                        _unparsable_line_message(header_line)
                    )

        self._name_unnamed_search_blocks()
        self._validate()
        logger.info('DISPLACEMENTS file read successfully using new parser.')

    def _has_search_blocks(self):
        """Return whether any search block has been read so far."""
        return any(isinstance(b, SearchBlock) for b in self.descendants)

    def _name_unnamed_search_blocks(self):
        """Give every search block a unique label.

        Labels given by the user are kept, except that a repeated one
        gets the first free suffix ' 2', ' 3', ... (as the old reader
        does). Blocks without a label are named search_<position>, with
        <position> their 1-based index among all search blocks. If that
        name is taken, <position> is increased until it is free.
        """
        blocks = [b for b in PreOrderIter(self) if isinstance(b, SearchBlock)]
        used = set()
        for block in blocks:
            if not block.label:
                continue
            label, suffix = block.label, 2
            while label in used:
                label = f'{block.label} {suffix}'
                suffix += 1
            block.label = label
            used.add(label)
        for position, block in enumerate(blocks, start=1):
            if block.label:
                continue
            number = position
            while f'search_{number}' in used:
                number += 1
            block.label = f'search_{number}'
            used.add(block.label)

    def _validate(self):
        """Validate the file structure."""
        if not self.children:
            raise DisplacementsSyntaxError('The file is empty.')

        offset_blocks = [
            b for b in self.descendants if isinstance(b, OffsetsBlock)
        ]
        if len(offset_blocks) > 1:
            raise OffsetsNotAtBeginningError(
                'There can only be one OFFSETS block at the beginning of the '
                'DISPLACEMENTS file.'
            )

        # Check if the first segment is an OFFSETS block
        if offset_blocks and offset_blocks[0] is not self.children[0]:
            raise OffsetsNotAtBeginningError(
                'The OFFSETS block must be at the beginning of the file.'
            )

        # Validate each segment
        for child in self.children:
            child.validate_segment()

    @property
    def _render_name(self):
        return 'DISPLACEMENTS'

    def __str__(self):
        """Return a string representation of the DISPLACEMENTS file."""
        return str(RenderTree(self, style=ContStyle()).by_attr('_render_name'))

    def next(self, current_rfac, r_fac_eps=1e-4):
        """Return the next search to be executed or raise StopIteration.

        Parameters
        ----------
        current_rfac : float
            The current R-factor value to check convergence criteria.
        r_fac_eps : float, optional
            The tolerance for the R-factor convergence check (default is 1e-4).

        Returns
        -------
        SearchBlock
            The next search block to be executed if there is one available.

        Raises
        ------
        StopIteration
            If there are no more search blocks to process.
        ValueError
            If called before the file has been read.
        """
        self._check_has_been_read()

        if self._child_id == 0 and isinstance(self.children[0], OffsetsBlock):
            # if the first child is an OFFSETS block, ignore it
            self._child_id += 1

        if self._child_id >= len(self.children):
            raise StopIteration

        # if the current child is a loop block, handle it
        if isinstance(self.children[self._child_id], LoopBlock):
            # call the loop block's next method to get the next search block
            try:
                return self.children[self._child_id].next(
                    current_rfac, r_fac_eps=r_fac_eps
                )
            except StopIteration:
                # if the loop block is exhausted, move to the next child
                self._child_id += 1
                return self.next(current_rfac, r_fac_eps=r_fac_eps)

        # otherwise, return the next search block
        search_block = self.children[self._child_id]
        self._child_id += 1
        return search_block


def _unparsable_line_message(line):
    """Return an error message for a `line` found at top level."""
    if isinstance(line, LoopMarkerLine) and line.kind == LOOP_END:
        return f'Found "{line}" without a matching "<loop>".'
    return f'Unable to parse line "{line}".'
