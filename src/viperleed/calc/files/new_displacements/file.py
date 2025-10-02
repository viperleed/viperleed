"""Module file of viperleed.calc.files.new_displacements."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2024-10-04'
__license__ = 'GPLv3+'

import logging

from anytree import NodeMixin, RenderTree
from anytree.render import ContStyle

from .errors import (
    AlreadyReadError,
    DisplacementsSyntaxError,
    OffsetsNotAtBeginningError,
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
        # an OFFSETS block is only allowed at the very beginning of the file
        # we check this by setting this flag to True after the first block

        # attributes related to iterating over the search blocks
        self._child_id = 0

    def _check_has_been_read(self):
        """Check if the file has been read."""
        if not self._has_been_read:
            raise ValueError('File has not been read yet. Call read() first.')

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
            while True:
                try:
                    header_line = next(reader)
                except StopIteration:
                    break
                processed = False
                for segment_class in TOP_LEVEL_SEGMENTS:
                    if not segment_class.is_my_header_line(header_line):
                        continue
                    new_segment = segment_class(header_line)
                    reader = new_segment.read_lines(reader)
                    new_segment.parent = self
                    processed = True
                if processed:
                    continue
                # if we reach here, the line was not processed
                msg = f'Unable to parse line: {header_line!r}.'
                raise DisplacementsSyntaxError(msg)

        self._validate()
        logger.info('DISPLACEMENTS file read successfully using new parser.')

    def _validate(self):
        """Validate the file structure."""
        if not self.children:
            raise DisplacementsSyntaxError('The file is empty.')

        # check if there is more than one OFFSETS block
        if (
            len(
                [
                    block
                    for block in self.descendants
                    if isinstance(block, OffsetsBlock)
                ]
            )
            > 1
        ):
            raise OffsetsNotAtBeginningError(
                'There can only be one OFFSETS block at the beginning of the '
                'DISPLACEMENTS file.'
            )

        # Check if the first segment is an OFFSETS block
        if any(
            isinstance(child, OffsetsBlock) for child in self.children
        ) and not isinstance(self.children[0], OffsetsBlock):
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
