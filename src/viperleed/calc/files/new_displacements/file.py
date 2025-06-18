"""Module file of viperleed.calc.files.displacements."""

__authors__ = ("Alexander M. Imre (@amimre)",)
__copyright__ = "Copyright (c) 2019-2025 ViPErLEED developers"
__created__ = "2024-10-04"
__license__ = "GPLv3+"

import logging

from anytree import NodeMixin, RenderTree
from anytree.render import ContStyle

from .errors import DisplacementsSyntaxError, OffsetsNotAtBeginningError
from .reader import DisplacementsReader
from .segments import OffsetsBlock, SearchBlock, LoopBlock

logger = logging.getLogger(__name__)

TOP_LEVEL_SEGMENTS = (
    OffsetsBlock,
    SearchBlock,
    LoopBlock,
)

class DisplacementsFile(NodeMixin):
    def __init__(self):
        self.name = self._render_name
        self.has_been_read = False
        #self.has_been_parsed = False
        # an OFFSETS block is only allowed at the very beginning of the file
        # we check this by setting this flag to False after the first block

    def offsets_block(self):
        """Return the OFFSETS block if present, else None."""
        if isinstance(self.children[0], OffsetsBlock):
            return self.children[0]
        return None

    def read(self, filename):
        """Read the file using the DisplacementsReader."""
        if self.has_been_read:
            raise ValueError('read() has already been called.')
        self.has_been_read = True

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
                    if segment_class is OffsetsBlock:
                        pass
                    new_segment = segment_class(header_line)
                    reader = new_segment.read_lines(reader)
                    new_segment.parent = self
                    processed = True
                if processed:
                    logger.debug('Parsed line: ' + str(new_segment))
                    continue
                # if we reach here, the line was not processed
                raise DisplacementsSyntaxError(
                    f'Unable to parse line: {header_line!r}.')

        self._validate()
        logger.info('DISPLACEMENTS file read successfully using new parser.')

    def _validate(self):
        """Validate the file structure."""
        if not self.children:
            raise DisplacementsSyntaxError('The file is empty.')

        # check if there is more than one OFFSETS block
        if len(list(block for block in self.descendants
                    if isinstance(block, OffsetsBlock))) > 1:
            raise OffsetsNotAtBeginningError(
                "There can only be one OFFSETS block at the beginning of the file."
            )

        # Check if the first segment is an OFFSETS block
        if any(isinstance(child, OffsetsBlock) for child in self.children):
            if not isinstance(self.children[0], OffsetsBlock):
                raise OffsetsNotAtBeginningError(
                    "The OFFSETS block must be at the beginning of the file."
                )

        # Validate each segment
        for child in self.children:
            child._validate_segment()

    @property
    def _render_name(self):
        return 'DISPLACEMENTS'

    def __str__(self):
        return str(RenderTree(self, style=ContStyle()).by_attr("_render_name"))
