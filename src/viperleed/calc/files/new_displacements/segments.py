"""Module segments of viperleed.calc.files.new_displacements."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__copyright__ = 'Copyright (c) 2019-2025 ViPErLEED developers'
__created__ = '2025-06-17'
__license__ = 'GPLv3+'


from abc import ABC, abstractmethod
import itertools

from anytree import NodeMixin, RenderTree
from anytree.render import ContStyle
import numpy as np

from .errors import DisplacementsSyntaxError, UnknownDisplacementsSegmentError
from .lines import (
    ConstraintLine,
    GeoDeltaLine,
    LoopMarkerLine,
    OccDeltaLine,
    OffsetsHeaderLine,
    OffsetsLine,
    SearchHeaderLine,
    SectionHeaderLine,
    VibDeltaLine,
)


class DisplacementsSegmentABC(ABC, NodeMixin):
    """Base class for logical segments in a DISPLACEMENTS file."""

    _subclasses = set()

    def __init_subclass__(cls):
        """Register known search segments."""
        cls._subclasses.add(cls)

    def __init__(self, header_line):
        """Initialize instance."""
        self._header = header_line
        self._lines = []
        self._subsegments = ()
        super().__init__()

    @property
    def subsegments(self):
        """Return the possible subsegments of this segment."""
        return self._subsegments

    @classmethod
    def from_header_line(cls, line):
        """Return a subclass instance for a given header `line`."""
        for subsegment in cls._subclasses:
            if subsegment.is_my_header_line(line):
                return subsegment(line)
        raise UnknownDisplacementsSegmentError(line)

    @staticmethod
    @abstractmethod
    def is_my_header_line(line):
        """Return whether `line` is a header for this segment."""

    @property
    def lines(self):
        """Return the lines belonging to this segment."""
        return tuple(self._lines)

    def read_lines(self, lines):
        while True:
            try:
                line = next(lines)
            except StopIteration:
                # No more lines to read; validate and return
                self.validate_segment()
                return iter([])

            if self._belongs_to_me(line):
                self._lines.append(line)
                continue
            # Check if this line starts one of *my* direct subsegments
            for subsegment in self.subsegments:
                if subsegment.is_my_header_line(line):
                    child = subsegment(line)
                    child.parent = self
                    lines = child.read_lines(lines)
                    # continue outer loop
                    # after child consumes what it needs
                    break
            # done reading this segment, return to parent
            self.validate_segment()
            return itertools.chain([line], lines)

    def __str__(self):
        return str(RenderTree(self, style=ContStyle()).by_attr('_render_name'))

    @property
    @abstractmethod
    def _render_name(self):
        pass

    @abstractmethod
    def validate_segment(self):
        """Check the segments contents run before returning to parent."""

    @abstractmethod
    def _belongs_to_me(self, line):
        """Return whether `line` is a header for this segment."""


class LineContainer(DisplacementsSegmentABC):
    """Base class for units that contain content lines."""

    header = None

    def __init_subclass__(cls, **kwargs):
        """Verify the presence of expected class attributes."""
        super().__init_subclass__(**kwargs)
        if not cls.header:
            msg = f"{cls.__name__} must define a 'header' class attribute"
            raise TypeError(msg)

    @property
    def _render_name(self):
        sep = '\n\t'
        return (
            f'{self.header}'
            + sep
            + f'{sep.join(line.raw_line for line in self._lines)}'
        )


class DeltaBlock(LineContainer):
    """Base class for blocks that contain delta lines."""

    line_type = None

    def __init_subclass__(cls, **kwargs):
        """Verify the presence of expected class attributes."""
        super().__init_subclass__(**kwargs)
        if not cls.line_type:
            msg = f"{cls.__name__} must define a 'line_type' class attribute"
            raise TypeError(msg)

    @classmethod
    def is_my_header_line(cls, line):
        """Check if the line is a header for this DELTA block."""
        return (
            isinstance(line, SectionHeaderLine) and line.section == cls.header
        )

    def _belongs_to_me(self, line):
        """Check if the line belongs to this delta block."""
        return isinstance(line, self.line_type) or self.is_my_header_line(line)

    def validate_segment(self):
        if not self._lines:
            raise DisplacementsSyntaxError(f'Found empty {self.header} block.')


class GeoDeltaBlock(DeltaBlock):
    """Class to hold information about a GEO_DELTA block in DISPLACEMENTS."""

    header = 'GEO_DELTA'
    line_type = GeoDeltaLine


class VibDeltaBlock(DeltaBlock):
    """Class to hold information about a VIB_DELTA block in DISPLACEMENTS."""
    header = 'VIB_DELTA'
    line_type = VibDeltaLine


class OccDeltaBlock(DeltaBlock):
    """Class to hold information about a OCC_DELTA block in DISPLACEMENTS."""
    header = 'OCC_DELTA'
    line_type = OccDeltaLine


class ConstraintBlock(DeltaBlock):
    """Class to hold information about a CONSTRAIN block in DISPLACEMENTS."""
    header = 'CONSTRAIN'
    line_type = ConstraintLine


class SearchBlock(DisplacementsSegmentABC):
    """Class to hold all information for a search block in the DISPLACEMENTS file."""

    subsegments = (
        GeoDeltaBlock,
        VibDeltaBlock,
        OccDeltaBlock,
        ConstraintBlock,
    )

    def __init__(self, header_line):
        """Initialize the SearchBlock with a header line."""
        super().__init__(header_line)
        self.label = header_line.label
        self.parent = None  # Parent can be set later if needed

    @staticmethod
    def is_my_header_line(line):
        """Return whether `line` is a header for this segment."""
        return isinstance(line, SearchHeaderLine)

    def _get_block_lines(self, block_type):
        """Return the first block of the given type."""
        for child in self.children:
            if isinstance(child, block_type):
                return child.lines
        return []

    @property
    def geo_delta_lines(self):
        """Return lines of the GeoDeltaBlock."""
        return self._get_block_lines(GeoDeltaBlock)

    @property
    def vib_delta_lines(self):
        """Return lines of the VibDeltaBlock."""
        return self._get_block_lines(VibDeltaBlock)

    @property
    def occ_delta_lines(self):
        """Return lines of the OccDeltaBlock."""
        return self._get_block_lines(OccDeltaBlock)

    @property
    def explicit_constraint_lines(self):
        """Return lines of the ConstraintBlock."""
        return self._get_block_lines(ConstraintBlock)

    def validate_segment(self):
        """Check the segments contents run before returning to parent."""
        if not self.children:
            msg = f'Empty search block: {self.label!r}.'
            raise DisplacementsSyntaxError(msg)

        # check that there is at most one block of each type
        for block_type in self.subsegments:
            blocks = [
                child
                for child in self.children
                if isinstance(child, block_type)
            ]
            if len(blocks) > 1:
                msg = (
                    f'Multiple {block_type.header} blocks found in search '
                    f"block: '{self.label}'."
                )
                raise DisplacementsSyntaxError(msg)

    def _belongs_to_me(self, line):
        """Check if the line belongs to this search block."""
        return False

    @property
    def _render_name(self):
        return str(self._header)


class LoopBlock(DisplacementsSegmentABC):
    """Class holding information about a loop block in DISPLACEMENTS."""

    def __init__(self, header_line):
        super().__init__(header_line)

        # attributes related to iterating over the search blocks
        # start with large value; loop will always decrease it
        self._last_rfac = np.inf
        self._current_child_id = 0
        self._subsegments = (SearchBlock, LoopBlock)

    @staticmethod
    def is_my_header_line(line):
        """Loops are started by Loop start lines."""
        return isinstance(line, LoopMarkerLine) and line.kind == 'start'

    def _belongs_to_me(self, line):
        """Loop end lines belong to this block."""
        # is the loop already closed?
        if self.lines and self.lines[-1].kind == 'end':
            return False
        # if not, check if this is the end line
        is_end_line = isinstance(line, LoopMarkerLine) and line.kind == 'end'
        if is_end_line:
            # accept no more subsegments after this
            self._subsegments = ()
        return is_end_line

    @property
    def _render_name(self):
        return 'Loop Block'

    def validate_segment(self):
        # loop block must contain exactly one line – the end loop marker
        if len(self.lines) < 1:
            raise DisplacementsSyntaxError('Unfinished loop block.')
        # it must contain at least one search block
        if not self.children:
            raise DisplacementsSyntaxError('Loop block cannot be empty.')
        # a loop cannot contain only one other loop block
        if len(self.children) == 1 and isinstance(self.children[0], LoopBlock):
            raise DisplacementsSyntaxError(
                'Loop block cannot contain only another loop block.'
            )

    def next(self, current_rfac, r_fac_eps=1e-4):
        """Check convergence criteria and return the next search block."""
        if self._current_child_id >= len(self.children):
            # if all children have been processed, check convergence
            # specifically do not use abs() here - break if we increased
            if self._last_rfac - current_rfac >= r_fac_eps:
                # if converged, reset the child index and call next()
                self._current_child_id = 0
                self._last_rfac = current_rfac
                return self.next(current_rfac, r_fac_eps)
            raise StopIteration

        next_block = self.children[self._current_child_id]
        # if the current child is a loop block, delegate to it
        if isinstance(next_block, LoopBlock):
            # call the loop block's next method to get the next search block
            try:
                return self.children[self._current_child_id].next(current_rfac)
            except StopIteration:
                # if the loop block is exhausted, move to the next child
                self._current_child_id += 1
                return self.next(current_rfac, r_fac_eps=r_fac_eps)

        # otherwise, return the next search block
        self._current_child_id += 1
        return next_block


class OffsetsBlock(LineContainer):
    """Base class for blocks that contain offsets lines."""

    @classmethod
    def is_my_header_line(cls, line):
        """Check if the line is a header for this DELTA block."""
        return isinstance(line, OffsetsHeaderLine)

    def _belongs_to_me(self, line):
        """Check if the line belongs to this delta block."""
        return isinstance(line, OffsetsLine)

    def validate_segment(self):
        """Validate that the offsets block is not empty."""
        if not self._lines:
            raise DisplacementsSyntaxError('Empty OFFSETS block.')
