"""Module file_segments of viperleed.calc.files.displacements."""

__authors__ = ("Alexander M. Imre (@amimre)",)
__copyright__ = "Copyright (c) 2019-2025 ViPErLEED developers"
__created__ = "2025-06-17"
__license__ = "GPLv3+"


from abc import ABC, abstractmethod
import itertools

import numpy as np

from anytree import NodeMixin
from anytree import RenderTree
from anytree.render import ContStyle

from .reader import (
    DISPLACEMENTS_FILE_SECTION,
    DisplacementsReader,
    LoopMarker,
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
    OffsetsHeaderLine,
)
from .errors import UnknownDisplacementsSegmentError, DisplacementsSyntaxError


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
                continue
            else:
                # done reading this segment, return to parent
                self._validate_segment()
                return itertools.chain([line], lines)

    def __str__(self):
        return str(RenderTree(self, style=ContStyle()).by_attr("_render_name"))

    @property
    @abstractmethod
    def _render_name(self):
        pass

    @abstractmethod
    def _validate_segment(self):
        """Check the segments contents run before returning to parent."""

class LineContainer(DisplacementsSegmentABC):
    """Base class for units that contain content lines."""

    SUBSEGMENTS = ()  # Line containers do not allow subsegments

    def __init__(self, header_line):
        self.header = header_line
        self._lines = []

    @abstractmethod
    def _belongs_to_me(self, line):
        """Check if the line belongs to this container."""

    @property
    def _render_name(self):
        sep = "\n\t"
        return (
            f"{self.header}"
            + sep
            + f"{sep.join(line.raw_line for line in self._lines)}"
        )


class DeltaBlock(LineContainer):
    """Base class for blocks that contain delta lines."""

    @property
    @abstractmethod
    def LINE_TYPE(self):
        """Return the line type for this delta block."""

    @property
    @abstractmethod
    def HEADER(self):
        """Return the header for this delta block."""

    @classmethod
    def is_my_header_line(cls, line):
        """Check if the line is a header for this DELTA block."""
        return isinstance(line, SectionHeaderLine) and line.section == cls.HEADER

    def _belongs_to_me(self, line):
        """Check if the line belongs to this delta block."""
        return isinstance(line, self.LINE_TYPE) or self.is_my_header_line(line)

    def _validate_segment(self):
        if not self._lines:
            raise DisplacementsSyntaxError(
                f"Found empty {self.HEADER} block.")


class GeoDeltaBlock(DeltaBlock):
    HEADER = "GEO_DELTA"
    LINE_TYPE = GeoDeltaLine


class VibDeltaBlock(DeltaBlock):
    HEADER = "VIB_DELTA"
    LINE_TYPE = VibDeltaLine


class OccDeltaBlock(DeltaBlock):
    HEADER = "OCC_DELTA"
    LINE_TYPE = OccDeltaLine


class ConstraintBlock(DeltaBlock):
    HEADER = "CONSTRAIN"
    LINE_TYPE = ConstraintLine


class SearchBlock(DisplacementsSegmentABC):
    """Class to hold all information for a search block in the DISPLACEMENTS file."""

    SUBSEGMENTS = (GeoDeltaBlock, VibDeltaBlock, OccDeltaBlock, ConstraintBlock)

    def __init__(self, header_line):
        """Initialize the SearchBlock with a header line."""
        super().__init__(header_line)
        self.label = header_line.label
        self.parent = None  # Parent can be set later if needed
        # self.sections = {s: [] for s in DISPLACEMENTS_FILE_SECTION}

    @staticmethod
    def is_my_header_line(line):
        """Return whether `line` is a header for this segment."""
        return isinstance(line, SearchHeaderLine)

    def _belongs_to_me(self, line):
        """Check if the line belongs to this search block."""
        return False

    @property
    def geo_delta(self):
        return self.sections["GEO_DELTA"]

    @property
    def vib_delta(self):
        return self.sections["VIB_DELTA"]

    @property
    def occ_delta(self):
        return self.sections["OCC_DELTA"]

    @property
    def explicit_constraints(self):
        return self.sections["CONSTRAIN"]

    @property
    def _render_name(self):
        return str(self._header)

    def _validate_segment(self):
        """Check the segments contents run before returning to parent."""
        if not self.children:
            raise DisplacementsSyntaxError(
                f"Empty search block: '{self.label}'."
            )


class LoopBlock(DisplacementsSegmentABC):
    """Class to hold all information for a search block in the DISPLACEMENTS file."""

    def __init__(self, header_line):
        super().__init__(header_line)
    
        # attributes related to iterating over the search blocks
        self._last_rfac = np.inf  # start with large value; loop will always
        self._current_child_id = 0

    @staticmethod
    def is_my_header_line(line):
        """Loops are started by Loop start lines."""
        return isinstance(line, LoopMarkerLine) and line.kind == 'start'

    def _belongs_to_me(self, line):
        """Loop end lines belong to this block."""
        return isinstance(line, LoopMarkerLine) and line.kind == 'end'

    @property
    def _render_name(self):
        return 'Loop Block'

    def _validate_segment(self):
        # loop block must contain exactly one line – the end loop marker
        if len(self._lines) != 1:
            raise DisplacementsSyntaxError(
                "Unfinished loop block."
            )
        # it must contain at least one search block
        if not self.children:
            raise DisplacementsSyntaxError(
                "Loop block must contain at least one search block."
            )
        # a loop cannot contain only one other loop block
        if len(self.children) == 1 and isinstance(self.children[0], LoopBlock):
            raise DisplacementsSyntaxError(
                "Loop block cannot contain only another loop block."
            )

    def next(self, current_rfac, r_fac_eps=1e-4):
        """Check convergence criteria and return the next search block."""

        if self._current_child_id >= len(self.children):
            # if all children have been processed, check convergence
            # specifically do not use abs() here – break if we increased
            if self._last_rfac - current_rfac <= r_fac_eps:
                # if converged, reset the child index and call next() again
                self._current_child_id = 0
                self._last_rfac = current_rfac
                return self.next(current_rfac, r_fac_eps)
            raise StopIteration

        next_block = self.children[self._current_child_id]
        # if the current child is a loop block, delegate to it
        if isinstance(next_block, LoopBlock):
            # call the loop block's next method to get the next search block
            try:
                return self.children[self.child_id].next(current_rfac)
            except StopIteration:
                # if the loop block is exhausted, move to the next child
                self.child_id += 1
                return self.next(current_rfac, r_fac_eps=r_fac_eps)

        # otherwise, return the next search block
        self._current_child_id += 1
        return next_block


# assign the subsegments for LoopBlock
LoopBlock.SUBSEGMENTS = (SearchBlock, LoopBlock)

class OffsetsBlock(LineContainer):
    """Base class for blocks that contain offsets lines."""

    @classmethod
    def is_my_header_line(cls, line):
        """Check if the line is a header for this DELTA block."""
        return isinstance(line, OffsetsHeaderLine)

    def _belongs_to_me(self, line):
        """Check if the line belongs to this delta block."""
        return isinstance(line, OffsetsLine)

    def _validate_segment(self):
        if not self._lines:
            raise DisplacementsSyntaxError("Empty OFFSETS block.")
