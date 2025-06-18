"""Module file_segments of viperleed.calc.files.displacements."""

__authors__ = ("Alexander M. Imre (@amimre)",)
__copyright__ = "Copyright (c) 2019-2025 ViPErLEED developers"
__created__ = "2025-06-17"
__license__ = "GPLv3+"


from abc import ABC, abstractmethod
import itertools


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
)
from .errors import UnknownDisplacementsSegmentError


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
            + "\n"
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


class LoopBlock(DisplacementsSegmentABC):
    """Class to hold all information for a search block in the DISPLACEMENTS file."""

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

# assign the subsegments for LoopBlock
LoopBlock.SUBSEGMENTS = (SearchBlock, LoopBlock)

class OffsetsBlock(LineContainer):
    """Base class for blocks that contain offsets lines."""

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
        return isinstance(line, SectionHeaderLine) and line.section == 'OFFSETS'

    def _belongs_to_me(self, line):
        """Check if the line belongs to this delta block."""
        return isinstance(line, OffsetsLine)

