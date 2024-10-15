from collections import namedtuple


from .targeting import BSTarget
from .direction import Direction
from .range import DisplacementsRange


LoopMarkerLine = namedtuple("LoopMarkerLine", ["type"])
SearchHeaderLine = namedtuple("SearchHeaderLine", ["label"])
SectionHeaderLine = namedtuple("SectionHeaderLine", ["section"])


def _get_target(label, which):
    if which is None:
        return BSTarget(label)
    else:
        return BSTarget(f"{label} {which}")


class GeoDeltaLine:
    def __init__(self, line, label, which, direction, start, stop, step):
        self.line = line
        self.target = _get_target(label, which)
        self.direction = Direction(direction)
        self.range = DisplacementsRange(start, stop, step)

    def __eq__(self, other):
        if isinstance(other, GeoDeltaLine):
            return (
                self.target == other.target
                and self.direction == other.direction
                and self.range == other.range
            )
        return False

    def __repr__(self):
        return self.line

class VibDeltaLine:
    def __init__(self, line, label, which, start, stop, step):
        self.line = line
        self.target = _get_target(label, which)
        self.range = DisplacementsRange(start, stop, step)

    def __eq__(self, other):
        if isinstance(other, VibDeltaLine):
            return (
                self.target == other.target
                and self.direction == other.direction
                and self.range == other.range
            )
        return False

    def __repr__(self):
        return self.line


class OccDeltaLine:
    def __init__(self, line, label, which, chem_blocks):
        self.line = line
        self.target = _get_target(label, which)
        self.chem_blocks = chem_blocks

    def __eq__(self, other):
        if isinstance(other, OccDeltaLine):
            return (
                self.target == other.target
                and self.chem_blocks == other.chem_blocks
            )
        return False

    def __repr__(self):
        return self.line


class ConstraintLine:
    def __init__(self, line, constraint_type, targets, direction, value):
        self.line = line
        self.constraint_type = constraint_type
        self.parameters = parameters
        self.value = value

    def __eq__(self, other):
        if isinstance(other, ConstraintLine):
            return (
                self.constraint_type == other.constraint_type
                and self.parameters == other.parameters
                and self.value == other.value
            )
        return False

    def __repr__(self):
        return self.line


class OffsetsLine:
    def __init__(self, line, offset_type, parameters, value):
        self.line = line
        self.offset_type = offset_type
        self.parameters = parameters
        self.value = value

    def __eq__(self, other):
        if isinstance(other, ConstraintLine):
            return (
                self.offset_type == other.constraint_type
                and self.parameters == other.parameters
                and self.value == other.value
            )
        return False

    def __repr__(self):
        return self.line
