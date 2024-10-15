from collections import namedtuple


from .targeting import BSTarget
from .direction import Direction


LoopMarkerLine = namedtuple("LoopMarkerLine", ["type"])
SearchHeaderLine = namedtuple("SearchHeaderLine", ["label"])
SectionHeaderLine = namedtuple("SectionHeaderLine", ["section"])


def _get_target(label, which):
    if which is None:
        return BSTarget(label)
    else:
        return BSTarget(f"{label} {which}")


class GeoDeltaLine:
    def __init__(self, label, which, direction, start, stop, step):
        self.target = _get_target(label, which)
        self.direction = Direction(direction)
        self.start = start
        self.stop = stop
        self.step = step

    def __eq__(self, other):
        if isinstance(other, GeoDeltaLine):
            return (
                self.target == other.target
                and self.direction == other.direction
                and self.start == other.start
                and self.stop == other.stop
                and self.step == other.step
            )
        return False

    def __repr__(self):
        return (
            f"GeoDeltaLine(target={self.target}, "
            f"direction={self.direction}, start={self.start}, stop={self.stop}, step={self.step})"
        )


class VibDeltaLine:
    def __init__(self, label, which, start, stop, step):
        self.target = _get_target(label, which)
        self.start = start
        self.stop = stop
        self.step = step

    def __eq__(self, other):
        if isinstance(other, VibDeltaLine):
            return (
                self.target == other.target
                and self.direction == other.direction
                and self.start == other.start
                and self.stop == other.stop
                and self.step == other.step
            )
        return False

    def __repr__(self):
        return (
            f"VibDeltaLine(target={self.target}, "
            f"start={self.start}, stop={self.stop}, step={self.step})"
        )


class OccDeltaLine:
    def __init__(self, label, which, chem_blocks):
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
        return (
            f"OccDeltaLine(target={self.target}, "
            f"chem_blocks={self.chem_blocks})"
        )


class ConstraintLine:
    def __init__(self, constraint_type, parameters, value):
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
        return (
            f"ConstraintLine(constraint_type={self.constraint_type}, "
            f"parameters={self.parameters}, value={self.value})"
        )


class OffsetsLine:
    def __init__(self, offset_type, parameters, value):
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
        return (
            f"ConstraintLine(constraint_type={self.offset_type}, "
            f"parameters={self.parameters}, value={self.value})"
        )
