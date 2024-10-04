from collections import namedtuple

LoopMarkerLine = namedtuple('LoopMarkerLine', ['type'])
SearchHeaderLine = namedtuple('SearchHeaderLine', ['label'])
SectionHeaderLine = namedtuple('SectionHeaderLine', ['section'])

class GeoDeltaLine:
    def __init__(self, label, which, direction, start, stop, step):
        self.label = label
        self.which = which
        self.direction = direction
        self.start = start
        self.stop = stop
        self.step = step

    def __eq__(self, other):
        if isinstance(other, GeoDeltaLine):
            return (self.label == other.label and
                    self.which == other.which and
                    self.direction == other.direction and
                    self.start == other.start and
                    self.stop == other.stop and
                    self.step == other.step)
        return False

    def __repr__(self):
        return (f"GeoDeltaLine(label={self.label}, which={self.which}, "
                f"direction={self.direction}, start={self.start}, stop={self.stop}, step={self.step})")

class VibDeltaLine:
    def __init__(self, label, which, start, stop, step):
        self.label = label
        self.which = which
        self.start = start
        self.stop = stop
        self.step = step

    def __eq__(self, other):
        if isinstance(other, VibDeltaLine):
            return (self.label == other.label and
                    self.which == other.which and
                    self.start == other.start and
                    self.stop == other.stop and
                    self.step == other.step)
        return False

    def __repr__(self):
        return (f"VibDeltaLine(label={self.label}, which={self.which}, "
                f"start={self.start}, stop={self.stop}, step={self.step})")

class OccDeltaLine:
    def __init__(self, label, which, chem_blocks):
        self.label = label
        self.which = which
        self.chem_blocks = chem_blocks

    def __eq__(self, other):
        if isinstance(other, OccDeltaLine):
            return (self.label == other.label and
                    self.which == other.which and
                    self.chem_blocks == other.chem_blocks)
        return False

    def __repr__(self):
        return (f"OccDeltaLine(label={self.label}, which={self.which}, "
                f"chem_blocks={self.chem_blocks})")

class ConstraintLine:
    def __init__(self, constraint_type, parameters, value):
        self.constraint_type = constraint_type
        self.parameters = parameters
        self.value = value

    def __eq__(self, other):
        if isinstance(other, ConstraintLine):
            return (self.constraint_type == other.constraint_type and
                    self.parameters == other.parameters and
                    self.value == other.value)
        return False

    def __repr__(self):
        return (f"ConstraintLine(constraint_type={self.constraint_type}, "
                f"parameters={self.parameters}, value={self.value})")

