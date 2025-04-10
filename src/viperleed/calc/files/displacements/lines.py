"""Module lines."""

__authors__ = ('Alexander M. Imre (@amimre)',)
__created__ = '2024-10-04'

from collections import namedtuple

from .direction import Direction
from .range import DisplacementsRange
from .targeting import BSTarget

LoopMarkerLine = namedtuple('LoopMarkerLine', ['type'])
SearchHeaderLine = namedtuple('SearchHeaderLine', ['label'])
SectionHeaderLine = namedtuple('SectionHeaderLine', ['section'])


def _get_target(label, which):
    if which is None:
        return BSTarget(label)
    return BSTarget(f'{label} {which}')


class GeoDeltaLine:
    def __init__(self, label, which, direction, start, stop, step,
                 line=None):
        self._line= line
        self.targets = _get_target(label, which)
        self.direction = Direction(direction)
        self.range = DisplacementsRange(start, stop, step)

    def __eq__(self, other):
        if isinstance(other, GeoDeltaLine):
            return (
                self.targets == other.targets
                and self.direction == other.direction
                and self.range == other.range
            )
        return False

    def __repr__(self):
        """Return the string representation of the line."""
        if self._line is None:
            line = f'{self.label} {self.which}'
            if self.direction is not None:
                line += f' {self.direction}'
            line += f' = {self.range.start}'
            if self.range.stop is not None:
                line += f' {self.range.stop}'
            if self.range.step is not None:
                line += f' {self.range.step}'
        else:
            line = self._line
        return line


class VibDeltaLine:
    def __init__(self, label, which, start, stop, step, line=None):
        self._line = line
        self.targets = _get_target(label, which)
        self.range = DisplacementsRange(start, stop, step)

    def __eq__(self, other):
        if isinstance(other, VibDeltaLine):
            return (
                self.targets == other.targets
                and self.range == other.range
            )
        return False

    def __repr__(self):
        """Return the string representation of the line."""
        if self._line is None:
            line = f'{self.label} {self.which}'
            line += f' = {self.range.start}'
            if self.range.stop is not None:
                line += f' {self.range.stop}'
            if self.range.step is not None:
                line += f' {self.range.step}'
        else:
            line = self._line
        return line


class OccDeltaLine:
    def __init__(self, label, which, chem_blocks, line=None):
        self._line = line
        self.targets = _get_target(label, which)
        self.chem_blocks = chem_blocks

    def __eq__(self, other):
        if isinstance(other, OccDeltaLine):
            return (
                self.targets == other.targets
                and self.chem_blocks == other.chem_blocks
            )
        return False

    def __repr__(self):
        """Return the string representation of the line."""
        if self._line is None:
            line = f'{self.label} {self.which}'
            line += f' = {self.chem_blocks}'
            line = self._line
        else:
            line = self._line
        return line


class ConstraintLine:
    def __init__(self, constraint_type, targets, direction, value, line=None):
        self._line = line
        self.constraint_type = constraint_type
        if self.constraint_type == 'geo':
            if direction is not None:
                self.direction = Direction(direction)
            else:
                self.direction = None
        elif self.constraint_type != 'geo' and direction is not None:
            msg = (
                f'A direction is not allowed for {self.constraint_type} '
                f'constraints.'
            )
            raise ValueError(msg)
        self.targets = BSTarget(targets)
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
        """Return the string representation of the line."""
        if self._line is None:
            line = f'{self.constraint_type} {self.targets}'
            if self.direction is not None:
                line += f' {self.direction}'
            line += f' = {self.value}'
        else:
            line = self._line
        return line


class OffsetsLine:
    def __init__(self, offset_type, targets, direction, value, line=None):
        self._line= line
        self.offset_type = offset_type
        self.targets = BSTarget(targets)
        self.value = value
        if self.offset_type == 'geo':
            if direction is not None:
                self.direction = Direction(direction)
            else:
                self.direction = None
        elif self.offset_type != 'geo' and direction is not None:
            msg = (
                f'A direction is not allowed for {self.offset_type} '
                f'offsets.'
            )
            raise ValueError(msg)

    def __eq__(self, other):
        if isinstance(other, ConstraintLine):
            return (
                self.offset_type == other.constraint_type
                and self.parameters == other.parameters
                and self.value == other.value
            )
        return False

    def __repr__(self):
        """Return the string representation of the line."""
        if self._line is None:
            line = f'{self.offset_type} {self.targets}'
            if self.direction is not None:
                line += f' {self.direction}'
            line += f' = {self.value}'
        else:
            line = self._line
        return line
